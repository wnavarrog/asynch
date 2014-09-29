#include <stdio.h>
#include <time.h>
#include <libpq-fe.h>
#include <string.h>
#include <unistd.h>
#include "asynch_interface.h"

int my_rank;
int np;

//!!!! The routines and structs here should probably go into a separate file !!!!

typedef struct CustomParams
{
	unsigned int ID;
	int offset;
} CustomParams;

typedef struct CustomParamsMaps
{
	unsigned int forecast_time;
	unsigned int period;
} CustomParamsMaps;

typedef struct ForecastData
{
	char* model_name;
} ForecastData;

//void CreateHaltFile(char* filename);
//short int CheckFinished(char* filename);
//void PerformTableMaintainance(ConnData* conninfo_stages,ConnData* conninfo_peakflows,ConnData* conninfo_maps,UnivVars* GlobalVars,ForecastData* Forecaster,short int* vac,short unsigned int hr1,unsigned int num_tables);
//void CheckPeakforecastTable(ConnData* conninfo,UnivVars* GlobalVars,ForecastData* Forecaster,unsigned int num_tables);
void CheckPartitionedTable(ConnData* conninfo,UnivVars* GlobalVars,ForecastData* Forecaster,unsigned int num_tables,char* tablename,char* colname);
void UploadPeakflows(asynchsolver* asynch,unsigned int wait_time);

int Output_Linkid(double t,VEC* y_i,VEC* global_params,VEC* params,IVEC* iparams,int state,void* user);
int Output_Timestamp(double t,VEC* y_i,VEC* global_params,VEC* params,IVEC* iparams,int state,void* user);
int Output_ForecastTime(double t,VEC* y_i,VEC* global_params,VEC* params,IVEC* iparams,int state,void* user);
void OutputPeakflow_Forecast_Maps(unsigned int ID,double peak_time,VEC* peak_value,VEC* params,VEC* global_params,double conversion,unsigned int area_idx,void* user,char* buffer);
void Init_Output_User_forecastparams(asynchsolver* asynch);
void Free_Output_User_forecastparams(asynchsolver* asynch);
void Set_Output_User_forecastparams(asynchsolver* asynch,unsigned int offset);

void Init_Output_PeakflowUser_Offset(asynchsolver* asynch);
void Free_Output_PeakflowUser_Offset(asynchsolver* asynch);
void Set_Output_PeakflowUser_Offset(asynchsolver* asynch,unsigned int forecast_time,unsigned int period);

ForecastData* Init_ForecastData(char* fcst_filename);
void Free_ForecastData(ForecastData** Forecaster);


int main(int argc,char* argv[])
{
	//Initialize MPI stuff
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
	MPI_Comm_size(MPI_COMM_WORLD,&np);

	//Parse input
	if(argc < 3)
	{
		if(my_rank == 0)	printf("Command line parameter required:\nA universal variable file (.gbl),\nA forecast file (.fcst)\n");
		MPI_Finalize();
		return 1;
	}

	//Declare variables
	unsigned int i,j,k,current_offset;
	double total_time = 0.0;
	time_t start,start2,stop;
	asynchsolver* asynch;
	PGresult *res;
	MPI_Status status;
	char* query = (char*) malloc(1024*sizeof(char));
	Link* current;

	if(my_rank == 0)
		printf("\nBeginning initialization...\n*****************************\n");
	MPI_Barrier(MPI_COMM_WORLD);
	start = time(NULL);

	//Init asynch object and the river network
	asynch = Asynch_Init(MPI_COMM_WORLD);
	Asynch_Parse_GBL(asynch,argv[1]);
	Asynch_Load_System(asynch);

	//Load Forecast related data
	ForecastData* Forecaster = Init_ForecastData(argv[2]);
	if(!Forecaster)
		MPI_Abort(MPI_COMM_WORLD,1);

	//Setup output for link id, if needed
	int setup_id = Asynch_Check_Output(asynch,"LinkID");
	int setup_timestamp = Asynch_Check_Output(asynch,"Timestamp");
	int setup_forecasttime = Asynch_Check_Output(asynch,"ForecastTime");
	if( (setup_id || setup_timestamp || setup_forecasttime) != 0 )
	{
		if(my_rank == 0)	printf("[%i]: Forecaster needs LinkID (%i), Timestamp (%i), Forecasttime (%i).\n",my_rank,setup_id,setup_timestamp,setup_forecasttime);
		MPI_Abort(MPI_COMM_WORLD,1);
	}
	Init_Output_User_forecastparams(asynch);
	Asynch_Set_Output(asynch,"LinkID",ASYNCH_INT,(void (*)(double,VEC*,VEC*,VEC*,IVEC*,int,void*)) &Output_Linkid,NULL,0);
	Asynch_Set_Output(asynch,"Timestamp",ASYNCH_INT,(void (*)(double,VEC*,VEC*,VEC*,IVEC*,int,void*)) &Output_Timestamp,NULL,0);
	Asynch_Set_Output(asynch,"ForecastTime",ASYNCH_INT,(void (*)(double,VEC*,VEC*,VEC*,IVEC*,int,void*)) &Output_ForecastTime,NULL,0);

	//Setup the peakflow information for maps
	int setup_peakflow_maps = Asynch_Check_Peakflow_Output(asynch,"Forecast_Maps");
	if(setup_peakflow_maps != 0)
	{
		if(my_rank == 0)	printf("[%i]: Forecaster with maps needs Forecast_Maps (%i) for peakflows.\n",my_rank,setup_peakflow_maps);
		MPI_Abort(MPI_COMM_WORLD,1);
	}
	Init_Output_PeakflowUser_Offset(asynch);
	Asynch_Set_Peakflow_Output(asynch,"Forecast_Maps",(void (*)(unsigned int,double,VEC*,VEC*,VEC*,double,unsigned int,void*,char*)) &OutputPeakflow_Forecast_Maps);

	//Get some values about the river system
	unsigned int N = Asynch_Get_Number_Links(asynch);
	unsigned int my_N = Asynch_Get_Local_Number_Links(asynch);
	char dump_filename[asynch->GlobalVars->string_size];

	//Create halt file
	//CreateHaltFile(asynch->GlobalVars->halt_filename);

	//Find the index of the forcing to use for forecasting
	unsigned int forecast_idx = asynch->GlobalVars->num_forcings;
	for(i=0;i<asynch->GlobalVars->num_forcings;i++)
	{
		if(asynch->forcings[i]->flag == 5)
		{
			if(forecast_idx == asynch->GlobalVars->num_forcings)
				forecast_idx = i;
			else
				if(my_rank == 0)	printf("[%i]: Warning: Multiple forecasting forcings may not be setup yet...\n",my_rank);
		}
	}
	if(forecast_idx == asynch->GlobalVars->num_forcings)
	{
		printf("[%i]: Error: No forecasting forcing set.\n",my_rank);
		MPI_Abort(MPI_COMM_WORLD,1);
	}

	//Reserve space for backups
	VEC** backup = (VEC**) malloc(N*sizeof(VEC*));
	for(i=0;i<N;i++)
	{
		if(asynch->assignments[i] == my_rank || asynch->getting[i] == 1)
			backup[i] = v_get(asynch->GlobalVars->dim);
		else
			backup[i] = NULL;
	}

	if(my_rank == 0)
	{
		printf("\nModel type is %u.\nGlobal parameters are:\n",asynch->GlobalVars->type);
		Print_Vector(asynch->GlobalVars->global_params);
		printf("\n");
	}

	MPI_Barrier(MPI_COMM_WORLD);
	stop = time(NULL);
	total_time = difftime(stop,start);
	if(my_rank == 0)	printf("Finished initializations. Total time: %f\n",total_time);
	MPI_Barrier(MPI_COMM_WORLD);
	sleep(1);

	//Make sure everyone is good before getting down to it...
	printf("Process %i (%i total) is good to go with %i links.\n",my_rank,np,my_N);
	MPI_Barrier(MPI_COMM_WORLD);
	start = time(NULL);

	//Make the initial solve
	Asynch_Advance(asynch,0);

	//Stop the clock
	MPI_Barrier(MPI_COMM_WORLD);
	stop = time(NULL);
	total_time += difftime(stop,start);

	//Output some data
	if(my_rank == 0)
	{
		printf("%i: The answer at ID %i at time %.12f is\n",my_rank,asynch->sys[asynch->my_sys[0]]->ID,asynch->sys[asynch->my_sys[0]]->last_t);
		Print_Vector(asynch->sys[asynch->my_sys[0]]->list->tail->y_approx);
		printf("Total time for calculations: %f\n",difftime(stop,start));
	}

	//Begin persistent calculations
	if(my_rank == 0)
		printf("\n\n===================================\nBeginning persistent calculations\n===================================\n");
	fflush(stdout);
	MPI_Barrier(MPI_COMM_WORLD);

	//Make some initializations and checks
	unsigned int history_time = 5*24*60*60;	//Amount of history to store for hydrographs and peakflows
	short unsigned int hr1 = 0;	//Hour of the day to perform maintainance on database
	short unsigned int hr2 = 12;	//Hour of the day to perform maintainance on database
	unsigned int wait_time = 120;	//Time to sleep if no rainfall data is available
	double forecast_time = 10.0*24*60;	//Time (mins) in future to make forecasts
	unsigned int num_tables = 10;
	unsigned int db_retry_time = 5;	//Time (secs) to wait if a database error occurs
	unsigned int num_future_peakflow_times = 9;
	double future_peakflow_times[] = {60.0, 180.0, 360.0, 720.0, 1440.0, 2880.0, 4320.0, 5760.0, 7200.0};
	//unsigned int num_rainsteps = 3;	//Number of rainfall intensities to use for the next forecast
	unsigned int num_rainsteps = asynch->forcings[forecast_idx]->num_rainsteps;	//Number of rainfall intensities to use for the next forecast
	//if(my_rank == 0 && asynch->GlobalVars->increment < num_rainsteps + 3)
	if(my_rank == 0 && asynch->forcings[forecast_idx]->increment < num_rainsteps + 3)
		printf("Warning: Increment for rain should probably be %u.\n",num_rainsteps + 3);
	asynch->forcings[forecast_idx]->increment = num_rainsteps;	//!!!! Not necessary, but makes me feel better. The solvers should really not do the last step where they download nothing. !!!!
	double db_stepsize = asynch->forcings[forecast_idx]->file_time,t;

	unsigned int nextraintime,repeat_for_errors,nextforcingtime;
	//short int halt = 0;
	int isnull;
	int message_buffer[1 + asynch->GlobalVars->num_forcings];
	short int vac = 0;	//0 if no vacuum has occured, 1 if vacuum has occured (during a specific hour)
	unsigned int last_file = asynch->forcings[forecast_idx]->last_file;
	unsigned int first_file = asynch->forcings[forecast_idx]->first_file;
	k = 0;
	for(i=0;i<N;i++)
		if(backup[i] != NULL)	v_copy(asynch->sys[i]->list->tail->y_approx,backup[i]);

	double simulation_time_with_data = 0.0;
	//for(i=0;i<asynch->GlobalVars->num_forcings;i++)
	//	if(asynch->forcings[i]->flag == 5)	simulation_time_with_data = max(simulation_time_with_data,asynch->forcings[i]->file_time * asynch->forcings[i]->num_rainsteps);
	simulation_time_with_data = max(simulation_time_with_data,asynch->forcings[forecast_idx]->file_time * asynch->forcings[forecast_idx]->num_rainsteps);

	//Setup temp files
	Set_Output_User_forecastparams(asynch,first_file);
	//Set_Output_PeakflowUser_Offset(asynch,first_file);
	Asynch_Set_Total_Simulation_Time(asynch,forecast_time);
	Asynch_Prepare_Temp_Files(asynch);

/*
	//Make some initializations to the database
	if(my_rank == 0)
	{
		printf("Making initializations to tables.\n");
		start = time(NULL);

		//Connect to hydrograph database
		ConnectPGDB(asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT]);

		//Make sure the hydrographs table exists
		sprintf(query,"SELECT 1 FROM pg_class WHERE relname='%s';",asynch->GlobalVars->hydro_table);
		res = PQexec(asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT]->conn,query);
		if(!PQntuples(res))
		{
			PQclear(res);
			sprintf(query,"CREATE TABLE %s(link_id int,time int,ratio real,discharge real); ALTER TABLE %s SET (autovacuum_enabled = false, toast.autovacuum_enabled = false);",asynch->GlobalVars->hydro_table,asynch->GlobalVars->hydro_table);
			res = PQexec(asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT]->conn,query);
			CheckResError(res,"creating hydrographs table");
		}
		else
		{
			PQclear(res);
			sprintf(query,"TRUNCATE %s;",asynch->GlobalVars->hydro_table);
			res = PQexec(asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT]->conn,query);
			CheckResError(res,"truncating hydrographs table");
		}
		PQclear(res);

		//Clear the future hydrographs in archive
		sprintf(query,"DELETE FROM master_archive_hydroforecast_%s WHERE forecast_time >= %u;",Forecaster->model_name,first_file);
		res = PQexec(asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT]->conn,query);
		CheckResError(res,"deleting future hydroforecasts");
		PQclear(res);

		//Make sure the hydroforecast tables are set correctly
		CheckPartitionedTable(asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT],asynch->GlobalVars,Forecaster,num_tables,"archive_hydroforecast","forecast_time");

		//Disconnect from hydrograph database, connect to peakflow database
		DisconnectPGDB(asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT]);
		ConnectPGDB(asynch->db_connections[ASYNCH_DB_LOC_PEAK_OUTPUT]);

		//Clear all future peakflows
		sprintf(query,"DELETE FROM master_archive_peakflows_%s WHERE forecast_time >= %u;",Forecaster->model_name,first_file);
		res = PQexec(asynch->db_connections[ASYNCH_DB_LOC_PEAK_OUTPUT]->conn,query);
		CheckResError(res,"deleting future peakforecast");
		PQclear(res);

		//Make sure the peakflow tables are set correctly
		CheckPartitionedTable(asynch->db_connections[ASYNCH_DB_LOC_PEAK_OUTPUT],asynch->GlobalVars,Forecaster,num_tables,"archive_peakflows","forecast_time");

		//Disconnect from peakflow database, connect to snapshot database
		DisconnectPGDB(asynch->db_connections[ASYNCH_DB_LOC_PEAK_OUTPUT]);
		ConnectPGDB(asynch->db_connections[ASYNCH_DB_LOC_SNAPSHOT_OUTPUT]);

		//Clear all future maps
		sprintf(query,"DELETE FROM master_archive_maps_%s WHERE forecast_time > %u;",Forecaster->model_name,first_file);
		res = PQexec(asynch->db_connections[ASYNCH_DB_LOC_SNAPSHOT_OUTPUT]->conn,query);
		CheckResError(res,"deleting future snapshots");
		PQclear(res);

		//Make sure the map tables are set correctly
		CheckPartitionedTable(asynch->db_connections[ASYNCH_DB_LOC_SNAPSHOT_OUTPUT],asynch->GlobalVars,Forecaster,num_tables,"archive_maps","forecast_time");

		//Disconnect from snapshot database
		DisconnectPGDB(asynch->db_connections[ASYNCH_DB_LOC_SNAPSHOT_OUTPUT]);

		stop = time(NULL);
		printf("Total time to initialize tables: %.2f.\n",difftime(stop,start));
	}
*/
	MPI_Barrier(MPI_COMM_WORLD);

	//Start the main loop
	//while(!halt)
	int passes = 48;
	for(k=0;k<passes;k++)
	{
		if(my_rank == 0)	printf("\n\nPass %u\n",k);

		//Clear buffers
		Flush_TransData(asynch->my_data);

		//Make some initializations
		//asynch->forcings[forecast_idx]->raindb_start_time = last_file;								//!!!! This all assumes one forcing from db !!!!
		first_file = last_file;
		last_file = last_file + (unsigned int) asynch->forcings[forecast_idx]->file_time * 60 * num_rainsteps;
		nextforcingtime = first_file + 60 * (unsigned int) rint(asynch->forcings[forecast_idx]->file_time) * (num_rainsteps-1);	//This is the actual timestamp of the last needed forcing data. This will be downloaded (unlike last_file)
		//nextforcingtime = first_file + 60 * (unsigned int) rint(asynch->forcings[forecast_idx]->file_time) * num_rainsteps;

		//Reset each link
		Asynch_Set_System_State(asynch,0.0,backup);
		Set_Output_User_forecastparams(asynch,first_file);
		Set_Output_PeakflowUser_Offset(asynch,first_file,first_file);
		Asynch_Write_Current_Step(asynch);
		Asynch_Set_Forcing_State(asynch,forecast_idx,0.0,first_file,last_file);

		//Check if a vacuum should be done
		//This will happen at hr1
		//if(my_rank == 0)	PerformTableMaintainance(asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT],asynch->db_connections[ASYNCH_DB_LOC_PEAK_OUTPUT],asynch->db_connections[ASYNCH_DB_LOC_SNAPSHOT_OUTPUT],asynch->GlobalVars,Forecaster,&vac,hr1,num_tables);

		//Make sure all buffer flushing is done
		MPI_Barrier(MPI_COMM_WORLD);

/*
		//Dump data for debugging and recovery
		if(k % 96 == 0)
		{
			//sprintf(dump_filename,"%u",last_file);
			sprintf(dump_filename,"/%u",first_file);
			Asynch_Take_System_Snapshot(asynch,dump_filename);
		}
*/

/*
		//Find the next time where rainfall occurs
		do
		{
			if(my_rank == 0)
			{
				//Init message_buffer
				message_buffer[0] = 1;

				//Check if connection to SQL database is still good
				for(i=0;i<asynch->GlobalVars->num_forcings;i++)
				{
					message_buffer[i+1] = 0x7FFFFFFF;	//Largest 32 bit int	!!!! I don't think this is necessary !!!!

					if(asynch->forcings[i]->flag == 5)
					{
						ConnectPGDB(asynch->db_connections[ASYNCH_DB_LOC_FORCING_START + i]);

						//Find the next rainfall time
						time(&start);
						sprintf(query,"SELECT min(unix_time) FROM rain_maps5_index WHERE unix_time >= %u AND link_count > -1;",nextforcingtime);	//!!!! Should vary with forcing (i) !!!!
						res = PQexec(asynch->db_connections[ASYNCH_DB_LOC_FORCING_START + i]->conn,query);
						CheckResError(res,"checking for new rainfall data");
						time(&stop);
						printf("Total time to check for new rainfall data: %f.\n",difftime(stop,start));

						isnull = PQgetisnull(res,0,0);
						if(!isnull)
						{
							message_buffer[0] = 0;
							int value = atoi(PQgetvalue(res,0,0));
							message_buffer[i+1] = (value < message_buffer[i+1]) ? value : message_buffer[i+1];
							PQclear(res);
							DisconnectPGDB(asynch->db_connections[ASYNCH_DB_LOC_FORCING_START + i]);
						}
						else
						{
							message_buffer[0] = 1;
							message_buffer[i+1] = 0;
							PQclear(res);
							DisconnectPGDB(asynch->db_connections[ASYNCH_DB_LOC_FORCING_START + i]);
							break;
						}
					}
				}
			}
			MPI_Bcast(message_buffer,1+asynch->GlobalVars->num_forcings,MPI_INT,0,MPI_COMM_WORLD);
			isnull = message_buffer[0];

			if(!isnull)	//!!!! Seems I don't need nextraintime. Do I need to even send more than just isnull? !!!!
				nextraintime = message_buffer[1];
			else
			{
				if(my_rank == 0)
				{
					printf("No rainfall values returned from SQL database for forcing %u. %u %u\n",i,last_file,isnull);
					//PerformTableMaintainance(asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT],asynch->db_connections[ASYNCH_DB_LOC_PEAK_OUTPUT],asynch->db_connections[ASYNCH_DB_LOC_SNAPSHOT_OUTPUT],asynch->GlobalVars,Forecaster,&vac,hr1,num_tables);
				}

				halt = CheckFinished(asynch->GlobalVars->halt_filename);
				if(halt)
				{
					//sprintf(dump_filename,"/%u",first_file);
					//Asynch_Take_System_Snapshot(asynch,dump_filename);
				}
				else
				{
					fflush(stdout);
					sleep(wait_time);
				}
			}
		//} while(isnull && !halt);
		} while(isnull);
*/
		//if(halt)	break;

		//Read in next set of rainfall data

		//Initialize some data for the first phase of calculations
		//GlobalVars->maxtime = GlobalVars->file_time * num_rainsteps;
		Asynch_Set_Total_Simulation_Time(asynch,simulation_time_with_data);		// !!!! This may not work for multiple forcings for forecasting. How do you handle different time resolutions? !!!!
		//conninfo->time_offset = first_file;
		current_offset = first_file;
		Set_Output_User_forecastparams(asynch,current_offset);
		Set_Output_PeakflowUser_Offset(asynch,current_offset,current_offset);

		MPI_Barrier(MPI_COMM_WORLD);
		time(&start);
if(my_rank == 0)
printf("first: %u last: %u\n",first_file,last_file);

		Asynch_Advance(asynch,1);

		MPI_Barrier(MPI_COMM_WORLD);
		if(my_rank == 0)
		{
			time(&stop);
			printf("Time for first phase calculations: %.2f\n",difftime(stop,start));
		}

		//Flush communication buffers
		Flush_TransData(asynch->my_data);

		//Reset the links (mostly) and make a backup for the second phase
		//!!!! Need routine for this !!!!
		for(i=0;i<N;i++)	//Set time to 0.0
		{
			current = asynch->sys[i];
			if(current->list != NULL)
			{
				while(current->current_iterations > 1)
				{
					Remove_Head_Node(current->list);
					(current->current_iterations)--;
				}

				current->steps_on_diff_proc = 1;
				current->iters_removed = 0;
				current->rejected = 0;
				if(current->numparents == 0)	current->ready = 1;
				else				current->ready = 0;

				v_copy(current->list->head->y_approx,backup[i]);
			}
		}

		//Upload a snapshot to the database
		sprintf(dump_filename,"%u",first_file);
		Asynch_Take_System_Snapshot(asynch,dump_filename);

		//Make second phase calculations. Peakflow data will be uploaded several times.
		MPI_Barrier(MPI_COMM_WORLD);
		time(&start);

		Asynch_Deactivate_Forcing(asynch,forecast_idx);

		for(i=0;i<num_future_peakflow_times;i++)
		{
			//t = future_peakflow_times[i] + db_stepsize*num_rainsteps;
			t = asynch->sys[asynch->my_sys[0]]->last_t;
			//Asynch_Set_Total_Simulation_Time(asynch,t);
			Asynch_Set_Total_Simulation_Time(asynch,future_peakflow_times[i] + db_stepsize*num_rainsteps);
			Asynch_Reset_Peakflow_Data(asynch);
			Set_Output_PeakflowUser_Offset(asynch,current_offset,current_offset + (unsigned int) (60.0*t+0.1));
			Asynch_Advance(asynch,1);
			UploadPeakflows(asynch,db_retry_time);
		}

		Asynch_Reset_Peakflow_Data(asynch);
		Asynch_Set_Total_Simulation_Time(asynch,forecast_time);
		Asynch_Advance(asynch,1);

		Asynch_Activate_Forcing(asynch,forecast_idx);
/*
		Asynch_Set_Total_Simulation_Time(asynch,forecast_time);
		Asynch_Deactivate_Forcing(asynch,forecast_idx);
		Asynch_Advance(asynch,1);
		Asynch_Activate_Forcing(asynch,forecast_idx);
*/
		MPI_Barrier(MPI_COMM_WORLD);
		if(my_rank == 0)
		{
			time(&stop);
			printf("Time for second phase calculations: %.2f\n",difftime(stop,start));
		}

		//Output some data
		if(my_rank == 0)
		{
			printf("[%i]: The answer at ID %i at time %.12f is\n",my_rank,asynch->sys[asynch->my_sys[0]]->ID,asynch->sys[asynch->my_sys[0]]->last_t);
			Print_Vector(asynch->sys[asynch->my_sys[0]]->list->tail->y_approx);
		}

/*
		//Upload the peak data to the database **********************************************************************************************
		MPI_Barrier(MPI_COMM_WORLD);
		start = time(NULL);

		UploadPeakflows(asynch,db_retry_time);

		MPI_Barrier(MPI_COMM_WORLD);
		if(my_rank == 0)
		{
			stop = time(NULL);
			printf("[%i]: Total time to transfer peak flow data: %.2f\n",my_rank,difftime(stop,start));
		}
*/

		//Upload the hydrographs to the database ********************************************************************************************
		MPI_Barrier(MPI_COMM_WORLD);
		start = time(NULL);
/*
		//Adjust the table hydrographs
		if(my_rank == 0)
		{
			//Make sure database connection is still good
			ConnectPGDB(asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT]);

			sprintf(query,"TRUNCATE %s;",asynch->GlobalVars->hydro_table);
			res = PQexec(asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT]->conn,query);
			CheckResError(res,"deleting hydrographs");
			PQclear(res);

			DisconnectPGDB(asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT]);
		}
		MPI_Barrier(MPI_COMM_WORLD);
*/
		repeat_for_errors = Asynch_Create_Output(asynch);
		while(repeat_for_errors > 0)
		{
			if(my_rank == 0)	printf("[%i]: Attempting resend of hydrographs data.\n",my_rank);
			sleep(5);
			repeat_for_errors = Asynch_Create_Output(asynch);
		}
/*
		//Call functions *********************************************************************************************************************
		if(my_rank == 0)
		{
			//Connect to database
			ConnectPGDB(asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT]);

			//Stage
			repeat_for_errors = 1;
			while(repeat_for_errors)
			{
				repeat_for_errors = 0;
				sprintf(query,"SELECT get_stages_%s();",Forecaster->model_name);
				res = PQexec(asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT]->conn,query);
				repeat_for_errors = repeat_for_errors || CheckResError(res,"calling stage function");
				PQclear(res);
				if(repeat_for_errors)
				{
					printf("[%i]: Attempting to call stage function again...\n",my_rank);
					sleep(5);
					CheckConnConnection(asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT]);
				}
			}

			//Warnings
			repeat_for_errors = 1;
			while(repeat_for_errors)
			{
				repeat_for_errors = 0;
				sprintf(query,"SELECT update_warnings_%s();",Forecaster->model_name);
				res = PQexec(asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT]->conn,query);
				repeat_for_errors = repeat_for_errors || CheckResError(res,"calling warnings function");
				PQclear(res);
				if(repeat_for_errors)
				{
					printf("[%i]: Attempting to call warning function again...\n",my_rank);
					sleep(5);
					CheckConnConnection(asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT]);
				}
			}

			//Stage archive
			repeat_for_errors = 1;
			while(repeat_for_errors)
			{
				repeat_for_errors = 0;
				sprintf(query,"ALTER TABLE master_archive_hydroforecast_%s ALTER COLUMN forecast_time SET DEFAULT %u;",Forecaster->model_name,current_offset);
				res = PQexec(asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT]->conn,query);
				repeat_for_errors = repeat_for_errors || CheckResError(res,"setting default value");
				PQclear(res);

				sprintf(query,"SELECT copy_to_archive_hydroforecast_%s();",Forecaster->model_name);
				res = PQexec(asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT]->conn,query);
				repeat_for_errors = repeat_for_errors || CheckResError(res,"calling stage archive function");
				PQclear(res);

				sprintf(query,"ALTER TABLE master_archive_hydroforecast_%s ALTER COLUMN forecast_time DROP DEFAULT;",Forecaster->model_name);
				res = PQexec(asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT]->conn,query);
				repeat_for_errors = repeat_for_errors || CheckResError(res,"dropping default value");
				PQclear(res);

				if(repeat_for_errors)
				{
					printf("[%i]: Attempting to call stage archive function again...\n",my_rank);
					sleep(5);
					CheckConnConnection(asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT]);
				}
			}

			//Disconnect
			DisconnectPGDB(asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT]);
		}
*/
		if(my_rank == 0)
		{
			time(&stop);
			printf("[%i]: Total time to transfer hydrograph data: %.2f\n",my_rank,difftime(stop,start));
		}
		fflush(stdout);
		MPI_Barrier(MPI_COMM_WORLD);

		//Check if program has received a terminate signal **********************************************************************************
		//k++;
		//halt = CheckFinished(asynch->GlobalVars->halt_filename);
	}

	//Clean up **********************************************************************************************************************************
	free(query);
	for(i=0;i<N;i++)	v_free(backup[i]);
	free(backup);
	Free_ForecastData(&Forecaster);
	Asynch_Delete_Temporary_Files(asynch);
	Free_Output_PeakflowUser_Offset(asynch);
	Free_Output_User_forecastparams(asynch);
	Asynch_Free(asynch);
	return 0;
}

/*
//Checks if the time is right to perform maintainance on the database.
void PerformTableMaintainance(ConnData* conninfo_stages,ConnData* conninfo_peakflows,ConnData* conninfo_maps,UnivVars* GlobalVars,ForecastData* Forecaster,short int* vac,short unsigned int hr1,unsigned int num_tables)
{
	int j;
	time_t start,stop;
	PGresult* res;
	char query[GlobalVars->string_size];
	struct tm* timeinfo;
	time(&start);
	timeinfo = localtime(&start);

	if(timeinfo->tm_hour == hr1 && *vac == 0)	//Main maintainance time
	{
		printf("[%i]: Performing maintainance. Current time is %s",my_rank,asctime(timeinfo));

		//Adjust partitioned hydroforecast tables
		ConnectPGDB(conninfo_stages);
		sprintf(query,"DROP TABLE archive_hydroforecast_%s_%u;",Forecaster->model_name,num_tables-1);
		res = PQexec(conninfo_stages->conn,query);
		CheckResError(res,"dropping end archive table");
		PQclear(res);

		for(j=num_tables-2;j>=0;j--)
		{
			sprintf(query,"ALTER TABLE archive_hydroforecast_%s_%i RENAME TO archive_hydroforecast_%s_%i;",Forecaster->model_name,j,Forecaster->model_name,j+1);
			res = PQexec(conninfo_stages->conn,query);
			CheckResError(res,"renaming table");
			PQclear(res);
		}

		sprintf(query,"CREATE TABLE archive_hydroforecast_%s_0 ( ) INHERITS (master_archive_hydroforecast_%s);",Forecaster->model_name,Forecaster->model_name);
		res = PQexec(conninfo_stages->conn,query);
		CheckResError(res,"creating archive_hydroforecast_0");
		PQclear(res);
		DisconnectPGDB(conninfo_stages);

		//Adjust peakflow table
		ConnectPGDB(conninfo_peakflows);
		sprintf(query,"DROP TABLE archive_peakflows_%s_%u;",Forecaster->model_name,num_tables-1);
		res = PQexec(conninfo_peakflows->conn,query);
		CheckResError(res,"dropping end peakflows archive table");
		PQclear(res);

		for(j=num_tables-2;j>=0;j--)
		{
			sprintf(query,"ALTER TABLE archive_peakflows_%s_%i RENAME TO archive_peakflows_%s_%i;",Forecaster->model_name,j,Forecaster->model_name,j+1);
			res = PQexec(conninfo_peakflows->conn,query);
			CheckResError(res,"renaming peakflows table");
			PQclear(res);
		}

		sprintf(query,"CREATE TABLE archive_peakflows_%s_0 ( ) INHERITS (master_archive_peakflows_%s);",Forecaster->model_name,Forecaster->model_name);
		res = PQexec(conninfo_peakflows->conn,query);
		CheckResError(res,"creating archive_peakflows_0");
		PQclear(res);
		DisconnectPGDB(conninfo_peakflows);

		//Adjust maps table
		ConnectPGDB(conninfo_maps);
		sprintf(query,"DROP TABLE archive_maps_%s_%u;",Forecaster->model_name,num_tables-1);
		res = PQexec(conninfo_maps->conn,query);
		CheckResError(res,"dropping end maps archive table");
		PQclear(res);

		for(j=num_tables-2;j>=0;j--)
		{
			sprintf(query,"ALTER TABLE archive_maps_%s_%i RENAME TO archive_maps_%s_%i;",Forecaster->model_name,j,Forecaster->model_name,j+1);
			res = PQexec(conninfo_maps->conn,query);
			CheckResError(res,"renaming maps table");
			PQclear(res);
		}

		sprintf(query,"CREATE TABLE archive_maps_%s_0 ( ) INHERITS (master_archive_maps_%s);",Forecaster->model_name,Forecaster->model_name);
		res = PQexec(conninfo_maps->conn,query);
		CheckResError(res,"creating archive_maps_0");
		PQclear(res);
		DisconnectPGDB(conninfo_maps);

		//Set flag, print time, disconnect
		*vac = 1;
		time(&stop);
		printf("[%i]: Database cleanup complete. Total time %.2f.\n\n",my_rank,difftime(stop,start));
	}
	else	if(timeinfo->tm_hour != hr1)	*vac = 0;
}
*/

//Checks that the timestamps of a partitioned table match up correctly with the trigger.
//If not, the child tables are adjusted.
//void CheckPeakforecastTable(ConnData* conninfo,UnivVars* GlobalVars,ForecastData* Forecaster,unsigned int num_tables)
void CheckPartitionedTable(ConnData* conninfo,UnivVars* GlobalVars,ForecastData* Forecaster,unsigned int num_tables,char* tablename,char* colname)
{
	unsigned int i,current_time,table_time,table_index,correct_table_index,diff_table_index;
	int diff_time,j,last_table_index;
	PGresult *res;
	char* query = conninfo->query;

	//Grab the current time from the database
	sprintf(query,"SELECT EXTRACT('epoch' FROM current_date AT time zone 'UTC');");
	res = PQexec(conninfo->conn,query);
	current_time = (unsigned int) rint(atof(PQgetvalue(res,0,0)));
	PQclear(res);

	//Find the first table with something in it
	table_index = num_tables;
	for(i=0;i<num_tables;i++)
	{
		//sprintf(query,"SELECT forecast_time FROM archive_hydroforecast_%s_%u LIMIT 1;",Forecaster->model_name,i);
		sprintf(query,"SELECT %s FROM %s_%s_%u LIMIT 1;",colname,tablename,Forecaster->model_name,i);
		res = PQexec(conninfo->conn,query);
		CheckResError(res,"checking table contents");
		if(PQntuples(res))
		{
			table_time = (unsigned int) atoi(PQgetvalue(res,0,0));
			diff_time = table_time - current_time;
			table_index = i;
			PQclear(res);
			break;
		}
		PQclear(res);
	}

	//If no data exists in any table, then all is well.
	if(table_index == num_tables)	return;

	//If table_time is in the correct table, return.
	if(-86400*(int)table_index < (int) diff_time && (int) diff_time <= -86400*((int)table_index-1))	return;

	//Otherwise, find the correct table
	correct_table_index = num_tables;
	for(i=table_index+1;i<num_tables;i++)
	{
		if(-86400*(int)i < (int)diff_time && (int)diff_time <= -86400*((int)i-1))
		{
			correct_table_index = i;
			break;
		}
	}
	diff_table_index = correct_table_index - table_index;

	//Trash the tables at the end
	last_table_index = num_tables-diff_table_index - 1;
	for(j=num_tables-1;j>last_table_index;j--)
	{
		//sprintf(query,"DROP TABLE archive_stageforecast_%s_%u;",Forecaster->model_name,j);
		sprintf(query,"DROP TABLE %s_%s_%u;",tablename,Forecaster->model_name,j);
		res = PQexec(conninfo->conn,query);
		CheckResError(res,"dropping table");
		PQclear(res);
	}

	//Move all the tables
	for(j=last_table_index;j>=0;j--)
	{
		//sprintf(query,"ALTER TABLE archive_stageforecast_%s_%i RENAME TO archive_stageforecast_%s_%i;",Forecaster->model_name,j,Forecaster->model_name,j+diff_table_index);
		sprintf(query,"ALTER TABLE %s_%s_%i RENAME TO %s_%s_%i;",tablename,Forecaster->model_name,j,tablename,Forecaster->model_name,j+diff_table_index);
		res = PQexec(conninfo->conn,query);
		CheckResError(res,"renaming table");
		PQclear(res);
	}

	//Create new tables
	for(i=0;i<diff_table_index;i++)
	{
		//sprintf(query,"CREATE TABLE archive_stageforecast_%s_%u ( ) INHERITS (master_archive_stageforecast_%s);",Forecaster->model_name,i,Forecaster->model_name);
		sprintf(query,"CREATE TABLE %s_%s_%u ( ) INHERITS (%s_%s);",tablename,Forecaster->model_name,i,tablename,Forecaster->model_name);
		res = PQexec(conninfo->conn,query);
		CheckResError(res,"creating table");
		PQclear(res);
	}
}

/*
//Creates the halt file and sets the value to 0
void CreateHaltFile(char* filename)
{
	FILE* outputfile;

	if(my_rank == 0)
	{
		outputfile = fopen(filename,"w");
		if(!outputfile)	printf("Warning: Could not create halt file %s.\n",filename);
		fprintf(outputfile,"0");
		fclose(outputfile);
	}

	MPI_Barrier(MPI_COMM_WORLD);
}
*/

/*
//Returns 0 if the program should continue, 1 if the program should terminate.
short int CheckFinished(char* filename)
{
	FILE* inputfile;
	short int halt;
	time_t now;
	struct tm* timeinfo;

	if(my_rank == 0)
	{
		//Check terminate file
		inputfile = fopen(filename,"r");
		if(!inputfile)	halt = 0;
		else		fscanf(inputfile,"%hu",&halt);

		//Check if signal received
		if(halt)
		{
			time(&now);
			timeinfo = localtime(&now);
			printf("\nReceived halt signal on %s",asctime(timeinfo));
		}

		fclose(inputfile);
	}

	//Notify all processes of how to proceed
	MPI_Bcast(&halt,1,MPI_SHORT,0,MPI_COMM_WORLD);

	return halt;
}
*/

//Calls the function to create peakflows. The function is called repeatedly until the data is sent.
void UploadPeakflows(asynchsolver* asynch,unsigned int wait_time)
{
	int repeat_for_errors;
	MPI_Barrier(MPI_COMM_WORLD);

	repeat_for_errors = Asynch_Create_Peakflows_Output(asynch);
	while(repeat_for_errors > 0)
	{
		if(my_rank == 0)	printf("[%i]: Attempting resend of peakflow data.\n",my_rank);
		sleep(wait_time);
		repeat_for_errors = Asynch_Create_Peakflows_Output(asynch);
	}
}

//Output functions ****************************************************************************
int Output_Linkid(double t,VEC* y_i,VEC* global_params,VEC* params,IVEC* iparams,int state,void* user)
{
	CustomParams* forecastparams = (CustomParams*) user;
	return forecastparams->ID;
}

int Output_Timestamp(double t,VEC* y_i,VEC* global_params,VEC* params,IVEC* iparams,int state,void* user)
{
	CustomParams* forecastparams = (CustomParams*) user;
	return (int)(round(t * 60.0 + forecastparams->offset) + 0.1);
}

int Output_ForecastTime(double t,VEC* y_i,VEC* global_params,VEC* params,IVEC* iparams,int state,void* user)
{
	CustomParams* forecastparams = (CustomParams*) user;
	return forecastparams->offset;
}

/*
double Output_Ratio(double t,VEC* y_i,VEC* global_params,VEC* params,IVEC* iparams,int state,void* user)
{
	CustomParams* forecastparams = (CustomParams*) user;
	return y_i->ve[0] / forecastparams->Q_TM;
}
*/
void OutputPeakflow_Forecast_Maps(unsigned int ID,double peak_time,VEC* peak_value,VEC* params,VEC* global_params,double conversion,unsigned int area_idx,void* user,char* buffer)
{
	CustomParamsMaps* forecastparams = (CustomParamsMaps*) user;
	sprintf(buffer,"%u,%u,%.6e,%u,%u\n",ID,forecastparams->forecast_time + (unsigned int)(peak_time*60 + .1),peak_value->ve[0],forecastparams->forecast_time,forecastparams->period);
}


//Custom parameters for forecasting ***********************************************************
void Init_Output_User_forecastparams(asynchsolver* asynch)
{
	unsigned int i,my_N = asynch->my_N,*my_sys = asynch->my_sys;
	Link** sys = asynch->sys;

	for(i=0;i<my_N;i++)
		sys[my_sys[i]]->output_user = malloc(sizeof(CustomParams));
}

void Free_Output_User_forecastparams(asynchsolver* asynch)
{
	unsigned int i,my_N = asynch->my_N,*my_sys = asynch->my_sys;
	Link** sys = asynch->sys;

	for(i=0;i<my_N;i++)
	{
		free(sys[my_sys[i]]->output_user);
		sys[my_sys[i]]->output_user = NULL;
	}
}

void Set_Output_User_forecastparams(asynchsolver* asynch,unsigned int offset)
{
	unsigned int i,my_N = asynch->my_N,*my_sys = asynch->my_sys;
	Link** sys = asynch->sys;
	CustomParams* forecastparams;

	for(i=0;i<my_N;i++)
	{
		forecastparams = (CustomParams*) sys[my_sys[i]]->output_user;
		forecastparams->ID = sys[my_sys[i]]->ID;
		forecastparams->offset = offset;
	}
}

void Init_Output_PeakflowUser_Offset(asynchsolver* asynch)
{
	unsigned int i,my_N = asynch->my_N,*my_sys = asynch->my_sys;
	Link** sys = asynch->sys;

	for(i=0;i<my_N;i++)
		sys[my_sys[i]]->peakoutput_user = malloc(sizeof(CustomParamsMaps));
}

void Free_Output_PeakflowUser_Offset(asynchsolver* asynch)
{
	unsigned int i,my_N = asynch->my_N,*my_sys = asynch->my_sys;
	Link** sys = asynch->sys;

	for(i=0;i<my_N;i++)
	{
		free(sys[my_sys[i]]->peakoutput_user);
		sys[my_sys[i]]->peakoutput_user = NULL;
	}
}

void Set_Output_PeakflowUser_Offset(asynchsolver* asynch,unsigned int forecast_time,unsigned int period)
{
	unsigned int i,my_N = asynch->my_N,*my_sys = asynch->my_sys;
	Link** sys = asynch->sys;
	CustomParamsMaps* peak_params;

	for(i=0;i<my_N;i++)
	{
		//*(unsigned int*)(sys[my_sys[i]]->peakoutput_user) = offset;
		peak_params = (CustomParamsMaps*) (sys[my_sys[i]]->peakoutput_user);
		peak_params->forecast_time = forecast_time;
		peak_params->period = period;
	}
}

ForecastData* Init_ForecastData(char* fcst_filename)
{
	FILE* inputfile = NULL;
	ForecastData* Forecaster;
	int errorcode,string_size = 128,length;

	MPI_Barrier(MPI_COMM_WORLD);

	if(my_rank == 0)
	{
		//Open file
		inputfile = fopen(fcst_filename,"r");
		errorcode = 0;
		if(!inputfile)
		{
			printf("[%i]: Error opening file %s.\n",fcst_filename);
			errorcode = 1;
		}
	}

	//Check if forecast file was openned
	MPI_Bcast(&errorcode,1,MPI_INT,0,MPI_COMM_WORLD);
	if(errorcode)	return NULL;

	//Reserve space
	Forecaster = (ForecastData*) malloc(sizeof(ForecastData));
	Forecaster->model_name = (char*) malloc(string_size*sizeof(char));

	//Read table name
	if(my_rank == 0)
	{
		fscanf(inputfile,"%s",Forecaster->model_name);
		length = strlen(Forecaster->model_name);
	}
	MPI_Bcast(&length,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(Forecaster->model_name,length+1,MPI_CHAR,0,MPI_COMM_WORLD);

	//Clean up
	if(my_rank == 0)
		fclose(inputfile);

	MPI_Barrier(MPI_COMM_WORLD);
	return Forecaster;
}

void Free_ForecastData(ForecastData** Forecaster)
{
	free((*Forecaster)->model_name);
	free(*Forecaster);
	*Forecaster = NULL;
}

