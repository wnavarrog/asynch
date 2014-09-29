#include "rkmethods.h"
#include "system.h"
#include "comm.h"
#include "riversys.h"
#include "processdata.h"
#include <time.h>
#include <libpq-fe.h>
#include "mpi.h"
#include "misc.h"
#include "rainfall.h"
#include <string.h>
#include "solvers.h"
#include <unistd.h>

int my_rank;
int np;

short int CheckFinished();
void PerformTableMaintainance(ConnData* conninfo,short int* vac,short unsigned int hr1,short unsigned int hr2,unsigned int archive_window,unsigned int num_tables,emaildata* erroremail);
void CheckPeakforecastTable(ConnData* conninfo,unsigned int num_tables,emaildata* erroremail);

int main(int argc,char* argv[])
{
	//Initialize MPI stuff
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
	MPI_Comm_size(MPI_COMM_WORLD,&np);

	//Parse input
	if(argc < 2)
	{
		if(my_rank == 0)
		{
			printf("Command line parameter required:  A universal variable file (.gbl).\n");
			printf("\n");
		}
		MPI_Finalize();
		return 1;
	}

	//Declare variables
	double total_time = 0.0;
	time_t start,start2,stop;
	unsigned int *my_sys,my_N,my_max_nodes,N,nummethods,my_save_size,save_size,i,j,k,l,*save_list,peaksave_size;
	int *assignments;
	char rkdfilename[256];
	FILE* outputfile = NULL;
	short int *getting = NULL;
	TransData* my_data = NULL;
	RKMethod** AllMethods;
	Link** sys;
	Link* current;
	unsigned int** id_to_loc;
	ErrorData* GlobalErrors;
	UnivVars* GlobalVars;
	TempStorage* workspace;
	ConnData* conninfo = CreateConnData("dbname=rm_model host=s-iihr58.iihr.uiowa.edu user=rainfall_feed password=r!Ain2012");
	emaildata* erroremail = InitializeEmail("3092315446@vtext.com","ASYNCHPERSIS");
	char* query = conninfo->query;
	MPI_Status status;
	PGresult *res;

	start = time(NULL);

	//Read in .gbl file
	GlobalVars = Read_Global_Data(argv[1],&GlobalErrors,conninfo,rkdfilename);
	if(GlobalVars == NULL)	return 1;

	//Read in remaining data from files
	sys = Create_River_System_parallel(rkdfilename,&N,&my_sys,&my_N,&my_max_nodes,&my_data,&assignments,&getting,&AllMethods,&nummethods,GlobalVars,GlobalErrors,&save_list,&my_save_size,&save_size,&peaksave_size,&id_to_loc,conninfo,&workspace);
	if(sys == NULL)		return 1;


	//Determine which links to store peakflow data, and reserve space for sending peak data to database
	#ifdef PRINTPEAKFLOW
		if(my_rank == 0)	printf("Finding links with order > 1...\n");
		start2 = time(NULL);

		unsigned int *my_links,*all_my_N;
		if(my_rank == 0)	all_my_N = (unsigned int*) malloc(np*sizeof(unsigned int));

		//Gather the total number of links each proc has
		MPI_Gather(&my_N,1,MPI_UNSIGNED,all_my_N,1,MPI_UNSIGNED,0,MPI_COMM_WORLD);

		//Have process 0 download and sort the Horton ordering
		if(my_rank == 0)
		{
			unsigned int* num_links = (unsigned int*) calloc(np,sizeof(unsigned int));
			unsigned int** link_locs = (unsigned int**) malloc(np*sizeof(unsigned int*));
			for(i=0;i<(unsigned int)np;i++)
				link_locs[i] = (unsigned int*) malloc(all_my_N[i]*sizeof(unsigned int));

			sprintf(query,"SELECT link_id FROM master_update WHERE h_order > 1;");
			CheckConnConnection(conninfo);
			res = PQexec(conninfo->conn,query);
			CheckResError(res,"checking Horton ordering",erroremail);

			//Find the locations
			j = PQntuples(res);
			for(i=0;i<j;i++)
			{
				k = find_link_by_idtoloc(atoi(PQgetvalue(res,i,0)),id_to_loc,N);
				l = assignments[k];
				link_locs[l][num_links[l]] = k;
				num_links[l]++;
			}

			//Send the links
			for(i=1;i<(unsigned int)np;i++)
			{
				MPI_Send(&(num_links[i]),1,MPI_UNSIGNED,i,0,MPI_COMM_WORLD);
				MPI_Send(link_locs[i],num_links[i],MPI_UNSIGNED,i,0,MPI_COMM_WORLD);
			}

			//Cleanup
			k = num_links[0];
			my_links = link_locs[0];
			for(i=1;i<(unsigned int)np;i++)
				free(link_locs[i]);
			free(link_locs);
			free(num_links);
			free(all_my_N);
		}
		else
		{
			MPI_Recv(&k,1,MPI_UNSIGNED,0,0,MPI_COMM_WORLD,&status);
			my_links = (unsigned int*) malloc(k*sizeof(unsigned int));
			MPI_Recv(my_links,k,MPI_UNSIGNED,0,0,MPI_COMM_WORLD,&status);
		}

		//Set the flags and peakflowdata buffer
		for(i=0;i<k;i++)
			sys[my_links[i]]->peak_flag = 1;
		char* peakflowdata = (char*) malloc(k*(10+10+20+10 + 4*1)*sizeof(char));
		free(my_links);

		stop = time(NULL);
		if(my_rank == 0)	printf("Time to find non-order 1 links: %f\n",difftime(stop,start2));
	#endif

	MPI_Barrier(MPI_COMM_WORLD);
	stop = time(NULL);
	total_time = difftime(stop,start);
	if(my_rank == 0)	printf("Finished reading files. Total time: %f\n",total_time);

	//Reserve space for backups
	VEC** backup = (VEC**) malloc(N*sizeof(VEC*));
	for(i=0;i<N;i++)
	{
		if(assignments[i] == my_rank || getting[i] == 1)
			backup[i] = v_get(GlobalVars->dim);
		else
			backup[i] = NULL;
	}

	//Make sure everyone is good before getting down to it...
	printf("Process %i (%i total) is good to go with %i links.\n",my_rank,np,my_N);
	MPI_Barrier(MPI_COMM_WORLD);
	start = time(NULL);
/*
	//Make the initial solve
	AsynchSolver(sys,N,my_sys,my_N,my_max_nodes,GlobalVars,assignments,id_to_loc,workspace,conninfo,my_data,0,outputfile);

	//Stop the clock
	MPI_Barrier(MPI_COMM_WORLD);
	stop = time(NULL);
	total_time += difftime(stop,start);

	//Output some data
	if(my_rank == 0)
	{
		printf("%i: The answer at ID %i at time %.12f is\n",my_rank,sys[my_sys[0]]->ID,sys[my_sys[0]]->last_t);
		Print_Vector(sys[my_sys[0]]->list->tail->y_approx);
		printf("Total time for calculations: %f\n",difftime(stop,start));
	}
*/
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
	double forecast_time = 5.0*24*60;	//Time (mins) in future to make forecasts
	unsigned int num_peakflow_tables = 10;

	unsigned int nextraintime,repeat_for_errors;
	short int halt = 0;
	unsigned int halt_time = GlobalVars->last_file;
	int isnull;
	int message_buffer[2];
	short int vac = 0;	//0 if no vacuum has occured, 1 if vacuum has occured (during a specific hour)
	short unsigned int rainfirst = 1;	//1 if rainfall is present at the first time, 0 if not !!!! Should this be 1 to start? !!!!
	//increment = 1;
	//unsigned int last_file = GlobalVars->last_file;
	unsigned int last_file = GlobalVars->first_file;
	unsigned int first_file = GlobalVars->first_file;
	k = 0;
	for(i=0;i<N;i++)
		if(backup[i] != NULL)	v_copy(sys[i]->list->tail->y_approx,backup[i]);

	char filename[256];
	char additional[256];
	//FILE* myoutputfile = NULL;
/*
	sprintf(filename,"TempData/persis/modeltest/output%.3i",my_rank);
	FILE* myoutputfile = fopen(filename,"w");
	if(!myoutputfile)
	{
		printf("%i: Could not open file %s.\n",my_rank,filename);
		abort();
	}
*/

	//Make some initializations to the database
/*
	if(my_rank == 0)
	{
		printf("Making initializations to tables.\n");
		start = time(NULL);

		//Make sure the hydrographs table exists
		sprintf(query,"SELECT 1 FROM pg_class WHERE relname='hydrographs';");
		res = PQexec(conninfo->conn,query);
		if(!PQntuples(res))
		{
			PQclear(res);
			//sprintf(query,"CREATE TABLE hydrographs(link_id int,time int,discharge real,primary key(link_id,time)); ALTER TABLE hydrographs SET (autovacuum_enabled = false, toast.autovacuum_enabled = false);");
			sprintf(query,"CREATE TABLE hydrographs(link_id int,time int,discharge real); ALTER TABLE hydrographs SET (autovacuum_enabled = false, toast.autovacuum_enabled = false);");
			res = PQexec(conninfo->conn,query);
			CheckResError(res,"creating hydrographs table",erroremail);
		}
		PQclear(res);

		//Make sure the peakforecast tables are set correctly
		CheckPeakforecastTable(conninfo,num_peakflow_tables,erroremail);

		//Clear all future peakflows
		sprintf(query,"DELETE FROM master_peakforecast WHERE start_time >= %u;",last_file);
		res = PQexec(conninfo->conn,query);
		CheckResError(res,"deleting future peakflows in archive",erroremail);
		PQclear(res);

		stop = time(NULL);
		printf("Total time to initialize tables: %.2f.\n",difftime(stop,start));
	}
*/
	MPI_Barrier(MPI_COMM_WORLD);

	//Start the main loop
	while(first_file < halt_time)
	{
		if(my_rank == 0)	printf("\n\nPass %u\n",k);

		//Flush communication buffers
		Flush_TransData(my_data);

		//Make some initializations to each link
		for(i=0;i<N;i++)	//Set time to 0.0
		{
			current = sys[i];
			if(current->list != NULL)
			{
				while(current->current_iterations > 1)
				{
					Remove_Head_Node(current->list);
					(current->current_iterations)--;
				}
				current->list->head->t = 0.0;
				current->last_t = 0.0;
				current->steps_on_diff_proc = 1;
				current->iters_removed = 0;
				current->rejected = 0;
				if(current->numparents == 0)	current->ready = 1;
				else				current->ready = 0;
				v_copy(backup[i],current->list->head->y_approx);

				//Reset the next_save time
				if(current->save_flag)
				{
					current->next_save = current->print_time;
					current->disk_iterations = 1;
				}

				//Reset peak flow information
				current->peak_time = 0.0;
				v_copy(current->list->head->y_approx,current->peak_value);
			}
		}

		//Check if a vacuum should be done
		//This will happen at hr1 (3 a.m. currently, regardless of DST, I hope)
		//if(my_rank == 0)	PerformTableMaintainance(conninfo,&vac,hr1,hr2,first_file - history_time,num_peakflow_tables,erroremail);

		//Make sure all buffer flushing is done
		//!!!! Is this really needed anymore? !!!!
		MPI_Barrier(MPI_COMM_WORLD);

		//Dump data for debugging and recovery
		//if(k % 288 == 0)		DataDump(sys,N,assignments,GlobalVars,last_file);

		//Make some initializations
		GlobalVars->raindb_start_time = last_file;
		first_file = last_file;
		last_file = last_file + (unsigned int) GlobalVars->file_time * 60;

		//Open the next output file
		sprintf(additional,"%u",first_file);
		GlobalVars->maxtime = forecast_time;	//!!!! This may not be needed... !!!!
		outputfile = PrepareTempFiles(sys,N,assignments,GlobalVars,save_list,save_size,my_save_size,filename,additional,id_to_loc);
/*
		if(myoutputfile)	fclose(myoutputfile);
		sprintf(filename,"TempData/persis/modeltest/output_%u_%.3i",first_file,my_rank);
		myoutputfile = fopen(filename,"w");
		if(!myoutputfile)
		{
			printf("%i: Could not open file %s.\n",my_rank,filename);
			abort();
		}
*/

		//Find the next time where rainfall occurs
		do
		{
			if(my_rank == 0)
			{
				//Check if connection to SQL database is still good
				CheckConnConnection(conninfo);

				//Find the next rainfall time
				time(&start);
				if(GlobalVars->outletlink == 0)
					sprintf(query,"SELECT min(unix_time) FROM link_rain5 where unix_time >= %u AND link_id > 1 AND link_id < %u;",last_file,N+2);
				else
					sprintf(query,"WITH subbasin AS (SELECT nodeX.link_id \
						FROM master_network AS nodeX, master_network AS parentX \
						WHERE (nodeX.left BETWEEN parentX.left AND parentX.right) AND parentX.link_id = %u ORDER BY nodeX.link_id) \
						SELECT min(unix_time) FROM link_rain5 WHERE unix_time >= %u \
							AND link_id IN (SELECT link_id FROM subbasin);",GlobalVars->outletlink,last_file);
				res = PQexec(conninfo->conn,query);
				CheckResError(res,"checking for new rainfall data",erroremail);
				time(&stop);
				printf("Total time to check for new rainfall data: %f.\n",difftime(stop,start));

				message_buffer[0] = PQgetisnull(res,0,0);
				if(!message_buffer[0])	message_buffer[1] = atoi(PQgetvalue(res,0,0));
				else			message_buffer[1] = 0;
				PQclear(res);
			}
			MPI_Bcast(message_buffer,2,MPI_INT,0,MPI_COMM_WORLD);
			isnull = message_buffer[0];

			if(!isnull)
				nextraintime = message_buffer[1];
			else
			{
/*
				if(my_rank == 0)
				{
					printf("No rainfall values returned from SQL database. %u %u\n",last_file,isnull);
					PerformTableMaintainance(conninfo,&vac,hr1,hr2,first_file - history_time,num_peakflow_tables,erroremail);
				}

				halt = CheckFinished();
				if(halt)
					DataDump(sys,N,assignments,GlobalVars,first_file);
				else
				{
					fflush(stdout);
					sleep(wait_time);
				}
*/
			}
		} while(isnull == 1 && !halt);

//		if(halt)	break;

		//Read in next set of rainfall data
		if(nextraintime == last_file)
		{
if(my_rank == 0)
printf("1first: %u last: %u\n",first_file,last_file);
			rainfirst = 1;
		}
		else
		{
			rainfirst = 0;
			SetRain0(sys,my_N,GlobalVars->maxtime,my_sys);
if(my_rank == 0)
printf("2first: %u last: %u\n",first_file,last_file);
		}

		//Initialize some data for the first phase of calculations
		GlobalVars->maxtime = GlobalVars->file_time;
		conninfo->time_offset = first_file;

		MPI_Barrier(MPI_COMM_WORLD);
		time(&start);
		Create_Rain_Database(sys,N,my_N,GlobalVars,my_sys,assignments,conninfo,first_file,last_file,id_to_loc,forecast_time);
		MPI_Barrier(MPI_COMM_WORLD);
		time(&stop);
		if(my_rank == 0)	printf("Total time to download rain data: %f.\n",difftime(stop,start));

		//Make some universal initializations
		current = sys[my_sys[my_N-1]];
		#ifdef PRINT2DATABASE
			conninfo->submission_content = 0;
			conninfo->submission[0] = '\0';
		#endif

		//Set a new step size
		for(i=0;i<my_N;i++)
			sys[my_sys[i]]->h = InitialStepSize(sys[my_sys[i]]->last_t,sys[my_sys[i]],GlobalVars,workspace);		

		//Make the first phase of calculations
		start = time(NULL);
		AsynchSolverPersis(sys,my_sys,my_N,my_max_nodes,GlobalVars,assignments,workspace,conninfo,my_data,1,outputfile);
		MPI_Barrier(MPI_COMM_WORLD);
		if(my_rank == 0)
		{
			stop = time(NULL);
			printf("Time for first phase calculations: %.2f\n",difftime(stop,start));
		}

		//Flush communication buffers
		Flush_TransData(my_data);

		//Reset the links (mostly) and make a backup for the second phase
		for(i=0;i<N;i++)	//Time does not change!
		{
			current = sys[i];
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

		//Make sure all buffer flushing is done
		//!!!! Do I really need this anymore? !!!!
		MPI_Barrier(MPI_COMM_WORLD);

		//Make second phase calculations
		start = time(NULL);
		GlobalVars->maxtime = forecast_time;
		AsynchSolverPersis(sys,my_sys,my_N,my_max_nodes,GlobalVars,assignments,workspace,conninfo,my_data,1,outputfile);
		MPI_Barrier(MPI_COMM_WORLD);
		if(my_rank == 0)
		{
			stop = time(NULL);
			printf("Time for second phase calculations: %.2f\n",difftime(stop,start));
		}

		//Output some data
		if(my_rank == 0)
		{
			printf("%i: The answer at ID %i at time %.12f is\n",my_rank,sys[my_sys[0]]->ID,sys[my_sys[0]]->last_t);
			Print_Vector(sys[my_sys[0]]->list->tail->y_approx);
		}

		//Upload the peak data to the database
/*
		#ifdef PRINTPEAKFLOW
		//if(rainfirst == 1)
		{
			start = time(NULL);

			//Adjust the table peakforecast
			if(my_rank == 0)
			{
				//Make sure database connection is still good
				CheckConnConnection(conninfo);

				//Clear the future peakflows
				//sprintf(query,"DELETE FROM peakforecast WHERE start_time >= %u;",first_file);
				//res = PQexec(conninfo->conn,query);
				//CheckResError(res,"deleting future peakflows");
				//PQclear(res);
//printf("%i: Deleted all entries in peakforecast with time after %u\n",my_rank,first_file);

				//Clear the past peakflows
				//This keeps hydrographs for 5 days in the past
				//sprintf(query,"DELETE FROM peakforecast WHERE start_time < %u;",first_file - history_time);
				//res = PQexec(conninfo->conn,query);
				//CheckResError(res,"deleting past peakflows");
				//PQclear(res);
//printf("%i: Deleted all entries in peakforecast with time before %u\n",my_rank,first_file - history_time);
			}
			MPI_Barrier(MPI_COMM_WORLD);

			//Loop over procs
			for(l=0;l<(unsigned int)np+20;l+=20)
			{
				if(l <= (unsigned int)my_rank && (unsigned int)my_rank < l + 20)
				{
					//Connect to the database
					if(my_rank != 0)	ConnectPGDB(conninfo);

					repeat_for_errors = 1;
					while(repeat_for_errors)
					{
						repeat_for_errors = 0;

						//Delete temp table
						sprintf(query,"SELECT 1 FROM pg_class WHERE relname='temp_prediction_%i';",my_rank);
						res = PQexec(conninfo->conn,query);
						repeat_for_errors = repeat_for_errors || CheckResError(res,"checking existence of temp_prediction table",erroremail);
						if(PQntuples(res))
						{
							PQclear(res);
							sprintf(query,"DROP TABLE temp_prediction_%i;",my_rank);
							res = PQexec(conninfo->conn,query);
							repeat_for_errors = repeat_for_errors || CheckResError(res,"dropping temp_prediction table",erroremail);
						}
						PQclear(res);
			
						//Create new temp table
						sprintf(query,"CREATE TABLE temp_prediction_%i(link_id int,peak_time int,peak_discharge double precision,start_time int);",my_rank);
						res = PQexec(conninfo->conn,query);
						repeat_for_errors = repeat_for_errors || CheckResError(res,"creating temp_prediction table",erroremail);
						PQclear(res);

						//Load peak data into a string
						peakflowdata[0] = '\0';
						for(i=0;i<my_N;i++)
						{
							current = sys[my_sys[i]];
							if(current->peak_flag)
							{
								double peakdis = (current->peak_value->ve[0] > 1e-36) ? current->peak_value->ve[0] : 0.0;
								sprintf(query,"%u,%u,%.12e,%u\n",current->ID,first_file + (unsigned int)current->peak_time*60,peakdis,first_file);
								strcat(peakflowdata,query);
							}
						}

						//Copy peakflowdata to database
						CheckConnConnection(conninfo);
						sprintf(query,"COPY temp_prediction_%i FROM STDIN WITH DELIMITER ',';",my_rank);
						printf("%s",PQresultErrorMessage(res));
						res = PQexec(conninfo->conn,query);
						PQclear(res);
						i = PQputCopyData(conninfo->conn,peakflowdata,strlen(peakflowdata));
						j = PQputCopyEnd(conninfo->conn,NULL);
						if(i != 1 || j != 1)
						{
							printf("%i: Returned i = %i j = %i\n",my_rank,i,j);
							repeat_for_errors = 1;
						}

						//Move data to permanent table
						sprintf(query,"INSERT INTO master_peakforecast SELECT * FROM temp_prediction_%i; DROP TABLE temp_prediction_%i;",my_rank,my_rank);
						res = PQexec(conninfo->conn,query);
						repeat_for_errors = repeat_for_errors || CheckResError(res,"copying/dropping temp_prediction table",erroremail);
						PQclear(res);

						if(repeat_for_errors)
						{
							printf("%i: Attempting resend of peakforecast data.\n",my_rank);
							CheckConnConnection(conninfo);
						}
					}

					//Disconnect from database
					if(my_rank != 0)	PQfinish(conninfo->conn);
				} //End if

				MPI_Barrier(MPI_COMM_WORLD);

			} //End for loop

			if(my_rank == 0)
			{
				stop = time(NULL);
				printf("Total time to transfer peak flow data: %.2f\n",difftime(stop,start));
			}
			fflush(stdout);
		}
		#endif
*/
/*
		//Upload the hydrographs to the database
		#ifdef PRINT2DATABASE
		//if(rainfirst == 1)
		{
			start = time(NULL);

			//Adjust the table hydrographs
			if(my_rank == 0)
			{
				//Make sure database connection is still good
				CheckConnConnection(conninfo);

				//Clear the future hydrographs
				//sprintf(query,"DELETE FROM hydrographs WHERE time > %u; ALTER SEQUENCE hydrographs_tid_seq RESTART;",first_file);
				sprintf(query,"DELETE FROM hydrographs WHERE time > %u;",first_file);
				res = PQexec(conninfo->conn,query);
				CheckResError(res,"deleting future hydrographs",erroremail);
				PQclear(res);
printf("%i: Deleted all entries with time after %u\n",my_rank,first_file);

				//Clear the past hydrographs
				sprintf(query,"DELETE FROM hydrographs WHERE time < %u;",first_file - history_time);
				res = PQexec(conninfo->conn,query);
				CheckResError(res,"deleting past hydrographs",erroremail);
				PQclear(res);
printf("%i: Deleted all entries with time before %u\n",my_rank,first_file - history_time);
			}
			MPI_Barrier(MPI_COMM_WORLD);

			//Loop over procs
			for(l=0;l<np+20;l+=20)
			{
				if(l <= my_rank && my_rank < l + 20)
				{
					//Connect to the database
					if(my_rank != 0)	ConnectPGDB(conninfo);

					repeat_for_errors = 1;
					while(repeat_for_errors)
					{
						repeat_for_errors = 0;

						//Delete temp table if it exists
						sprintf(query,"SELECT 1 FROM pg_class WHERE relname='temp_hydrographs_%i';",my_rank);
						res = PQexec(conninfo->conn,query);
						repeat_for_errors = repeat_for_errors || CheckResError(res,"checking existence of temp_hydrographs table",erroremail);
						if(PQntuples(res))
						{
							PQclear(res);
							sprintf(query,"DROP TABLE temp_hydrographs_%i;",my_rank);
							res = PQexec(conninfo->conn,query);
							repeat_for_errors = repeat_for_errors || CheckResError(res,"dropping temp_hydrographs table",erroremail);
						}
						PQclear(res);

						//Create new temp table
						sprintf(query,"CREATE TABLE temp_hydrographs_%i(link_id int,time int,discharge real,primary key(link_id,time));",my_rank);
						res = PQexec(conninfo->conn,query);
						repeat_for_errors = repeat_for_errors || CheckResError(res,"creating temp_hydrographs table",erroremail);
						PQclear(res);

						//Copy hydrographs to database
						CheckConnConnection(conninfo);
						sprintf(query,"COPY temp_hydrographs_%i FROM STDIN WITH DELIMITER ',';",my_rank);
						printf("%s",PQresultErrorMessage(res));
						res = PQexec(conninfo->conn,query);
						PQclear(res);
						i = PQputCopyData(conninfo->conn,conninfo->submission,conninfo->submission_content);
						j = PQputCopyEnd(conninfo->conn,NULL);
						if(i != 1 || j != 1)
						{
							printf("%i: Returned i = %i j = %i\n",my_rank,i,j);
							repeat_for_errors = 1;
						}

						//Move data to permanent table
						sprintf(query,"INSERT INTO hydrographs SELECT * FROM temp_hydrographs_%i; DROP TABLE temp_hydrographs_%i;",my_rank,my_rank);
						res = PQexec(conninfo->conn,query);
						repeat_for_errors = repeat_for_errors || CheckResError(res,"copying/dropping hydrographs",erroremail);
						PQclear(res);

						if(repeat_for_errors)
						{
							printf("%i: Attempting resend of hydrographs data.\n",my_rank);
							CheckConnConnection(conninfo);
						}
					}

					//Disconnect from database
					if(my_rank != 0)	PQfinish(conninfo->conn);
				} //End if

				MPI_Barrier(MPI_COMM_WORLD);

			} //End for loop

			if(my_rank == 0)
			{
				stop = time(NULL);
				printf("Total time to transfer hydrograph data: %.2f\n",difftime(stop,start));
			}
			fflush(stdout);
		}
		#endif
*/

		//Dump hydrographs to disk
		//fputs(conninfo->submission,myoutputfile);
		if(outputfile)	fclose(outputfile);
MPI_Barrier(MPI_COMM_WORLD);
		Process_Data(sys,GlobalVars,N,save_list,save_size,my_save_size,id_to_loc,assignments,additional);

		//Check if program has received a terminate signal
		k++;

/*
		halt = CheckFinished();

		//If stopping, make a .rec file
		if(halt)
		{
			for(i=0;i<N;i++)
			{
				current = sys[i];
				if(current->list != NULL)
					v_copy(backup[i],current->list->tail->y_approx);
			}

			DataDump(sys,N,assignments,GlobalVars,last_file);
		}
*/
	}

printf("%i: All done! Cleaning up...\n",my_rank);

	//Free some memory
	TransData_Free(my_data);

	//Cleanup
	ConnData_Free(conninfo);

	//fclose(myoutputfile);

	for(i=0;i<N;i++)	v_free(backup[i]);
	free(backup);
	Destroy_ErrorData(GlobalErrors);
	Destroy_Workspace(workspace,GlobalVars->max_s,GlobalVars->max_parents);
	free(workspace);
	free(assignments);
	free(getting);
	//if(my_save_size > 0)	fclose(outputfile);

	//Last bit of cleanup
	for(i=0;i<N;i++)	Destroy_Link(sys[i],GlobalVars->iter_limit,rkdfilename[0] != '\0',GlobalVars);
	free(sys);
	free(my_sys);
	for(i=0;i<nummethods;i++)	Destroy_RKMethod(AllMethods[i]);
	free(AllMethods);
	free(save_list);
	for(i=0;i<N;i++)
		free(id_to_loc[i]);
	free(id_to_loc);
	#ifdef PRINTPEAKFLOW
		free(peakflowdata);
	#endif
	Destroy_UnivVars(GlobalVars);

	MPI_Finalize();
	return 0;
}

//Checks if the time is right to perform maintainance on the database.
//!!!! Archive window is not currently being used. Make peakforecast match up with this. !!!!
void PerformTableMaintainance(ConnData* conninfo,short int* vac,short unsigned int hr1,short unsigned int hr2,unsigned int archive_window,unsigned int num_tables,emaildata* erroremail)
{
	int j;
	time_t start,start2,stop;
	PGresult* res;
	char* query = conninfo->query;
	struct tm* timeinfo;

	//Check if connection to SQL database is still good
	CheckConnConnection(conninfo);

	time(&start);
	timeinfo = localtime(&start);
	if( (timeinfo->tm_hour == hr1 || timeinfo->tm_hour == hr2) && *vac == 0)	//Main maintainance time
	{
		printf("%i: Performing maintainance. Current time is %s",my_rank,asctime(timeinfo));

		//Vacuum the hydrographs table
		time(&start2);
		sprintf(query,"VACUUM FULL ANALYZE hydrographs;");
		res = PQexec(conninfo->conn,query);
		CheckResError(res,"running hydrographs vacuum",erroremail);
		PQclear(res);
		time(&stop);
		printf("%i: Time to vacuum hydrographs %.2f.\n",my_rank,difftime(stop,start2));

		//Adjust peakforecast table
		if(timeinfo->tm_hour == hr1)
		{
			time(&start2);

			sprintf(query,"DROP TABLE peakforecast_%u;",num_tables-1);
			res = PQexec(conninfo->conn,query);
			CheckResError(res,"dropping end table",erroremail);
			PQclear(res);

			for(j=num_tables-2;j>=0;j--)
			{
				sprintf(query,"ALTER TABLE peakforecast_%i RENAME TO peakforecast_%i;",j,j+1);
				res = PQexec(conninfo->conn,query);
				CheckResError(res,"renaming table",erroremail);
				PQclear(res);
			}

			//sprintf(query,"CREATE TABLE peakforecast_0 ( ) INHERITS (master_peakforecast); ALTER TABLE peakforecast_0 SET (autovacuum_enabled = false, toast.autovacuum_enabled = false);");
			sprintf(query,"CREATE TABLE peakforecast_0 ( ) INHERITS (master_peakforecast);");
			res = PQexec(conninfo->conn,query);
			CheckResError(res,"creating table 0",erroremail);
			PQclear(res);

			time(&stop);
			printf("%i: Time to alter peakforecast table %.2f.\n",my_rank,difftime(stop,start2));
		}

		*vac = 1;
		time(&stop);
		printf("%i: Database cleanup complete. Total time %.2f.\n\n",my_rank,difftime(stop,start));
	}
	else	if(timeinfo->tm_hour != hr1 && timeinfo->tm_hour != hr2)	*vac = 0;
}

//Checks that the timestamps of the Peakforecast tables match up correctly with the trigger.
//If not, the tables are adjusted.
//Note: The total number of tables is hard coded below.
void CheckPeakforecastTable(ConnData* conninfo,unsigned int num_tables,emaildata* erroremail)
{
	unsigned int i,current_time,table_time,table_index,correct_table_index,diff_table_index,last_table_index;
	int diff_time,j;
	PGresult *res;
	char* query = conninfo->query;
	//unsigned int num_tables = 10;
	
	//Grab the current time from the database
	//sprintf(query,"SELECT EXTRACT('epoch' FROM now() AT time zone 'UTC');");
	sprintf(query,"SELECT EXTRACT('epoch' FROM current_date AT time zone 'UTC');");
	res = PQexec(conninfo->conn,query);
	CheckResError(res,"checking current date",erroremail);
	current_time = (unsigned int) rint(atof(PQgetvalue(res,0,0)));
	PQclear(res);


	//Find the first table with something in it
	table_index = num_tables;
	for(i=0;i<num_tables;i++)
	{
		sprintf(query,"SELECT start_time FROM peakforecast_%u LIMIT 1;",i);
		res = PQexec(conninfo->conn,query);
		CheckResError(res,"checking table contents",erroremail);
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
	for(i=num_tables-1;i>last_table_index;i--)
	{
		sprintf(query,"DROP TABLE peakforecast_%u;",i);
		res = PQexec(conninfo->conn,query);
		CheckResError(res,"dropping table",erroremail);
		PQclear(res);
	}

	//Move all the tables
	for(j=last_table_index;j>=0;j--)
	{
		sprintf(query,"ALTER TABLE peakforecast_%i RENAME TO peakforecast_%i;",j,j+diff_table_index);
		res = PQexec(conninfo->conn,query);
		CheckResError(res,"renaming table",erroremail);
		PQclear(res);
	}

	//Create new tables
	for(i=0;i<diff_table_index;i++)
	{
		//sprintf(query,"CREATE TABLE peakforecast_%u ( ) INHERITS (master_peakforecast); ALTER TABLE peakforecast_%u SET (autovacuum_enabled = false, toast.autovacuum_enabled = false);",i,i);
		sprintf(query,"CREATE TABLE peakforecast_%u ( ) INHERITS (master_peakforecast);",i);
		res = PQexec(conninfo->conn,query);
		CheckResError(res,"creating table",erroremail);
		PQclear(res);
	}
}
/*
//Adjusts the peakforecast tables by one.
void AdjustPeakforecastTable(ConnData* conninfo)
{
	int j;
	PGresult *res;
	char* query = conninfo->query;
	unsigned int num_tables = 10;

	//Trash the table at the end
	sprintf(query,"DROP TABLE peakforecast_%u;",num_tables-1);
	res = PQexec(conninfo->conn,query);
	CheckResError(res,"dropping end table",erroremail);
	PQclear(res);

	//Move all the tables
	for(j=num_tables-2;j>=0;j--)
	{
		sprintf(query,"ALTER TABLE peakforecast_%i RENAME TO peakforecast_%i;",j,j+1);
		res = PQexec(conninfo->conn,query);
		CheckResError(res,"renaming table",erroremail);
		PQclear(res);
	}

	//Create new table
	sprintf(query,"CREATE TABLE peakforecast_0 ( ) INHERITS (master_peakforecast);");
	res = PQexec(conninfo->conn,query);
	CheckResError(res,"creating table 0",erroremail);
	PQclear(res);
}
*/

//Returns 0 if the program should continue, 1 if the program should terminate.
short int CheckFinished()
{
	FILE* inputfile;
	short int halt;
	time_t now;
	struct tm* timeinfo;

	if(my_rank == 0)
	{
		//Check terminate file
		inputfile = fopen("terminate","r");
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

