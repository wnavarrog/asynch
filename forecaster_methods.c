#include "forecaster_methods.h"


//Deletes all future values from a set of partitioned tables.
int DeleteFutureValues(ConnData* conninfo,unsigned int num_tables,UnivVars* GlobalVars,char* table_name,char* model_name,unsigned int clear_after,unsigned int equal)
{
	PGresult* res;
	int i,del_table,error = 0,last_table;
	time_t current_time,clear_time = (time_t) clear_after;
	char* query = (char*) malloc(GlobalVars->query_size*sizeof(char));
	char operation[3];
	if(equal)	sprintf(operation,">=");
	else		sprintf(operation,">");

	//Find the last table with values to destroy
	time(&current_time);
	del_table = (int) (difftime(current_time,clear_time) / (60.0*60.0*24.0) + 1e-6);

	//Truncate any tables with index lower than del_table
	last_table = (num_tables < del_table) ? num_tables : del_table;
	for(i=0;i<last_table;i++)
	{
		sprintf(query,"TRUNCATE %s_%s_%i;",table_name,model_name,i);
		res = PQexec(conninfo->conn,query);
		error = CheckResError(res,"truncating table");
		PQclear(res);
		if(error)	return error;
	}

	//Delete del_table
	if(del_table < num_tables)
	{
		sprintf(query,"DELETE FROM %s_%s_%i WHERE forecast_time %s %u;",table_name,model_name,i,operation,clear_after);
		res = PQexec(conninfo->conn,query);
		error = CheckResError(res,"deleting from table");
		PQclear(res);
		if(error)	return error;
	}

	//Clean up
	free(query);
	return error;
}


//Checks if the time is right to perform maintainance on the database.
void PerformTableMaintainance(ConnData* conninfo_hydros,UnivVars* GlobalVars,ForecastData* Forecaster,short int* vac,short unsigned int hr1,unsigned int num_tables,char* tablename)
{
	int j;
	time_t start,stop;
	PGresult* res;
	char query[GlobalVars->query_size];
	struct tm* timeinfo;
	time(&start);
	timeinfo = localtime(&start);

	if(timeinfo->tm_hour == hr1 && *vac == 0)
	{
		printf("[%i]: Performing maintainance. Current time is %s",my_rank,asctime(timeinfo));

		//Adjust partitioned hydroforecast tables
		ConnectPGDB(conninfo_hydros);
		sprintf(query,"DROP TABLE %s_%s_%u;",tablename,Forecaster->model_name,num_tables-1);
		res = PQexec(conninfo_hydros->conn,query);
		CheckResError(res,"dropping end archive table");
		PQclear(res);

		for(j=num_tables-2;j>=0;j--)
		{
			sprintf(query,"ALTER TABLE %s_%s_%i RENAME TO %s_%s_%i;",tablename,Forecaster->model_name,j,tablename,Forecaster->model_name,j+1);
			res = PQexec(conninfo_hydros->conn,query);
			CheckResError(res,"renaming table");
			PQclear(res);

			sprintf(query,"ALTER INDEX idx_%s_%s_%i_forecast_time_link_id RENAME TO idx_%s_%s_%i_forecast_time_link_id;",tablename,Forecaster->model_name,j,tablename,Forecaster->model_name,j+1);
			res = PQexec(conninfo_hydros->conn,query);
			CheckResError(res,"renaming index on archive table");
			PQclear(res);

/*
			sprintf(query,"ALTER INDEX idx_%s_%s_%i_forecast_time RENAME TO idx_%s_%s_%i_forecast_time;",tablename,Forecaster->model_name,j,tablename,Forecaster->model_name,j+1);
			res = PQexec(conninfo_hydros->conn,query);
			CheckResError(res,"renaming forecast_time index on archive table");
			PQclear(res);

			sprintf(query,"ALTER INDEX idx_%s_%s_%i_link_id RENAME TO idx_%s_%s_%i_link_id;",tablename,Forecaster->model_name,j,tablename,Forecaster->model_name,j+1);
			res = PQexec(conninfo_hydros->conn,query);
			CheckResError(res,"renaming link_id index on archive table");
			PQclear(res);
*/
		}

		sprintf(query,"CREATE TABLE %s_%s_0 ( ) INHERITS (master_%s_%s);",tablename,Forecaster->model_name,tablename,Forecaster->model_name);
		res = PQexec(conninfo_hydros->conn,query);
		CheckResError(res,"creating archive table 0");
		PQclear(res);

		sprintf(query,"CREATE INDEX idx_%s_%s_0_forecast_time_link_id ON %s_%s_0 USING btree (forecast_time,link_id);",tablename,Forecaster->model_name,tablename,Forecaster->model_name);
		res = PQexec(conninfo_hydros->conn,query);
		CheckResError(res,"creating index on archive table 0");
		PQclear(res);

/*
		sprintf(query,"CREATE INDEX idx_%s_%s_0_forecast_time ON %s_%s_0 USING btree (forecast_time);",tablename,Forecaster->model_name,tablename,Forecaster->model_name);
		res = PQexec(conninfo_hydros->conn,query);
		CheckResError(res,"creating forecast_time index on archive table 0");
		PQclear(res);

		sprintf(query,"CREATE INDEX idx_%s_%s_0_link_id ON %s_%s_0 USING btree (link_id);",tablename,Forecaster->model_name,tablename,Forecaster->model_name);
		res = PQexec(conninfo_hydros->conn,query);
		CheckResError(res,"creating link_id index on archive table 0");
		PQclear(res);
*/

		DisconnectPGDB(conninfo_hydros);

		//Set flag and print the total time
		*vac = 1;
		time(&stop);
		printf("[%i]: Database cleanup complete. Total time %.2f.\n\n",my_rank,difftime(stop,start));
	}
	else	if(timeinfo->tm_hour != hr1)	*vac = 0;
}

//Checks that the timestamps of the Peakforecast tables match up correctly with the trigger.
//If not, the tables are adjusted.
//Note: The total number of tables is hard coded below.
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
		sprintf(query,"ALTER TABLE %s_%s_%i RENAME TO %s_%s_%i;",tablename,Forecaster->model_name,j,tablename,Forecaster->model_name,j+diff_table_index);
		res = PQexec(conninfo->conn,query);
		CheckResError(res,"renaming table");
		PQclear(res);

		sprintf(query,"ALTER INDEX idx_%s_%s_%i_forecast_time_link_id RENAME TO idx_%s_%s_%i_forecast_time_link_id;",tablename,Forecaster->model_name,j,tablename,Forecaster->model_name,j+diff_table_index);
		res = PQexec(conninfo->conn,query);
		CheckResError(res,"renaming index on archive table");
		PQclear(res);

/*
		sprintf(query,"ALTER INDEX idx_%s_%s_%i_forecast_time RENAME TO idx_%s_%s_%i_forecast_time;",tablename,Forecaster->model_name,j,tablename,Forecaster->model_name,j+diff_table_index);
		res = PQexec(conninfo->conn,query);
		CheckResError(res,"renaming forecast_time index on archive table");
		PQclear(res);

		sprintf(query,"ALTER INDEX idx_%s_%s_%i_link_id RENAME TO idx_%s_%s_%i_link_id;",tablename,Forecaster->model_name,j,tablename,Forecaster->model_name,j+diff_table_index);
		res = PQexec(conninfo->conn,query);
		CheckResError(res,"renaming link_id index on archive table");
		PQclear(res);
*/
	}

	//Create new tables
	for(i=0;i<diff_table_index;i++)
	{
		sprintf(query,"CREATE TABLE %s_%s_%u ( ) INHERITS (master_%s_%s);",tablename,Forecaster->model_name,i,tablename,Forecaster->model_name);
		res = PQexec(conninfo->conn,query);
		CheckResError(res,"creating table");
		PQclear(res);

		sprintf(query,"CREATE INDEX idx_%s_%s_%u_forecast_time_link_id ON %s_%s_%u USING btree (forecast_time,link_id);",tablename,Forecaster->model_name,i,tablename,Forecaster->model_name,i);
		res = PQexec(conninfo->conn,query);
		CheckResError(res,"creating index on archive table");
		PQclear(res);

/*
		sprintf(query,"CREATE INDEX idx_%s_%s_%u_forecast_time ON %s_%s_%u USING btree (forecast_time);",tablename,Forecaster->model_name,i,tablename,Forecaster->model_name,i);
		res = PQexec(conninfo->conn,query);
		CheckResError(res,"creating forecast_time index on archive table");
		PQclear(res);

		sprintf(query,"CREATE INDEX idx_%s_%s_%u_link_id ON %s_%s_%u USING btree (link_id);",tablename,Forecaster->model_name,i,tablename,Forecaster->model_name,i);
		res = PQexec(conninfo->conn,query);
		CheckResError(res,"creating link_id index on archive table");
		PQclear(res);
*/
	}
}

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



ForecastData* Init_ForecastData(char* fcst_filename,unsigned int string_size)
{
	FILE* inputfile = NULL;
	ForecastData* Forecaster;
	int errorcode,length;
	char end_char;

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
	MPI_Bcast(&length,1,MPI_UNSIGNED,0,MPI_COMM_WORLD);
	MPI_Bcast(Forecaster->model_name,length+1,MPI_CHAR,0,MPI_COMM_WORLD);

	//Read if data is displayed on ifis
	if(my_rank == 0)
		fscanf(inputfile,"%hi",&(Forecaster->ifis_display));
	MPI_Bcast(&(Forecaster->ifis_display),1,MPI_SHORT,0,MPI_COMM_WORLD);

	//Read which forcing index is used for forecasting
	if(my_rank == 0)
		fscanf(inputfile,"%u",&(Forecaster->forecasting_forcing));
	MPI_Bcast(&(Forecaster->forecasting_forcing),1,MPI_UNSIGNED,0,MPI_COMM_WORLD);

	//Read number of rainfall steps to use per forecast
	if(my_rank == 0)
		fscanf(inputfile,"%u",&(Forecaster->num_rainsteps));
	MPI_Bcast(&(Forecaster->num_rainsteps),1,MPI_UNSIGNED,0,MPI_COMM_WORLD);

	//Read halt filename
	Forecaster->halt_filename = (char*) malloc(string_size*sizeof(char));
	if(my_rank == 0)
	{
		fscanf(inputfile,"%s",Forecaster->halt_filename);
		length = strlen(Forecaster->halt_filename);
	}
	MPI_Bcast(&length,1,MPI_UNSIGNED,0,MPI_COMM_WORLD);
	MPI_Bcast(Forecaster->halt_filename,length+1,MPI_CHAR,0,MPI_COMM_WORLD);

	//Read ending mark
	if(my_rank == 0)
		fscanf(inputfile,"%s",&end_char);
	MPI_Bcast(&end_char,1,MPI_CHAR,0,MPI_COMM_WORLD);

	//Clean up
	if(my_rank == 0)
		fclose(inputfile);

	MPI_Barrier(MPI_COMM_WORLD);
	if(end_char != '#')
	{
		if(my_rank == 0)
			printf("[%i]: Error: Ending mark not seen in %s.\n",my_rank,fcst_filename);
		return NULL;
	}
	return Forecaster;
}

void Free_ForecastData(ForecastData** Forecaster)
{
	free((*Forecaster)->model_name);
	free((*Forecaster)->halt_filename);
	free(*Forecaster);
	*Forecaster = NULL;
}

