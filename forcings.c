#include "forcings.h"


Forcing* InitializeForcings()
{
	Forcing* forcing = (Forcing*) malloc(sizeof(Forcing));
	forcing->filename = NULL;
	forcing->GlobalForcing = NULL;
	return forcing;
}


void FreeForcing(Forcing** forcings)
{
	if(forcings && *forcings)
	{
		if((*forcings)->filename)	free((*forcings)->filename);
		Destroy_ForcingData(&((*forcings)->GlobalForcing));
		free(*forcings);
		*forcings = NULL;
	}
}

//GetPasses ************************************************************************************
//Forcings (0 = none, 1 = .str, 2 = binary, 3 = database, 4 = .ustr, 5 = forcasting, 6 = .gz binary, 7 = recurring)

//For flag = 0,1,4
unsigned int PassesOther(Forcing* forcing)
{
	return 1;
}

//For flag = 2,6
unsigned int PassesBinaryFiles(Forcing* forcing)
{
	unsigned int passes = (forcing->last_file - forcing->first_file + 1) / forcing->increment;
	if((forcing->last_file - forcing->first_file + 1)%forcing->increment != 0)	passes++;
	return passes;
}

//For flag = 3
unsigned int PassesDatabase(Forcing* forcing)
{
	return (unsigned int) ((forcing->last_file - forcing->first_file + 1)/(60.0*forcing->increment*forcing->file_time)) + 1;
}

//For flag = 7
unsigned int PassesRecurring(Forcing* forcing)
{
	return (forcing->last_file - forcing->first_file)/31536000 + 1;
}

//GetNextForcing ********************************************************************************
//Forcings (0 = none, 1 = .str, 2 = binary, 3 = database, 4 = .ustr, 5 = forcasting, 6 = .gz binary, 7 = recurring)

//For flag = 0,1,4
double NextForcingOther(Link** sys,unsigned int N,unsigned int* my_sys,unsigned int my_N,int* assignments,UnivVars* GlobalVars,Forcing* forcing,ConnData** db_connections,unsigned int** id_to_loc,unsigned int forcing_idx)
{
	return GlobalVars->maxtime;
}

//For flag = 2
double NextForcingBinaryFiles(Link** sys,unsigned int N,unsigned int* my_sys,unsigned int my_N,int* assignments,UnivVars* GlobalVars,Forcing* forcing,ConnData** db_connections,unsigned int** id_to_loc,unsigned int forcing_idx)
{
	unsigned int passes = forcing->passes, iteration = forcing->iteration;
	double maxtime;
	if(iteration == passes-1)	maxtime = GlobalVars->maxtime;
	else				maxtime = min(GlobalVars->maxtime,(iteration+1)*forcing->file_time*forcing->increment);
	int maxfileindex = (int) min((double) forcing->first_file+(iteration+1)*forcing->increment,(double) forcing->last_file);

	Create_Rain_Data_Par(sys,N,my_N,GlobalVars,my_sys,assignments,forcing->filename,forcing->first_file+iteration*forcing->increment,maxfileindex,iteration*forcing->file_time*forcing->increment,forcing->file_time,forcing,id_to_loc,forcing->increment+1,forcing_idx);

	return maxtime;
}

//For flag = 6
double NextForcingGZBinaryFiles(Link** sys,unsigned int N,unsigned int* my_sys,unsigned int my_N,int* assignments,UnivVars* GlobalVars,Forcing* forcing,ConnData** db_connections,unsigned int** id_to_loc,unsigned int forcing_idx)
{
	unsigned int passes = forcing->passes, iteration = forcing->iteration;
	double maxtime;
	if(iteration == passes-1)	maxtime = GlobalVars->maxtime;
	else				maxtime = min(GlobalVars->maxtime,(iteration+1)*forcing->file_time*forcing->increment);
	int maxfileindex = (int) min((double) forcing->first_file+(iteration+1)*forcing->increment,(double) forcing->last_file);

	Create_Rain_Data_GZ(sys,N,my_N,GlobalVars,my_sys,assignments,forcing->filename,forcing->first_file+iteration*forcing->increment,maxfileindex,iteration*forcing->file_time*forcing->increment,forcing->file_time,forcing,id_to_loc,forcing->increment+1,forcing_idx);

	return maxtime;
}

//For flag = 3
double NextForcingDatabase(Link** sys,unsigned int N,unsigned int* my_sys,unsigned int my_N,int* assignments,UnivVars* GlobalVars,Forcing* forcing,ConnData** db_connections,unsigned int** id_to_loc,unsigned int forcing_idx)
{
	unsigned int passes = forcing->passes, iteration = forcing->iteration,first_timestamp = 0;
	double maxtime;

	if( (int)(sys[my_sys[0]]->last_t*60 + .001) + forcing->raindb_start_time == (int)(forcing->first_file+iteration*forcing->file_time*60.0*forcing->increment+0.01) )
	{
		if(iteration == passes-1)	maxtime = GlobalVars->maxtime;	//!!!! Is this really needed? !!!!
		else				maxtime = min(GlobalVars->maxtime,(iteration+1)*forcing->file_time*forcing->increment);
		int maxfileindex = (int) min((double) forcing->first_file+(iteration+1)*60*forcing->file_time*forcing->increment,(double) forcing->last_file);

		first_timestamp = forcing->first_file + (unsigned int)(iteration*forcing->file_time*60.0*forcing->increment + 0.01);

		//If first_timestamp is off from the database increment, then read one previous timestamp
		if( (int)(first_timestamp - forcing->good_timestamp) % (int)(forcing->file_time*60+.01))
			first_timestamp -= (unsigned int)(forcing->file_time*60.0 + 0.01);
/*
if(forcing->flag == 3)
{
printf("Got real first = %u using first = %u good = %u\n",first_timestamp,forcing->first_file + (unsigned int)(iteration*forcing->file_time*60.0*forcing->increment + 0.01),forcing->good_timestamp);
printf("file_time = %f diff = %i mod = %i\n",forcing->file_time,(first_timestamp - forcing->good_timestamp),(int)(forcing->file_time*60+.01));
}
*/
		//Create_Rain_Database(sys,N,my_N,GlobalVars,my_sys,assignments,db_connections[ASYNCH_DB_LOC_FORCING_START + forcing_idx],forcing->first_file+iteration*forcing->file_time*60.0*forcing->increment,maxfileindex,forcing,id_to_loc,GlobalVars->maxtime,forcing_idx);
		Create_Rain_Database(sys,N,my_N,GlobalVars,my_sys,assignments,db_connections[ASYNCH_DB_LOC_FORCING_START + forcing_idx],first_timestamp,maxfileindex,forcing,id_to_loc,GlobalVars->maxtime,forcing_idx);
		(forcing->iteration)++;
	}
	else
		maxtime = min(GlobalVars->maxtime,(iteration)*forcing->file_time*forcing->increment);

/*
if(forcing->flag == 3)
{
int i,loc = find_link_by_idtoloc(318208,id_to_loc,N);
if(assignments[loc] == my_rank)
{
	printf("maxtime = %f real = %f ************************\n",maxtime,GlobalVars->maxtime);
	//printf("time = %e raindb_start_time = %u first = %u first used = %u\n",sys[my_sys[0]]->last_t,forcing->raindb_start_time,(int)(forcing->first_file+iteration*forcing->file_time*60.0*val+0.01),first_timestamp);
	//printf("time = %e raindb_start_time = %u first used = %u\n",sys[my_sys[0]]->last_t,forcing->raindb_start_time,first_timestamp);
	printf("current time = %e change time = %e current value = %e\n",sys[loc]->last_t,sys[loc]->forcing_change_times[2],sys[loc]->forcing_values[2]);
	for(i=0;i<sys[loc]->forcing_buff[2]->n_times-1;i++)
	{
		printf("%f %f\n",sys[loc]->forcing_buff[2]->rainfall[i][0],sys[loc]->forcing_buff[2]->rainfall[i][1]);
	}
	printf("\n");
}
MPI_Barrier(MPI_COMM_WORLD);
}
*/

/*
if(forcing->flag == 5)
{
int i,loc = find_link_by_idtoloc(318208,id_to_loc,N);
if(assignments[loc] == my_rank)
{
	printf("maxtime = %f real = %f ************************\n",maxtime,GlobalVars->maxtime);
	//printf("time = %e raindb_start_time = %u first = %u first used = %u\n",sys[my_sys[0]]->last_t,forcing->raindb_start_time,(int)(forcing->first_file+iteration*forcing->file_time*60.0*val+0.01),first_timestamp);
	//printf("time = %e raindb_start_time = %u first used = %u\n",sys[my_sys[0]]->last_t,forcing->raindb_start_time,first_timestamp);
	printf("current time = %e change time = %e current value = %e\n",sys[loc]->last_t,sys[loc]->forcing_change_times[0],sys[loc]->forcing_values[0]);
	for(i=0;i<sys[loc]->forcing_buff[0]->n_times-1;i++)
	{
		printf("%f %f\n",sys[loc]->forcing_buff[0]->rainfall[i][0],sys[loc]->forcing_buff[0]->rainfall[i][1]);
	}
	printf("\n");
}
MPI_Barrier(MPI_COMM_WORLD);
}
*/

	return maxtime;
}

//For flag = 7
double NextForcingRecurring(Link** sys,unsigned int N,unsigned int* my_sys,unsigned int my_N,int* assignments,UnivVars* GlobalVars,Forcing* forcing,ConnData** db_connections,unsigned int** id_to_loc,unsigned int forcing_idx)
{
	struct tm next_time,*first_time;
	time_t casted_first_file = (time_t) forcing->first_file;

	//Set current_epoch
	if(forcing->iteration)
	{
		//Get the year for the initial timestamp
		first_time = gmtime(&casted_first_file);
		next_time.tm_sec = 0;
		next_time.tm_min = 0;
		next_time.tm_hour = 0;
		next_time.tm_mday = 1;
		next_time.tm_mon = 0;
		next_time.tm_year = first_time->tm_year + forcing->iteration;
		next_time.tm_isdst = 0;
	}
	else
	{
		first_time = gmtime(&casted_first_file);
		copy_tm(first_time,&next_time);	//Yeah, this is lazy. But it's only done once...
	}

	return CreateForcing_Monthly(sys,my_N,my_sys,GlobalVars,forcing->GlobalForcing,forcing_idx,&next_time,forcing->first_file,forcing->last_file,sys[my_sys[0]]->last_t);
}


