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

int my_rank;
int np;

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
	time_t start,stop;
	unsigned int *my_sys,my_N,my_max_nodes,N,nummethods,my_save_size,save_size,i,j,k,*save_list;
	int *assignments;
	char filename[128],rkdfilename[128];
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
	//PGconn* conn = NULL;
	ConnData* conninfo = CreateConnData();

	start = time(NULL);

	//Read in .gbl file
	GlobalVars = Read_Global_Data(argv[1],&GlobalErrors,conninfo,rkdfilename);
	if(GlobalVars == NULL)	return 1;

	// **** Use additional input arguments here ****
	printf("\n!!!! Warning: Overwritting parameters in .gbl file !!!!\n\n");
	strcpy(GlobalVars->identifier,argv[3]);
	GlobalVars->global_params->ve[8] = atof(argv[2]);
	// **** ****

	//Make sure type is correct
	if(GlobalVars->type != 200)
		printf("Warning: Type should be 200 for dirac. %u\n",GlobalVars->type);

	//Read in remaining data from files
	sys = Create_River_System_parallel(rkdfilename,&N,&my_sys,&my_N,&my_max_nodes,&my_data,&assignments,&getting,&AllMethods,&nummethods,GlobalVars,GlobalErrors,&save_list,&my_save_size,&save_size,&id_to_loc,conninfo,&workspace);
	if(sys == NULL)		return 1;

	//Put together the output filename string
	char outputfilename[128];
	sprintf(outputfilename,"%s%s",GlobalVars->results_folder,GlobalVars->identifier);

	if(GlobalVars->print_par_flag == 1)
	{
		for(i=0;i<GlobalVars->global_params->dim;i++)
		{
			sprintf(filename,"_%.4e",GlobalVars->global_params->ve[i]);
			strcat(outputfilename,filename);
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);
	stop = time(NULL);
	total_time = difftime(stop,start);
	if(my_rank == 0)	printf("Finished reading files. Total time: %f\n",total_time);

	//Initialize remaining data
	short int* done = calloc(my_max_nodes,sizeof(short int));
	short int parentsval;
	unsigned int alldone;
	RKSolutionNode *roottail;
	current = sys[my_sys[my_N-1]];
	unsigned int last_idx,curr_idx,around;
	unsigned int two_my_N = 2*my_N;
	short int rain_flag = GlobalVars->rain_flag;

	//Free some memory
	TransData_Free(my_data);

	double scale = 2.0 / N * GlobalVars->global_params->ve[9];
	double pi = 3.141592653589;
	double DL = GlobalVars->global_params->ve[8];
	double U = GlobalVars->global_params->ve[7];
	double time_step = 0.1;
	double t,integral;
	double* length_to_outlet = (double*) calloc(my_N,sizeof(double));
	unsigned int times = (unsigned int)((GlobalVars->maxtime - 0.0) / time_step) + 1;
	double* my_total = (double*) calloc(times,sizeof(double));
	double* total = (double*) malloc(times*sizeof(double));
	double S_0;
	double my_integral = 0.0;
	time_t sub_start;

	//Make sure everyone is good before getting down to it...
	printf("Process %i (%i total) is good to go with %i links.\n",my_rank,np,my_N);
	fflush(stdout);
	MPI_Barrier(MPI_COMM_WORLD);

	//Calculate length
	if(my_rank == 0)
	{
		printf("\nCalculating the length to the outlet\n");
		sub_start = time(NULL);
	}
	for(i=0;i<my_N;i++)
	{
		current = sys[my_sys[i]];
		while(current != NULL)
		{
			length_to_outlet[i] += current->params->ve[0];
			current = current->c;
		}
	}
	if(my_rank == 0)
	{
		stop = time(NULL);
		printf("Time to calculate length %f\n\n",difftime(stop,sub_start));
	}

	//Calculate response
	if(my_rank == 0)
	{
		printf("Evaluating basin response at outlet\n");
		sub_start = time(NULL);
	}
	my_total[0] = 0.0;
	t = time_step;
	for(i=1;i<times;i++)
	{
		for(j=0;j<my_N;j++)
			my_total[i] += length_to_outlet[j] * exp( -sq(length_to_outlet[j]-U*t*60.0)/(4.0*DL*t*60.0) );

		my_total[i] *= scale * pow(16.0*pi*DL*t*t*t*60.0*60.0*60.0,-0.5);
		t += time_step;
	}
	if(my_rank == 0)
	{
		stop = time(NULL);
		printf("Time to calculate the response %f\n\n",difftime(stop,sub_start));
	}

	//Calculate storage
	if(my_rank == 0)
	{
		printf("Calculating the total initial storage\n");
		sub_start = time(NULL);
	}
	for(i=0;i<times-1;i++)
		my_integral += time_step / 2.0 * (my_total[i] + my_total[i+1]);
	my_integral *= 60.0;	//Since time is in minutes
	if(my_rank == 0)
	{
		stop = time(NULL);
		printf("Time to calculate length %f\n\n",difftime(stop,sub_start));
	}

	//Reduce the data
	if(my_rank == 0)
	{
		printf("Reducing calculations\n");
		sub_start = time(NULL);
	}
	MPI_Reduce(my_total,total,times,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Reduce(&my_integral,&integral,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	if(my_rank == 0)
	{
		stop = time(NULL);
		printf("Time for reduction %f\n\n",difftime(stop,sub_start));
	}

	//Stop the clock
	MPI_Barrier(MPI_COMM_WORLD);
	stop = time(NULL);
	total_time += difftime(stop,start);

	//Output some data
	if(my_rank == 0)
	{
		printf("%i: The answer at time %.12f is %.12e\n",my_rank,GlobalVars->maxtime,total[times-1]);
		printf("%i: The total initial storage is %.12e\n",my_rank,integral);
	}

	//Write data to disk
	if(my_rank == 0)
	{
		sprintf(filename,"%s.dat",outputfilename);
		FILE* outputfile = fopen(filename,"w");

		if(outputfile == NULL)
		{
			printf("Error creating output file\n");
			return 1;
		}

		//fprintf(outputfile,"1\n200\n\n%u %u\n",sys[my_sys[0]]->ID,times);
		t = 0.0;
		for(i=0;i<times;i++)
		{
			fprintf(outputfile,"%.4f %.12e\n",t,total[i]);
			t += time_step;
		}

		fclose(outputfile);
	}

	//Cleanup
	free(my_total);
	free(total);
	free(length_to_outlet);
	ConnData_Free(conninfo);
	Destroy_ErrorData(GlobalErrors);
	Destroy_Workspace(workspace,GlobalVars->max_s,GlobalVars->max_parents);
	free(workspace);
	free(done);
	free(getting);
	MPI_Barrier(MPI_COMM_WORLD);	//Very important: make sure all processes have flushed their buffers before data is processed.

	fflush(stdout);
	if(my_rank == 0)	printf("The total time for the entire program: %f\n\n",total_time);

	//Last bit of cleanup
	for(i=0;i<N;i++)	Destroy_Link(sys[i],GlobalVars->iter_limit,rkdfilename[0] != '\0',GlobalVars);
	free(sys);
	free(my_sys);
	free(assignments);
	for(i=0;i<nummethods;i++)	Destroy_RKMethod(AllMethods[i]);
	free(AllMethods);
	free(save_list);
	for(i=0;i<N;i++)
		free(id_to_loc[i]);
	free(id_to_loc);

	//!!!! Make method for this !!!!
	if(GlobalVars->rvr_filename != NULL)	free(GlobalVars->rvr_filename);
	if(GlobalVars->prm_filename != NULL)	free(GlobalVars->prm_filename);
	if(GlobalVars->init_filename != NULL)	free(GlobalVars->init_filename);
	if(GlobalVars->rain_filename != NULL)	free(GlobalVars->rain_filename);
	v_free(GlobalVars->global_params);
	free(GlobalVars);

	MPI_Finalize();
	return 0;
}

