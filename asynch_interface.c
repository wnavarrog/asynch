#include "asynch_interface.h"


//Initializes the asynch solver object. Also initializes MPI, if it's not set.
asynchsolver* Asynch_Init(MPI_Comm comm)
{
	unsigned int i;

	asynchsolver* asynch = (asynchsolver*) malloc(sizeof(asynchsolver));
	asynch->comm = comm;
	int init_flag;
	if(comm != MPI_COMM_WORLD)	printf("Warning: asynchsolver object my not work fully with in a comm other than MPI_COMM_WORLD.\n");

	//Initialize MPI stuff
	MPI_Initialized(&init_flag);
	if(!init_flag)	MPI_Init(NULL,NULL);
	MPI_Comm_rank(comm,&(asynch->my_rank));
	MPI_Comm_size(comm,&(asynch->np));

	//Sets the global variables
	np = asynch->np;
	my_rank = asynch->my_rank;

	//Initialize asynchsolver members
	asynch->outputfile = NULL;
	asynch->getting = NULL;
	asynch->my_data = NULL;
	asynch->peakfile = NULL;
	asynch->peakfilename = NULL;
	asynch->custom_model = NULL;
	asynch->ExternalInterface = NULL;
	for(i=0;i<ASYNCH_MAX_DB_CONNECTIONS;i++)	asynch->db_connections[i] = NULL;
	for(i=0;i<ASYNCH_MAX_DB_CONNECTIONS - ASYNCH_DB_LOC_FORCING_START;i++)	asynch->forcings[i] = InitializeForcings();

	//Initialize data types
	asynch->dt_info = Init_DataTypes();

	return asynch;
}

//Creates and loads a custom model
//Return 1 if there is an error, 0 if everything is ok
int Asynch_Custom_Model(asynchsolver* asynch,void (*SetParamSizes)(UnivVars*,void*),void (*Convert)(VEC*,unsigned int,void*),void (*Routines)(Link*,unsigned int,unsigned int,unsigned short int,void*),
	void (*Precalculations)(Link*,VEC*,VEC*,IVEC*,unsigned int,unsigned int,unsigned short int,unsigned int,void*),unsigned int (*InitializeEqs)(VEC*,VEC*,IVEC*,QVSData*,unsigned short int,VEC*,unsigned int,unsigned int,unsigned int,void*))
{
	asynch->custom_model = (model*) malloc(sizeof(model));
/*
	// !!!! Commented for now. I might actually want functions set to NULL for interfaces... !!!!
	if(!SetParamSizes || !Convert || !Routines || !Precalculations || !InitializeEqs)
	{
		if(my_rank == 0)
			printf("[%i]: Error creating custom model: Some routines to initialize the model are NULL.\n",my_rank);
		return 1;
	}
*/
	asynch->custom_model->SetParamSizes = SetParamSizes;
	asynch->custom_model->Convert = Convert;
	asynch->custom_model->Routines = Routines;
	asynch->custom_model->Precalculations = Precalculations;
	asynch->custom_model->InitializeEqs = InitializeEqs;

	return 0;
}


//Reads a .gbl file.
void Asynch_Parse_GBL(asynchsolver* asynch,char* gbl_filename)
{
	//Read in .gbl file
	asynch->GlobalVars = Read_Global_Data(gbl_filename,&(asynch->GlobalErrors),(Forcing**) &(asynch->forcings),asynch->db_connections,asynch->rkdfilename,asynch->custom_model,asynch->ExternalInterface);
	if(!asynch->GlobalVars)
	{
		printf("[%i]: An error occurred reading the .gbl file. See above messages for details.\n",my_rank);
		MPI_Abort(asynch->comm,1);
	}
}

//Loads the topology of the river network.
void Asynch_Load_System(asynchsolver* asynch)
{
	unsigned int i;

	//Read in remaining data from files
	asynch->sys = Create_River_System_parallel(asynch->rkdfilename,&(asynch->N),&(asynch->my_sys),&(asynch->my_N),&(asynch->my_max_nodes),&(asynch->my_data),&(asynch->assignments),&(asynch->getting),
		&(asynch->AllMethods),&(asynch->nummethods),asynch->GlobalVars,asynch->GlobalErrors,&(asynch->save_list),&(asynch->my_save_size),&(asynch->save_size),
		&(asynch->peaksave_size),&(asynch->id_to_loc),asynch->forcings,asynch->db_connections,&(asynch->workspace),asynch->custom_model,asynch->ExternalInterface);
	if(asynch->sys == NULL)
	{
		printf("[%i]: An error occurred setting up the river network. See above messages for details.\n",my_rank);
		MPI_Abort(asynch->comm,1);
	}

	//Put together the output filename string
	char filename[asynch->GlobalVars->string_size];
	if(asynch->GlobalVars->print_par_flag == 1)
	{
		for(i=0;i<asynch->GlobalVars->global_params->dim;i++)
			sprintf(filename,"_%.4e",asynch->GlobalVars->global_params->ve[i]);
	}
/*
	//Setup temporary output data file
	asynch->outputfile = PrepareTempFiles(asynch->sys,asynch->N,asynch->assignments,asynch->GlobalVars,asynch->save_list,asynch->save_size,asynch->my_save_size,NULL,asynch->id_to_loc);

	//Setup peakflow files
	if(PreparePeakFlowFiles(asynch->GlobalVars,asynch->peaksave_size))
		MPI_Abort(MPI_COMM_WORLD,1);
*/
	//Everything should be initialized...
	MPI_Barrier(asynch->comm);
}

//Trash an asynchsolver object
void Asynch_Free(asynchsolver* asynch)
{
	int i,finalized_flag;

	Free_DataTypes(&(asynch->dt_info));
	TransData_Free(asynch->my_data);
	for(i=0;i<ASYNCH_MAX_DB_CONNECTIONS;i++)	ConnData_Free(asynch->db_connections[i]);
	Destroy_ErrorData(asynch->GlobalErrors);
	Destroy_Workspace(asynch->workspace,asynch->GlobalVars->max_s,asynch->GlobalVars->max_parents);
	free(asynch->workspace);
	free(asynch->getting);
	if(asynch->outputfile)	fclose(asynch->outputfile);
	if(asynch->peakfilename)	free(asynch->peakfilename);

	for(i=0;i<asynch->N;i++)	Destroy_Link(asynch->sys[i],asynch->GlobalVars->iter_limit,asynch->rkdfilename[0] != '\0',asynch->forcings,asynch->GlobalVars);
	for(i=0;i<ASYNCH_MAX_DB_CONNECTIONS - ASYNCH_DB_LOC_FORCING_START;i++)	FreeForcing(&(asynch->forcings[i]));
	free(asynch->sys);
	free(asynch->my_sys);
	free(asynch->assignments);
	for(i=0;i<asynch->nummethods;i++)	Destroy_RKMethod(asynch->AllMethods[i]);
	free(asynch->AllMethods);
	if(asynch->save_list != NULL)	free(asynch->save_list);
	for(i=0;i<asynch->N;i++)
		free(asynch->id_to_loc[i]);
	free(asynch->id_to_loc);
	Destroy_UnivVars(asynch->GlobalVars);
	if(asynch->custom_model)	free(asynch->custom_model);

	free(asynch);

	//Finalize MPI
	MPI_Finalized(&finalized_flag);
	if(!finalized_flag)	MPI_Finalize();
}


void Asynch_Advance(asynchsolver* asynch,short int print_flag)
{
	if(print_flag && !OutputsSet(asynch->GlobalVars))
	{
		printf("[%i]: Warning: Solver advance requested with data output enabled, but not all outputs are initialized. Continuing solver without outputing data.\n",my_rank);
		print_flag = 0;
	}

	AsynchSolver(asynch->sys,asynch->N,asynch->my_sys,asynch->my_N,asynch->my_max_nodes,asynch->GlobalVars,asynch->assignments,
		asynch->id_to_loc,asynch->workspace,asynch->forcings,asynch->db_connections,asynch->my_data,print_flag,asynch->outputfile);
}


//Returns 0 if snapshot was made, 1 if an error was encountered, -1 if no snapshot was taken
int Asynch_Take_System_Snapshot(asynchsolver* asynch,char* name)
{
	//DataDump2(asynch->sys,asynch->N,asynch->assignments,asynch->GlobalVars,name);
	if(!(asynch->GlobalVars->output_data->CreateSnapShot))	return -1;
	return asynch->GlobalVars->output_data->CreateSnapShot(asynch->sys,asynch->N,asynch->assignments,asynch->GlobalVars,name,asynch->db_connections[ASYNCH_DB_LOC_SNAPSHOT_OUTPUT]);
}

void Asynch_Set_Database_Connection(asynchsolver* asynch,char* database_info,unsigned int conn_idx)
{
	if(asynch->db_connections[conn_idx])
	{
		ConnData_Free(asynch->db_connections[conn_idx]);
		asynch->db_connections[conn_idx] = NULL;
	}
	asynch->db_connections[conn_idx] = CreateConnData(database_info);
}

double Asynch_Get_Total_Simulation_Time(asynchsolver* asynch)
{
	return asynch->GlobalVars->maxtime;
}

void Asynch_Set_Total_Simulation_Time(asynchsolver* asynch,double new_time)
{
	asynch->GlobalVars->maxtime = new_time;
}

unsigned int Asynch_Get_Last_Rainfall_Timestamp(asynchsolver* asynch,unsigned int forcing_idx)
{
	return asynch->forcings[forcing_idx]->first_file + (unsigned int) (60.0 * asynch->forcings[forcing_idx]->maxtime);
	//return asynch->GlobalVars->first_file + (unsigned int) (60.0 * asynch->GlobalVars->maxtime);
}

void Asynch_Set_First_Rainfall_Timestamp(asynchsolver* asynch,unsigned int epoch_timestamp,unsigned int forcing_idx)
{
	asynch->forcings[forcing_idx]->first_file = epoch_timestamp;
	//asynch->GlobalVars->first_file = epoch_timestamp;
}

void Asynch_Set_Last_Rainfall_Timestamp(asynchsolver* asynch,unsigned int epoch_timestamp,unsigned int forcing_idx)
{
	asynch->forcings[forcing_idx]->last_file = epoch_timestamp;
	//asynch->GlobalVars->last_file = epoch_timestamp;
}

unsigned int Asynch_Get_First_Rainfall_Timestamp(asynchsolver* asynch,unsigned int forcing_idx)
{
	return asynch->forcings[forcing_idx]->first_file;
	//return asynch->GlobalVars->first_file;
}

void Asynch_Set_RainDB_Starttime(asynchsolver* asynch,unsigned int epoch_timestamp,unsigned int forcing_idx)
{
	asynch->forcings[forcing_idx]->raindb_start_time = epoch_timestamp;
	//asynch->GlobalVars->raindb_start_time = epoch_timestamp;
}

//Returns 0 if everything is good. Returns 1 if init_flag does not support a timestamp.
int Asynch_Set_Init_Timestamp(asynchsolver* asynch,unsigned int epoch_timestamp)
{
	if(asynch->GlobalVars->init_flag != 3)	return 1;
	asynch->GlobalVars->init_timestamp = epoch_timestamp;
	return 0;
}

unsigned int Asynch_Get_Init_Timestamp(asynchsolver* asynch)
{
	return asynch->GlobalVars->init_timestamp;
}

void Asynch_Set_Init_File(asynchsolver* asynch,char* filename)
{
	sprintf(asynch->GlobalVars->init_filename,"%s",filename);

	int length;
	for(length=0;length<256;length++)
		if(filename[length] == '\0')	break;

	if(length < 3 || length == 256)
	{
		if(my_rank == 0)	printf("Error: Bad init filename: %s.\n",filename);
		MPI_Abort(asynch->comm,1);
	}

	//Check what type of file
	if(filename[length-1] == 'c')
	{
		if(filename[length-2] == 'e' && filename[length-3] == 'r' && filename[length-4] == '.')
		{
			asynch->GlobalVars->init_flag = 2;
			return;
		}
	}

	if(filename[length-1] == 'i' && filename[length-2] == 'n' && filename[length-3] == 'i')
	{
		if(filename[length-4] == '.')
		{
			asynch->GlobalVars->init_flag = 0;
			return;
		}

		if(filename[length-4] == 'u' && filename[length-5] == '.')
		{
			asynch->GlobalVars->init_flag = 1;
			return;
		}
	}

	if(my_rank == 0)	printf("Error: Bad init filename: %s.\n",filename);
	MPI_Abort(asynch->comm,1);
}


void Asynch_Prepare_Output(asynchsolver* asynch)
{
	if(asynch->GlobalVars->output_data->PrepareOutput)	asynch->GlobalVars->output_data->PrepareOutput(asynch->GlobalVars,asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT]);
}


void Asynch_Prepare_Temp_Files(asynchsolver* asynch)
{
	if(asynch->GlobalVars->output_data->PrepareTempOutput)
		asynch->outputfile = asynch->GlobalVars->output_data->PrepareTempOutput(asynch->sys,asynch->N,asynch->assignments,asynch->GlobalVars,asynch->save_list,asynch->save_size,asynch->my_save_size,NULL,asynch->id_to_loc);
	else
		asynch->outputfile = NULL;
}

//Writes the current state to the temp files, if the current state is at a print time.
//Returns 0 if ok, 1 if no temp file is open.
int Asynch_Write_Current_Step(asynchsolver* asynch)
{
	Link **sys = asynch->sys,*current;
	unsigned int i,N = asynch->N,loc,save_size = asynch->save_size;
	int *assignments = asynch->assignments;
	double time_diff;

	if(asynch->my_save_size && !(asynch->outputfile))
	{
		printf("[%i]: Error writting step. No temporary file is open.\n",my_rank);
		return 1;
	}

	for(i=0;i<save_size;i++)
	{
		loc = find_link_by_idtoloc(asynch->save_list[i],asynch->id_to_loc,N);
		current = sys[loc];
		if(assignments[loc] == my_rank)
		{
			time_diff = fabs(current->next_save - current->last_t);
			if( time_diff/current->next_save < 1e-12 || ( (fabs(current->next_save) < 1e-12) ? (time_diff < 1e-12) : 0 ) )
			{
//if(current->ID == 318208)
//printf("Init step %e\n",current->list->head->y_approx->ve[0]);
				//WriteStep(current->last_t,current->list->head->y_approx,asynch->GlobalVars,current->params,current->iparams,current->state,asynch->outputfile,current->output_user,&(current->pos));
				WriteStep(current->last_t,current->list->head->y_approx,asynch->GlobalVars,current->params,current->iparams,current->state,asynch->outputfile,current->output_user,&(current->pos_offset));	//!!!! Should be tail? !!!!
				current->next_save += current->print_time;
				(current->disk_iterations)++;
			}
		}
	}

	return 0;
}

void Asynch_Prepare_Peakflow_Output(asynchsolver* asynch)
{
	if(asynch->GlobalVars->output_data->PreparePeakflowOutput)	asynch->GlobalVars->output_data->PreparePeakflowOutput(asynch->GlobalVars,asynch->peaksave_size);
}

//Return 0 means ok, -1 means no data to output
int Asynch_Create_Output(asynchsolver* asynch,char* additional_out)
{
	if(asynch->GlobalVars->output_data->CreateOutput && asynch->GlobalVars->hydrosave_flag)
		return asynch->GlobalVars->output_data->CreateOutput(asynch->sys,asynch->GlobalVars,asynch->N,asynch->save_list,asynch->save_size,asynch->my_save_size,asynch->id_to_loc,asynch->assignments,NULL,additional_out,asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT],&(asynch->outputfile));
	return -1;
}

//Return 0 means ok, -1 means no peakflows to output
int Asynch_Create_Peakflows_Output(asynchsolver* asynch)
{
	if(asynch->GlobalVars->output_data->CreatePeakflowOutput && asynch->GlobalVars->peaksave_flag)
		return asynch->GlobalVars->output_data->CreatePeakflowOutput(asynch->sys,asynch->GlobalVars,asynch->my_sys,asynch->my_N,asynch->peaksave_size,asynch->db_connections[ASYNCH_DB_LOC_PEAK_OUTPUT]);
	return -1;
}

//Return 0 means ok, 1 means no peakflows to output
int Asynch_Set_Peakflow_Output_Name(asynchsolver* asynch,char* peakflowname)
{
	if(asynch->GlobalVars->peaks_loc_filename)
	{
		sprintf(asynch->GlobalVars->peaks_loc_filename,peakflowname);
		return 0;
	}
	else
		return 1;
}

//Return 0 means ok, 1 means no peakflows to output
int Asynch_Get_Peakflow_Output_Name(asynchsolver* asynch,char* peakflowname)
{
	if(asynch->GlobalVars->peaks_loc_filename)
	{
		sprintf(peakflowname,asynch->GlobalVars->peaks_loc_filename);
		return 0;
	}
	else
		return 1;
}


unsigned int Asynch_Get_Number_Links(asynchsolver* asynch)
{
	if(!asynch)	return 0;
	return asynch->N;
}

unsigned int Asynch_Get_Local_Number_Links(asynchsolver* asynch)
{
	if(!asynch)	return 0;
	return asynch->my_N;
}
/*
int Asynch_Upload_Hydrographs_Database(asynchsolver* asynch)
{
	return UploadHydrosDB(asynch->sys,asynch->GlobalVars,asynch->N,asynch->save_list,asynch->save_size,asynch->my_save_size,asynch->id_to_loc,asynch->assignments,NULL,asynch->db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT]);
}
*/

int Asynch_Set_Temp_Files(asynchsolver* asynch,double set_time,void* set_value,unsigned int output_idx)
{
	return SetTempFiles(set_time,set_value,asynch->GlobalVars->output_types[output_idx],output_idx,asynch->sys,asynch->N,asynch->outputfile,asynch->GlobalVars,asynch->my_save_size,asynch->id_to_loc,asynch->dt_info);
	//return SetTempFiles(set_time,asynch->sys,asynch->N,asynch->outputfile,asynch->GlobalVars,asynch->my_save_size,asynch->id_to_loc);
}

int Asynch_Reset_Temp_Files(asynchsolver* asynch,double set_time)
{
	return ResetTempFiles(set_time,asynch->sys,asynch->N,asynch->outputfile,asynch->GlobalVars,asynch->my_save_size,asynch->id_to_loc);
}

int Asynch_Set_Output(asynchsolver* asynch,char* name,short int data_type,void (*func)(double,VEC*,VEC*,VEC*,IVEC*,int,void*),int* used_states,int num_states)
{
	unsigned int i,j,idx;
	UnivVars* GlobalVars = asynch->GlobalVars;

	//Find index
	for(i=0;i<asynch->GlobalVars->num_print;i++)
	{
		if(strcmp(name,asynch->GlobalVars->output_names[i]) == 0)
		{
			idx = i;
			break;
		}
	}

	if(i == asynch->GlobalVars->num_print)
	{
		printf("[%i]: Output %s not set.\n",my_rank,name);
		return 0;
	}

	//Add new output
	switch(data_type)
	{
		case ASYNCH_DOUBLE:
			asynch->GlobalVars->output_types[idx] = ASYNCH_DOUBLE;
			asynch->GlobalVars->outputs_d[idx] = (double (*)(double,VEC*,VEC*,VEC*,IVEC*,int,void*)) *func;
			break;
		case ASYNCH_INT:
			asynch->GlobalVars->output_types[idx] = ASYNCH_INT;
			asynch->GlobalVars->outputs_i[idx] = (int (*)(double,VEC*,VEC*,VEC*,IVEC*,int,void*)) *func;
			break;
		default:
			printf("[%i]: Error: Cannot set output. Bad function type (%hi).\n",my_rank,data_type);
			MPI_Abort(MPI_COMM_WORLD,1);
	}

	asynch->GlobalVars->output_sizes[idx] = GetByteSize(data_type);
	GetSpecifier(asynch->GlobalVars->output_specifiers[idx],data_type);

	//Check if anything should be added to dense_indices from used_states
	unsigned int num_to_add = 0;
	unsigned int* states_to_add = (unsigned int*) malloc(num_states*sizeof(unsigned int));
	for(i=0;i<num_states;i++)
	{
		for(j=0;j<GlobalVars->num_dense;j++)
		{
			if(used_states[i] == GlobalVars->dense_indices[j])
			{
				states_to_add[num_to_add++] = used_states[i];
				break;
			}
		}
	}

	if(num_to_add)
	{
		GlobalVars->dense_indices = (unsigned int*) realloc(GlobalVars->dense_indices,(GlobalVars->num_dense + num_to_add)*sizeof(unsigned int));
		for(i=0;i<num_to_add;i++)
			GlobalVars->dense_indices[i+GlobalVars->num_dense] = states_to_add[i];
		GlobalVars->num_dense += num_to_add;
		merge_sort_1D(GlobalVars->dense_indices,GlobalVars->num_dense);
	}

	free(states_to_add);

	return 1;
}

//Returns 1 if output function set successfully, 0 if there was a problem
int Asynch_Set_Peakflow_Output(asynchsolver* asynch,char* name,void (*func)(unsigned int,double,VEC*,VEC*,VEC*,double,unsigned int,void*,char*))
{
	if(!asynch->GlobalVars->peakflow_function_name)	asynch->GlobalVars->peakflow_function_name = (char*) malloc(asynch->GlobalVars->string_size*sizeof(char));

	if(strcmp(name,asynch->GlobalVars->peakflow_function_name))
	{
		printf("[%i]: Error setting peakflow output function to %s. Function %s was previously specified.\n",my_rank,asynch->GlobalVars->peakflow_function_name,name);
		return 0;
	}
	
	asynch->GlobalVars->peakflow_output = func;

	return 1;
}

//Returns the link id of the link at my_sys[location]
unsigned int Asynch_Get_Local_LinkID(asynchsolver* asynch,unsigned int location)
{
	return asynch->sys[asynch->my_sys[location]]->ID;
}

//Returns 1 if output name is set,
//0 if output name is not set,
//-1 if output name is not present
int Asynch_Check_Output(asynchsolver* asynch,char* name)
{
	unsigned int i;

	for(i=0;i<asynch->GlobalVars->num_print;i++)
	{
		if(strcmp(name,asynch->GlobalVars->output_names[i]) == 0)
		{
			if(asynch->GlobalVars->output_types[i] == ASYNCH_BAD_TYPE)	return 0;
			else								return 1;
		}
	}

	return -1;
}

//Returns 1 if using the peakflow output name is set,
//0 if the peakflow output name is not specified,
//-1 if the peakflow output name is not present
int Asynch_Check_Peakflow_Output(asynchsolver* asynch,char* name)
{
	if(!asynch->GlobalVars->peakflow_function_name || !name || strcmp(name,asynch->GlobalVars->peakflow_function_name))	return -1;
	return PeakflowOutputsSet(asynch->GlobalVars);
}


int Asynch_Delete_Temporary_Files(asynchsolver* asynch)
{
	int ret_val = RemoveTemporaryFiles(asynch->GlobalVars,asynch->my_save_size,NULL);
	//if(ret_val == 1)	printf("[%i]: Error deleting temp file. File does not exist.\n");
	if(ret_val == 2)	printf("[%i]: Error deleting temp file.\n",my_rank);
	return ret_val;
}

//Return 0 if ok, 1 if error
//!!!! This should set the current forcing to something. But what if maxtime is beyond the ceiling term? !!!!
int Asynch_Activate_Forcing(asynchsolver* asynch,unsigned int idx)
{
	unsigned int i,j,l,my_N = asynch->my_N;
	Link* current;

	if(idx >= asynch->GlobalVars->num_forcings)
	{
		printf("[%i]: Cannot activate forcing %u. Not enough forcings.\n",my_rank,idx);
		return 1;
	}

	asynch->forcings[idx]->active = 1;

/*
printf("Here\n");
getchar();
	//Set the forcing value at each link
	for(i=0;i<my_N;i++)
	{
		current = asynch->sys[asynch->my_sys[i]];
printf("rainfall buffer  t = %e\n",current->last_t);
printf("*********\n");
for(l=0;l<current->forcing_buff[idx]->n_times;l++)
{
printf("%e %e\n",current->forcing_buff[idx]->rainfall[l][0],current->forcing_buff[idx]->rainfall[l][1]);
getchar();
}
printf("*********\n");

		//Find the right index in rainfall
		for(l=0;l<current->forcing_buff[idx]->n_times;l++)
{
printf("l = %u/%u\n",l,current->forcing_buff[idx]->n_times);
printf("last t = %e\n",current->last_t);
printf("buffer time = %e\n",current->forcing_buff[idx]->rainfall[l][0]);
getchar();
			if( fabs(current->last_t - current->forcing_buff[idx]->rainfall[l][0]) < 1e-8 )	break;
}

		double forcing_buffer = current->forcing_buff[idx]->rainfall[l][1];
		current->forcing_values[idx] = forcing_buffer;

		//Find and set the new change in rainfall
		for(j=l+1;j<current->forcing_buff[idx]->n_times;j++)
		{
			if( fabs(current->forcing_buff[idx]->rainfall[j][1] - forcing_buffer) > 1e-8 )
			{
				current->forcing_change_times[idx] = current->forcing_buff[idx]->rainfall[j][0];
				break;
			}
		}
		if(j == current->forcing_buff[idx]->n_times)
			current->forcing_change_times[idx] = current->forcing_buff[idx]->rainfall[j-1][0];
	}
*/

	return 0;
}

//Return 0 if ok, 1 if error
int Asynch_Deactivate_Forcing(asynchsolver* asynch,unsigned int idx)
{
	unsigned int i,my_N = asynch->my_N,*my_sys = asynch->my_sys;
	Link** sys = asynch->sys;

	if(idx >= asynch->GlobalVars->num_forcings)
	{
		printf("[%i]: Cannot deactivate forcing %u. Not enough forcings.\n",my_rank,idx);
		return 1;
	}

	//Deactivate forcing
	asynch->forcings[idx]->active = 0;

	//Clear forcing values from links
	for(i=0;i<my_N;i++)
	{
		sys[my_sys[i]]->forcing_values[idx] = 0.0;
		sys[my_sys[i]]->forcing_change_times[idx] = asynch->GlobalVars->maxtime + 1.0;
	}

	return 0;
}

//Sets information for a forcing.
//Returns 0 if everything is ok, 1 if there is an error.
int Asynch_Set_Forcing_State(asynchsolver* asynch,unsigned int idx,double t_0,unsigned int first_file,unsigned int last_file)
{
	if(asynch->GlobalVars->num_forcings <= idx)
	{
		printf("[%i]: Error setting forcing state. Bad forcing index %u / %u.\n",my_rank,idx,asynch->GlobalVars->num_forcings);
		return 1;
	}

	asynch->forcings[idx]->raindb_start_time = first_file;
	asynch->forcings[idx]->maxtime = t_0;
	asynch->forcings[idx]->iteration = 0;
	asynch->forcings[idx]->first_file = first_file;
	asynch->forcings[idx]->last_file = last_file;

	return 0;
}

//Resets the peakflow information at each link.
//The current time is used for the time to peak, and the last state is used for the values.
void Asynch_Reset_Peakflow_Data(asynchsolver* asynch)
{
	unsigned int i,my_N = asynch->my_N,*my_sys = asynch->my_sys;
	Link *current,**sys = asynch->sys;
	double t_0 = sys[my_sys[0]]->last_t;

	for(i=0;i<my_N;i++)
	{
		current = sys[my_sys[i]];
		current->peak_time = t_0;
		v_copy(current->list->tail->y_approx,current->peak_value);
	}
}

void Asynch_Set_System_State(asynchsolver* asynch,double t_0,VEC** backup)
{
	unsigned i,j,k,l;
	Link* current;

	//Unpack asynch
	Link** sys = asynch->sys;
	unsigned int N = asynch->N,problem_dim = asynch->GlobalVars->problem_dim,num_forcings = asynch->GlobalVars->num_forcings;
	//Forcing** forcings = asynch->forcings;
	UnivVars* GlobalVars = asynch->GlobalVars;

	//Reset temp file
	Asynch_Reset_Temp_Files(asynch,t_0);

	//Reset links
	for(i=0;i<N;i++)
	{
		current = sys[i];
		if(current->list != NULL)
		{
			while(current->current_iterations > 1)
			{
				Remove_Head_Node(current->list);
				(current->current_iterations)--;
			}
			current->list->head->t = t_0;
			current->last_t = t_0;
			current->steps_on_diff_proc = 1;
			current->iters_removed = 0;
			current->rejected = 0;
			if(current->numparents == 0)	current->ready = 1;
			else				current->ready = 0;
			for(j=0;j<problem_dim;j++)	current->list->head->y_approx->ve[j] = backup[i]->ve[j];
			v_copy(backup[i],current->list->head->y_approx);

/*
			//Reset the next_save time
			if(current->save_flag)
			{
				current->next_save = t_0;		//!!!! This forces the print times to match up with the assimilation times !!!!
				current->disk_iterations = 1;
			}
*/

			//Reset peak flow information
			current->peak_time = t_0;
			v_copy(current->list->head->y_approx,current->peak_value);

			//Reset current state
			if(current->state_check != NULL)
				current->state = current->state_check(current->list->head->y_approx,GlobalVars->global_params,current->params,current->qvs,current->dam);
			current->list->head->state = current->state;

			//Write initial state
			//if(current->save_flag)
			//	WriteStep(t_0,current->list->head->y_approx,asynch->GlobalVars,current->params,current->iparams,current->state,asynch->outputfile,current->output_user,&(current->pos));

			//Set forcings
			//!!!! This block was not here before I started toying with data assimilation stuff. Perhaps it causes problems with the forecasters... !!!!
			if(current->forcing_buff)
			{
				for(k=0;k<num_forcings;k++)
				{
					if(!(asynch->forcings[k]->flag))	continue;

					//Find the right index in rainfall
					for(l=0;l<current->forcing_buff[k]->n_times-1;l++)
						if( current->forcing_buff[k]->rainfall[l][0] <= t_0 && t_0 < current->forcing_buff[k]->rainfall[l+1][0] )	break;
					double rainfall_buffer = current->forcing_buff[k]->rainfall[l][1];
					current->forcing_values[k] = rainfall_buffer;
					current->forcing_indices[k] = l;

					//Find and set the new change in rainfall
					for(j=l+1;j<current->forcing_buff[k]->n_times;j++)
					{
						if( fabs(current->forcing_buff[k]->rainfall[j][1] - rainfall_buffer) > 1e-12 )
						{
							current->forcing_change_times[k] = current->forcing_buff[k]->rainfall[j][0];
							break;
						}
					}
					if(j == current->forcing_buff[k]->n_times)
						current->forcing_change_times[k] = current->forcing_buff[k]->rainfall[j-1][0];

					//Select new step size
					//current->h = InitialStepSize(current->last_t,current,GlobalVars,workspace);
				}
			}
		}
	}
}

//Allocates space for output_user at each link.
//Returns 0 if everything is good, 1 if space is already allocated, 2 if an error occurred.
int Asynch_Create_OutputUser_Data(asynchsolver* asynch,unsigned int data_size)
{
	if(!asynch)	return 2;

	unsigned int my_N = asynch->my_N,i,*my_sys = asynch->my_sys;
	Link** sys = asynch->sys;
	
	if(sys[my_sys[0]]->output_user)	return 1;

	for(i=0;i<my_N;i++)
		sys[my_sys[i]]->output_user = malloc(data_size);

	return 0;
}

//Deallocates space for output_user at each link.
//Returns 0 if everything is good, 1 if space is already deallocated, 2 if an error occurred.
int Asynch_Free_OutputUser_Data(asynchsolver* asynch)
{
	if(!asynch)	return 2;

	unsigned int my_N = asynch->my_N,i,*my_sys = asynch->my_sys;
	Link** sys = asynch->sys;
	
	if(sys[my_sys[0]]->output_user == NULL)	return 1;

	for(i=0;i<my_N;i++)
	{
		free(sys[my_sys[i]]->output_user);
		sys[my_sys[i]]->output_user = NULL;
	}

	return 0;
}

//Copies source into output_user for the link with location my_sys[location]
void Asynch_Copy_Local_OutputUser_Data(asynchsolver* asynch,unsigned int location,void* source,unsigned int size)
{
	//printf("Copying size %u %u %p\n",size,location,asynch->sys[asynch->my_sys[location]]->output_user);
	memcpy(asynch->sys[asynch->my_sys[location]]->output_user,source,size);
}

//Allocates space for output_user for the link with location my_sys[location]
void Asynch_Set_Size_Local_OutputUser_Data(asynchsolver* asynch,unsigned int location,unsigned int size)
{
	asynch->sys[asynch->my_sys[location]]->output_user = malloc(size);
	//printf("Setting to size %u %p\n",size,asynch->sys[asynch->my_sys[location]]->output_user);
}


