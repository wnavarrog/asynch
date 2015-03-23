#include "riversys.h"

//Reads in all a river system from several files. Used for a parallel implementation.
//char *_filename[]: Strings of filenames to load (river => .rvr, rk => .rkd, init => .ini, param => .prm, storm => .str, save => .sav).
//int* N (set by this method): Will be the number of links in the returned system.
//int** my_sys (set by this method): An array that will contain the location in the system of each link assigned to this process.
//int* my_N (set by this method): Will be the number of links assigned to this process.
//int* my_max_nodes (set by this method): Will contain the number of entries of my_sys (not all will be used).
//TransData* my_data (set by this method): Information about the links that each process will communicate information about will be stored in my_data.
//int** assignments (set by this method): Will be an array with N entries. assignments[i] will the process link sys[i] is assigned to.
//int** getting (set by this method): Will be an array of N entries. getting[i] = 1 if this process is receiving information about link sys[i], else 0.
//RKMethod*** AllMethods (set by this method): Will be a list of all available RK methods.
//int* nummethods (set by this method): Will be the number of methods in AllMethods.
//UnivVars* GlobalVars: Contains all the information that is shared by every link in the system.
//ErrorData* GlobalErrors: Contains information for the global error tolerances, if not using a .rkd file.
//int* my_save_size (set by this method): Will be the number of links assigned to this process for which information will be written to file.
//int* save_size (set by this method): Will be the number of links for which information will be written to file.
//int*** id_to_loc (set by this method): Will be an array with N rows and 2 columns, sorted by first col. First col is a link id and second is
//				the location of the id in sys.
//PGconn* conn:	Pointer for an SQL database. Will be NULL if no database is needed.
//TempStorage** workspace:	Will be set to memory pointers for temporary calculations.
//Returns an array (sys) of links created by reading in data from the files in the arguments.
Link** Create_River_System_parallel(char rk_filename[],unsigned int* N,unsigned int** my_sys,unsigned int* my_N,unsigned int* my_max_nodes,TransData** my_data,int** assignments,short int** getting,RKMethod*** AllMethods,unsigned int* nummethods,UnivVars* GlobalVars,ErrorData* GlobalErrors,unsigned int** save_list,unsigned int* my_save_size,unsigned int* save_size,unsigned int* peaksave_size,unsigned int*** id_to_loc,Forcing** forcings,ConnData** db_connections,TempStorage** workspace,model* custom_model,void* external)
{
	Link** system;
	Link *current,*prev;
	time_t start,stop;
	int ii,data1,crap2,inittype,min,max,db;
	unsigned int i,j,k,l,m,id,param_start,curr_loc,parentsval;
	unsigned int dim = GlobalVars->dim;
	int type = GlobalVars->type;
	Link** upstream_order;
	double t_0,crap1,trash,univ_forcing_change_time[ASYNCH_MAX_DB_CONNECTIONS - ASYNCH_DB_LOC_FORCING_START];
	//RainData* GlobalRain = NULL;
	//double* GlobalInit = NULL;
	VEC *y_0,*uniform_y_0;
	FILE* riverdata = NULL;
	if(GlobalVars->rvr_flag == 0)	riverdata = fopen(GlobalVars->rvr_filename,"r");
	FILE* rkdata;
	if(rk_filename[0] != '\0')	rkdata = fopen(rk_filename,"r");
	FILE* initdata = (GlobalVars->init_filename) ? fopen(GlobalVars->init_filename,"r") : NULL;
	FILE* paramdata = NULL;
	if(GlobalVars->prm_flag == 0)	paramdata = fopen(GlobalVars->prm_filename,"r");
	//FILE* forcingdata;
	PGresult *mainres,*res;
	VEC** db_init_buffer = NULL;
	//char* query = conninfo->query;

	if(GlobalVars->rvr_flag == 0 && riverdata == NULL)
	{
		printf("Error: file %s not found for .rvr file\n",GlobalVars->rvr_filename);
		*N = 0;
		return NULL;
	}
	if(rk_filename[0] != '\0' && rkdata == NULL)
	{
		printf("Error: file %s not found for .rkd file\n",rk_filename);
		*N = 0;
		return NULL;
	}
	if(GlobalVars->init_flag != 3 && initdata == NULL)
	{
		printf("Error: file %s not found for .ini or .uini file\n",GlobalVars->init_filename);
		*N = 0;
		return NULL;
	}
	if(GlobalVars->prm_flag == 0 && paramdata == NULL)
	{
		printf("Error: file %s not found for .prm file\n",GlobalVars->prm_filename);
		*N = 0;
		return NULL;
	}

	//if(rk_filename[0] != '\0')	fscanf(rkdata,"%i",&rktype);
	//else 				rktype = -1;
	if(initdata)	fscanf(initdata,"%i",&inittype);
	else		inittype = GlobalVars->type;
/*
	//If type is 1, then open the rainfall data
	if(GlobalVars->rain_flag == 1 || GlobalVars->rain_flag == 4)
	{
		//stormdata = fopen(storm_filename,"r");
		stormdata = fopen(GlobalVars->rain_filename,"r");
		if(stormdata == NULL)
		{
			printf("Error: file %s not found\n",GlobalVars->rain_filename);
			return NULL;
		}
	}
*/

	//Set the dimensions of the problem the .rkd and .ini files use
	//SetRKInitDim(rktype,inittype,&rkdim,&initdim);
	//rkdim = initdim = GlobalVars->problem_dim;		//!!!! Is this ok? !!!!

	//Get information from database, if needed
	unsigned int *db_link_id,*dbres_link_id,*dbres_parent,sizeres;
	float **db_params = NULL;
	db_link_id = dbres_link_id = dbres_parent = NULL;
	sizeres = 0;

	if(GlobalVars->rvr_flag == 1 || GlobalVars->prm_flag == 1)
	{
		if(my_rank == 0)	printf("\nTransferring data from database...\n");
		MPI_Barrier(MPI_COMM_WORLD);
		start = time(NULL);

		if(my_rank == 0)
		{
			if(GlobalVars->outletlink == 0)	//Grab entire network
			{
				if(GlobalVars->rvr_flag == 1 && GlobalVars->prm_flag == 1)
				{
					//Note: parent_id on the SQL database is like the child id here
					db = ConnectPGDB(db_connections[ASYNCH_DB_LOC_TOPO]);
					if(db)
					{
						printf("[%i]: Error connecting to the topology database.\n",my_rank);
						return NULL;
					}
					db = ConnectPGDB(db_connections[ASYNCH_DB_LOC_PARAMS]);
					if(db)
					{
						printf("[%i]: Error connecting to the parameter database.\n",my_rank);
						return NULL;
					}
					mainres = PQexec(db_connections[ASYNCH_DB_LOC_PARAMS]->conn,db_connections[ASYNCH_DB_LOC_PARAMS]->queries[0]);
					res = PQexec(db_connections[ASYNCH_DB_LOC_TOPO]->conn,db_connections[ASYNCH_DB_LOC_TOPO]->queries[1]);
					//mainres = PQexec(conninfo->conn,"select link_id,up_area,length,area from master_km where link_id < 620174 and link_id > 1 ORDER BY link_id");
					//res = PQexec(conninfo->conn,"WITH alllinks(id) AS (SELECT link_id FROM master_km) SELECT alllinks.id,master_km.link_id FROM master_km,alllinks WHERE alllinks.id > 1 AND master_km.parent_link = alllinks.id ORDER BY alllinks.id");
					CheckResError(res,"querying connectivity");
					CheckResError(mainres,"querying DEM data");
					param_start = 1;
					DisconnectPGDB(db_connections[ASYNCH_DB_LOC_TOPO]);
					DisconnectPGDB(db_connections[ASYNCH_DB_LOC_PARAMS]);
				}
				else if(GlobalVars->rvr_flag == 1)
				{
					//Note: parent_id on the SQL database is like the child id here
					db = ConnectPGDB(db_connections[ASYNCH_DB_LOC_TOPO]);
					if(db)
					{
						printf("[%i]: Error connecting to the topology database.\n",my_rank);
						return NULL;
					}
					mainres = PQexec(db_connections[ASYNCH_DB_LOC_TOPO]->conn,db_connections[ASYNCH_DB_LOC_TOPO]->queries[0]);
					res = PQexec(db_connections[ASYNCH_DB_LOC_TOPO]->conn,db_connections[ASYNCH_DB_LOC_TOPO]->queries[1]);
					//mainres = PQexec(conninfo->conn,"select link_id from master_km where link_id < 620174 and link_id > 1 ORDER BY link_id");
					//res = PQexec(conninfo->conn,"WITH alllinks(id) AS (SELECT link_id FROM master_km) SELECT alllinks.id,master_km.link_id FROM master_km,alllinks WHERE alllinks.id > 1 AND master_km.parent_link = alllinks.id ORDER BY alllinks.id");
					CheckResError(res,"querying connectivity");
					CheckResError(mainres,"querying DEM data");
					DisconnectPGDB(db_connections[ASYNCH_DB_LOC_TOPO]);
				}
				else if(GlobalVars->prm_flag == 1)
				{
					db = ConnectPGDB(db_connections[ASYNCH_DB_LOC_PARAMS]);
					if(db)
					{
						printf("[%i]: Error connecting to the parameter database.\n",my_rank);
						return NULL;
					}
					mainres = PQexec(db_connections[ASYNCH_DB_LOC_PARAMS]->conn,db_connections[ASYNCH_DB_LOC_PARAMS]->queries[0]);
					//mainres = PQexec(conninfo->conn,"select link_id,up_area,length,area from master_km where link_id < 620174 and link_id > 1 ORDER BY link_id");
					CheckResError(mainres,"querying DEM data");
					param_start = 1;
					DisconnectPGDB(db_connections[ASYNCH_DB_LOC_PARAMS]);
				}
			}
			else	//Grab a sub basin
			{
				//Currently, we assume both rvr and prm data is read from the database
				if( !(GlobalVars->rvr_flag == 1 && GlobalVars->prm_flag == 1) )
				{
					printf("Error: Database selection requires both rvr and prm data to come from database. %hu %hu\n",GlobalVars->rvr_flag,GlobalVars->prm_flag);
					return NULL;
				}

				db = ConnectPGDB(db_connections[ASYNCH_DB_LOC_TOPO]);
				if(db)
				{
					printf("[%i]: Error connecting to the topology database.\n",my_rank);
					return NULL;
				}
				db = ConnectPGDB(db_connections[ASYNCH_DB_LOC_PARAMS]);
				if(db)
				{
					printf("[%i]: Error connecting to the parameter database.\n",my_rank);
					return NULL;
				}

				//Make the queries
				//Be careful to not overload the database
				sprintf(db_connections[ASYNCH_DB_LOC_PARAMS]->query,db_connections[ASYNCH_DB_LOC_PARAMS]->queries[1],GlobalVars->outletlink);
				mainres = PQexec(db_connections[ASYNCH_DB_LOC_PARAMS]->conn,db_connections[ASYNCH_DB_LOC_PARAMS]->query);
				//sprintf(query,"SELECT nodeX.link_id,nodeX.up_area,nodeX.length,nodeX.area FROM master_km AS nodeX, master_km AS parentX WHERE (nodeX.left BETWEEN parentX.left AND parentX.right) AND parentX.link_id = %u ORDER BY link_id;",GlobalVars->outletlink);
				//mainres = PQexec(conninfo->conn,query);
				sprintf(db_connections[ASYNCH_DB_LOC_TOPO]->query,db_connections[ASYNCH_DB_LOC_TOPO]->queries[2],GlobalVars->outletlink);
				res = PQexec(db_connections[ASYNCH_DB_LOC_TOPO]->conn,db_connections[ASYNCH_DB_LOC_TOPO]->query);
/*
				sprintf(query,"WITH alllinks(id) AS (SELECT link_id FROM master_km) \
					SELECT alllinks.id,master_km.link_id FROM master_km,alllinks \
					WHERE (alllinks.id IN \
						(SELECT nodeX.link_id FROM master_km AS nodeX, master_km AS parentX \
						WHERE (nodeX.left BETWEEN parentX.left AND parentX.right) AND parentX.link_id = %u) ) \
					AND master_km.parent_link = alllinks.id ORDER BY alllinks.id",GlobalVars->outletlink);
*/
				//res = PQexec(conninfo->conn,query);
				CheckResError(res,"querying connectivity");
				CheckResError(mainres,"querying DEM data");
				param_start = 1;

				DisconnectPGDB(db_connections[ASYNCH_DB_LOC_TOPO]);
				DisconnectPGDB(db_connections[ASYNCH_DB_LOC_PARAMS]);
			}

			*N = PQntuples(mainres);

			//Allocate space and load buffers
			db_link_id = (unsigned int*) malloc(*N*sizeof(unsigned int));
			for(i=0;i<*N;i++)	db_link_id[i] = atoi(PQgetvalue(mainres,i,0));

			if(GlobalVars->rvr_flag == 1)
			{
				sizeres = PQntuples(res);
				dbres_link_id = (unsigned int*) malloc(sizeres*sizeof(unsigned int));
				dbres_parent = (unsigned int*) malloc(sizeres*sizeof(unsigned int));

				for(i=0;i<sizeres;i++)
				{
					dbres_link_id[i] = atoi(PQgetvalue(res,i,0));
					dbres_parent[i] = atoi(PQgetvalue(res,i,1));
				}

				PQclear(res);
			}
			if(GlobalVars->prm_flag == 1)
			{
				db_params = (float**) malloc(GlobalVars->disk_params * sizeof(float*));
				for(i=0;i<GlobalVars->disk_params;i++)
					db_params[i] = (float*) malloc(*N*sizeof(float));

				for(i=0;i<GlobalVars->disk_params;i++)
					for(j=0;j<*N;j++)
						db_params[i][j] = atof(PQgetvalue(mainres,j,param_start+i));
			}

			//Send sizes and data to other processes
			MPI_Bcast(&sizeres,1,MPI_INT,0,MPI_COMM_WORLD);
			MPI_Bcast(N,1,MPI_INT,0,MPI_COMM_WORLD);

			//Cleanup
			PQclear(mainres);
		}
		else	//Receiving DEM data from process 0
		{
			//Receive sizes
			MPI_Bcast(&sizeres,1,MPI_INT,0,MPI_COMM_WORLD);
			MPI_Bcast(N,1,MPI_INT,0,MPI_COMM_WORLD);

			//Allocate space
			db_link_id = (unsigned int*) malloc(*N*sizeof(unsigned int));
			dbres_link_id = (unsigned int*) malloc(sizeres*sizeof(unsigned int));
			dbres_parent = (unsigned int*) malloc(sizeres*sizeof(unsigned int));
			db_params = (float**) malloc(GlobalVars->disk_params * sizeof(float*));
			for(i=0;i<GlobalVars->disk_params;i++)
				db_params[i] = (float*) malloc(*N*sizeof(float));
		}

		//Broadcast the data
		MPI_Bcast(db_link_id,*N,MPI_INT,0,MPI_COMM_WORLD);
		if(GlobalVars->rvr_flag == 1)
		{
			MPI_Bcast(dbres_link_id,sizeres,MPI_INT,0,MPI_COMM_WORLD);
			MPI_Bcast(dbres_parent,sizeres,MPI_INT,0,MPI_COMM_WORLD);
		}
		if(GlobalVars->prm_flag == 1)
		{
			for(i=0;i<GlobalVars->disk_params;i++)
				MPI_Bcast(db_params[i],*N,MPI_FLOAT,0,MPI_COMM_WORLD);
		}

		MPI_Barrier(MPI_COMM_WORLD);
		stop = time(NULL);
		if(my_rank == 0)	printf("Time to receive data: %f\n",difftime(stop,start));
	}

	//Integrity Checks
	if(GlobalVars->rvr_flag == 0)	fscanf(riverdata,"%i",N);
	if(rk_filename[0] != '\0')	fscanf(rkdata,"%i",&j);
	else				j = *N;
	if(GlobalVars->prm_flag == 0)	fscanf(paramdata,"%i",&l);
	else				l = *N;
	if(GlobalVars->init_flag != 1 && initdata)	fscanf(initdata,"%i",&i);
	else 				i = l;
	if(initdata)	fscanf(initdata,"%lf",&t_0);
	else		t_0 = 0.0;
	GlobalVars->t_0 = t_0;
	if(rk_filename[0] != '\0')
	{
		fscanf(rkdata,"%i",&j);
		if(j != dim)
		{
			printf("[%i]: Error reading .rkd file: file has %u states; expected %u.\n",my_rank,j,dim);
			MPI_Abort(MPI_COMM_WORLD,1);
		}
	}

	if(my_rank == 0 && *N < (unsigned int)np)
		printf("\nWarning: using more processes (%i) than links (%u).\n",np,*N);

	//Adjust error and dimension
	//!!!! Enable AssimError, and upstream link calculations below !!!!
	if(GlobalVars->assim_flag == 1)
	{
		if(GlobalVars->type >= 300)	AssimError(*N,GlobalVars,GlobalErrors);
		else				printf("Notice: Skipping AssimError.\n");
		dim = GlobalVars->dim;
	}

	//Allocate some space
	system = (Link**) malloc(*N*sizeof(Link*));
	for(i=0;i<*N;i++)	system[i] = (Link*) malloc(sizeof(Link));
	unsigned int** loc_to_children = malloc(*N*sizeof(unsigned int*));	//This holds the ID of the children
	for(i=0;i<*N;i++)	loc_to_children[i] = malloc(10*sizeof(unsigned int));	//This assumes there are no more than 10 children
	upstream_order = malloc(*N*sizeof(Link*));
	y_0 = v_get(dim);
	uniform_y_0 = NULL;

	//Build all the RKMethods
	*nummethods = 4;
	*AllMethods = malloc(*nummethods * sizeof(RKMethod*));
	(*AllMethods)[0] = RKDense3_2();
	(*AllMethods)[1] = TheRKDense4_3();
	(*AllMethods)[2] = DOPRI5_dense();
	(*AllMethods)[3] = RadauIIA3_dense();
	GlobalVars->max_localorder = (*AllMethods)[0]->localorder;
	GlobalVars->max_s = (*AllMethods)[0]->s;
	for(i=1;i<*nummethods;i++)
	{
		GlobalVars->max_localorder = (GlobalVars->max_localorder < (*AllMethods)[i]->localorder) ? (*AllMethods)[i]->localorder : GlobalVars->max_localorder;
		GlobalVars->max_s = (GlobalVars->max_s > (*AllMethods)[i]->s) ? GlobalVars->max_s : (*AllMethods)[i]->s;
		//!!!! Note: Use a +1 for Radau solver? !!!!
	}

	//First pass: Build tree structure, save index, find root
	if(my_rank == 0)	printf("\nFirst read...\n");
	MPI_Barrier(MPI_COMM_WORLD);
	start = time(NULL);

	//Topology
	if(GlobalVars->rvr_flag == 0)
	{
		for(i=0;i<*N;i++)
		{
			fscanf(riverdata,"%u %i",&id,&data1);
			system[i]->location = i;
			system[i]->ID = id;

			//Set the parents
			system[i]->numparents = data1;
			system[i]->parents = (Link**) malloc(system[i]->numparents * sizeof(Link*));
			system[i]->c = NULL;
			GlobalVars->max_parents = (GlobalVars->max_parents > data1) ? GlobalVars->max_parents : data1;

			//Set a few other data
			system[i]->output_user = NULL;
			system[i]->f = NULL;
			system[i]->dam = 0;
			system[i]->qvs = NULL;
			system[i]->discont = NULL;
			system[i]->discont_start = 0;
			system[i]->discont_end = GlobalVars->discont_size-1;
			system[i]->discont_count = 0;
			system[i]->discont_send = NULL;
			system[i]->discont_order_send = NULL;
			system[i]->discont_send_count = 0;
			system[i]->res = 0;
			//system[i]->equations = NULL;

			//Set the information to find the parent links
			for(j=0;j<system[i]->numparents;j++)
			{
				fscanf(riverdata,"%i",&data1);
				loc_to_children[i][j] = data1;
			}
		}
	}
	else if(GlobalVars->rvr_flag == 1)
	{
		//Setup location and id info
		for(i=0;i<*N;i++)
		{
			id = db_link_id[i];

			system[i]->location = i;
			system[i]->ID = id;

			//Set a few other data
			system[i]->output_user = NULL;
			system[i]->numparents = 0;
			system[i]->parents = NULL;
			system[i]->c = NULL;
			system[i]->f = NULL;
			system[i]->dam = 0;
			system[i]->qvs = NULL;
			system[i]->discont = NULL;
			system[i]->discont_start = 0;
			system[i]->discont_end = GlobalVars->discont_size-1;
			system[i]->discont_count = 0;
			system[i]->discont_send = NULL;
			system[i]->discont_order_send = NULL;
			system[i]->discont_send_count = 0;
			system[i]->res = 0;
			//system[i]->equations = NULL;
		}
/*
		if(GlobalVars->outletlink == 0)	//Entire network
		{
			//Setup parents
			for(i=0;i<sizeres;i+=j)
			{
				//!!!! Assuming location = id - 2 !!!!
				id = dbres_link_id[i];
				current = system[id-2];

				for(j=0;i+j<sizeres;j++)
					if(dbres_link_id[i+j] != id)	break;

				//Set the parents
				current->numparents = j;
				current->parents = (Link**) malloc(current->numparents * sizeof(Link*));
				GlobalVars->max_parents = (GlobalVars->max_parents > current->numparents) ? GlobalVars->max_parents : current->numparents;

				//Set the information to find the child links !!!! Change to a query? !!!!
				for(k=0;k<j;k++)
					loc_to_children[id-2][k] = dbres_parent[i+k];
			}
		}
		else	//Sub basin
		{
*/
			//Setup parents
			curr_loc = 0;
			for(i=0;i<sizeres;i+=j)
			{
				//Select the next link
				current = system[curr_loc];
				id = current->ID;

				//Count the parents
				for(j=0;i+j<sizeres;j++)
					if(dbres_link_id[i+j] != id)	break;

				//Set the parents
				current->numparents = j;
				current->parents = (Link**) malloc(current->numparents * sizeof(Link*));
				GlobalVars->max_parents = (GlobalVars->max_parents > current->numparents) ? GlobalVars->max_parents : current->numparents;

				//Set the information to find the child links !!!! Change to a query? !!!!
				for(k=0;k<j;k++)
					loc_to_children[curr_loc][k] = dbres_parent[i+k];

				//Set the next location
				curr_loc++;
			}
//		}
	}
	else
	{
		printf("Error: Invalid topology flag (%hu).\n",GlobalVars->rvr_flag);
		return NULL;
	}

	//Parameters
	for(i=0;i<*N;i++)
	{
		if(GlobalVars->prm_flag == 0)
			fscanf(paramdata,"%u",&id);
		else
			if(db_link_id[i] != system[i]->ID)	printf("Something is wrong!!\n");

		if(GlobalVars->rvr_flag == 0 && GlobalVars->prm_flag == 0 && id != system[i]->ID)
		{
			printf("Error: Links must be listed in the same order in every data file %u %u %u\n",id,system[i]->ID,i);
			*N=0;
			return NULL;
		}

		//Read in the parameters
		system[i]->params = v_get(GlobalVars->params_size);
		system[i]->iparams = iv_get(GlobalVars->iparams_size);
		if(GlobalVars->prm_flag == 0)
		{
			for(j=0;j<GlobalVars->disk_params;j++)
				fscanf(paramdata,"%lf",&(system[i]->params->ve[j]));
		}
		else
		{
			for(j=0;j<GlobalVars->disk_params;j++)
				system[i]->params->ve[j] = db_params[j][i];
		}

		if(custom_model)	custom_model->Convert(system[i]->params,type,external);
		else	ConvertParams(system[i]->params,type,external);
	}

	//Cleanup database connection
	if(GlobalVars->rvr_flag == 1 || GlobalVars->prm_flag == 1)
	{
		if(db_link_id != NULL)		free(db_link_id);
		if(dbres_link_id != NULL)	free(dbres_link_id);
		if(dbres_parent != NULL)	free(dbres_parent);
		if(db_params != NULL)
		{
			for(i=0;i<GlobalVars->disk_params;i++)	free(db_params[i]);
			free(db_params);
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);
	stop = time(NULL);
	if(my_rank == 0)	printf("Time for first read: %f\n",difftime(stop,start));

	//Order the links by upstream area
	//Uses merge sort
	if(my_rank == 0)	printf("Sorting the system by upstream area...\n");
	MPI_Barrier(MPI_COMM_WORLD);
	start = time(NULL);

	for(i=0;i<*N;i++)
		upstream_order[i] = system[i];
	merge_sort(upstream_order,*N,GlobalVars->area_idx);

	MPI_Barrier(MPI_COMM_WORLD);
	stop = time(NULL);
	if(my_rank == 0)	printf("Time to sort by upstream area: %f\n",difftime(stop,start));

	//Make a list of ids and locations, sorted by id. Also make a list so id_to_loc can be modified easily.
	//Uses merge sort
	if(my_rank == 0)	printf("Sorting ids...\n");
	MPI_Barrier(MPI_COMM_WORLD);
	start = time(NULL);
	*id_to_loc = malloc(*N*sizeof(int*));
	for(i=0;i<*N;i++)	(*id_to_loc)[i] = malloc(2*sizeof(int));
	//int* id_tracker = malloc(*N*sizeof(int));	//More like loc_tracker. Takes locations into their corresponding row in id_to_loc.

	for(i=0;i<*N;i++)
	{
		(*id_to_loc)[i][0] = system[i]->ID;
		(*id_to_loc)[i][1] = i;
	}

	merge_sort_ids(*id_to_loc,*N);
	//for(i=0;i<*N;i++)	id_tracker[(*id_to_loc)[i][1]] = i;
	MPI_Barrier(MPI_COMM_WORLD);
	stop = time(NULL);
	if(my_rank == 0)	printf("Time to sort ids: %f\n",difftime(stop,start));

	//Setup Parent Data
	if(my_rank == 0)	printf("Setting up parent info...\n");
	MPI_Barrier(MPI_COMM_WORLD);
	start = time(NULL);

	for(i=0;i<*N;i++)
	{
		for(j=0;j<system[i]->numparents;j++)
		{
			k = *N/2;
			min = 0;
			max = *N;
			while((*id_to_loc)[k][0] != loc_to_children[i][j])
			{
				if((*id_to_loc)[k][0] > loc_to_children[i][j])	max = k;
				else						min = k;
				k = (max + min)/2;
			}

			system[i]->parents[j] = system[(*id_to_loc)[k][1]];
			system[(*id_to_loc)[k][1]]->c = system[i];	//!!!! New for Morgan's problem !!!!
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	stop = time(NULL);
	if(my_rank == 0)	printf("Time to setup parents: %f\n",difftime(stop,start));

	//Perform a DFS to sort the leaves and find every link's child
	if(my_rank == 0)	printf("Performing DFS...\n");
	MPI_Barrier(MPI_COMM_WORLD);
	start = time(NULL);

	Link** stack = malloc(*N * sizeof(Link*)); //Holds the index in system
	int stack_size = 0;

	Link** leaves = malloc(*N * sizeof(Link*));
	unsigned int leaves_size = 0;
	unsigned short int numparents;

	for(j=0;j<*N;j++)	//!!!! Iterate over a list of roots? !!!!
	{
		if(upstream_order[j]->c == NULL)	//If upstream_order[j] doesn't have a child set yet, it must be an outlet
		{
			stack[0] = upstream_order[j];
			stack_size = 1;

			while(stack_size > 0)
			{
				current = stack[stack_size-1];	//Top of stack
				numparents = current->numparents;

				if(numparents == 0)
				{
					stack_size--;
					leaves[leaves_size] = current;
					leaves_size++;
				}
				else
				{
					//If current is not a leaf, replace it with its parents
					for(i=0;i<numparents;i++)
					{
						stack[stack_size - 1 + i] = current->parents[numparents - 1 - i];
						//stack[stack_size - 1 + i]->c = current;
					}
					stack_size += numparents - 1;
				}
			}
		}
	}

/*
	for(j=0;j<*N;j++)	//!!!! Iterate over a list of roots? !!!!
	{
		if(upstream_order[j]->c == NULL)	//If upstream_order[j] doesn't have a child set yet, it must be an outlet
		{
			stack[0] = upstream_order[j];
			stack_size = 1;

			while(stack_size > 0)
			{
				current = stack[stack_size-1];	//Top of stack
				numparents = current->numparents;

				if(numparents == 0)
				{
					stack_size--;
					leaves[leaves_size] = current;
					leaves_size++;
				}
				else
				{
					if(current->parents[numparents-1]->c != NULL)	//If not null, this link has been visited already
					{
						stack_size--;
					}
					else
					{
						//If current is not a leaf, replace it with its parents
						for(i=0;i<numparents;i++)
						{
							stack[stack_size - 1 + i] = current->parents[numparents - 1 - i];
							stack[stack_size - 1 + i]->c = current;
						}
						stack_size += numparents - 1;
					}
				}
			}
		}
	}
*/
	free(upstream_order);
/*
for(i=0;i<*N;i++)
{
	if(system[i]->numparents == 0)
	{
		for(j=0;j<leaves_size;j++)
		{
			if(leaves[j]->ID == system[i]->ID)	break;
		}
		if(j == leaves_size)
			printf("Leaf %u not identified.\n",system[i]->ID);
	}
}

printf("***********\n");
printf("number of leaves = %u\n",leaves_size);
printf("***********\n");

for(i=0;i<*N;i++)
{
if(system[i]->c == NULL)
	printf("%u is an outlet.\n",system[i]->ID);
}
getchar();
*/
	MPI_Barrier(MPI_COMM_WORLD);
	stop = time(NULL);
	if(my_rank == 0)	printf("Time for DFS: %f\n",difftime(stop,start));

	//Cleanup a bit
	for(i=0;i<*N;i++)	free(loc_to_children[i]);
	free(loc_to_children);

	//Calculate the distance and number of upstream links for each link
	if(my_rank == 0)	printf("Calculating distances...\n");
	MPI_Barrier(MPI_COMM_WORLD);
	start = time(NULL);

	//Distance
	for(i=0;i<*N;i++)	system[i]->distance = 0;

	for(i=0;i<leaves_size;i++)
	{
		prev = leaves[i];
		prev->distance = 1;
		for(current = prev->c; current != NULL; current = current->c)
		{
			if(current->distance > prev->distance+1)	break;
			else						current->distance = prev->distance + 1;
			prev = current;
		}
	}

	for(i=0;i<*N;i++)	system[i]->upstream = NULL;

	//For extended data assimilation
	//!!!! Upstream link calculations, and AssimError above. !!!!
	if(GlobalVars->assim_flag == 1)
	{
		unsigned int* temp_numupstream = (unsigned int*) calloc(*N,sizeof(unsigned int));
		for(i=0;i<leaves_size;i++)
			temp_numupstream[leaves[i]->location] = 1;

		//Count upstream links
		for(i=0;i<leaves_size;i++)
		{
			for(current = leaves[i]->c; current != NULL; current = current->c)
			{
				parentsval = 0;
				for(j=0;j<current->numparents;j++)	parentsval += (temp_numupstream[current->parents[j]->location] > 0);

				if(parentsval == current->numparents)	//All parents have temp_numupstream set
				{
					temp_numupstream[current->location] = 1;
					for(j=0;j<current->numparents;j++)
						temp_numupstream[current->location] += temp_numupstream[current->parents[j]->location];
				}
				else
					break;
			}
		}

		//Set the upstream links
		unsigned int** temp_upstream = (unsigned int**) malloc(*N*sizeof(unsigned int*));	//temp_upstream[i] is list of all links upstream from link i
		for(i=0;i<*N;i++)
			temp_upstream[i] = (unsigned int*) malloc(temp_numupstream[system[i]->location] * sizeof(unsigned int));
		unsigned int* counter = (unsigned int*) calloc(*N,sizeof(unsigned int));

		stack_size = leaves_size;
		for(i=0;i<leaves_size;i++)	stack[i] = leaves[i];

		while(stack_size > 0)
		{
			current = stack[stack_size-1];
			l = current->location;

			//Add this link to its own upstream list
			temp_upstream[l][counter[l]] = l;
			counter[l]++;

			//Add each parents upstream list
			for(i=0;i<current->numparents;i++)
			{
				m = current->parents[i]->location;
				for(j=0;j<counter[m];j++)
					temp_upstream[l][counter[l]+j] = temp_upstream[m][j];
				counter[l] += counter[m];
			}

			stack_size--;

			//If every parent of current's child has an upstream list determined, add it to the stack
			if(current->c != NULL)
			{
				parentsval = 0;
				for(i=0;i<current->c->numparents;i++)
				{
					m = current->c->parents[i]->location;
					parentsval += (counter[m] > 0);
				}

				if(parentsval == current->c->numparents)
				{
					stack[stack_size] = current->c;
					stack_size++;
				}
			}
		}

		//Move the data from temp_upstream into the child upstream
		for(i=0;i<*N;i++)
		{
			system[i]->upstream = (unsigned int**) malloc(system[i]->numparents * sizeof(unsigned int*));
			for(j=0;j<system[i]->numparents;j++)
				system[i]->upstream[j] = temp_upstream[system[i]->parents[j]->location];
		}

		//Set number of upstream link to parents
		for(i=0;i<*N;i++)
		{
			system[i]->numupstream = (unsigned int*) malloc(system[i]->numparents * sizeof(unsigned int));
			for(j=0;j<system[i]->numparents;j++)
				system[i]->numupstream[j] = temp_numupstream[system[i]->parents[j]->location];
		}

		//Cleanup
		for(i=0;i<*N;i++)
			if(system[i]->c == NULL)	free(temp_upstream[i]);
		free(temp_upstream);
		free(temp_numupstream);
		free(counter);
/*
		//Try removing low order links from the upstream lists
		printf("!!!! Removing low order links from upstream list...!!!!\n");
		unsigned int cut_off = 2,drop;
		unsigned int* order = (unsigned int*) malloc(*N*sizeof(unsigned int));
		unsigned short int* complete = (unsigned short int*) malloc(*N*sizeof(unsigned short int));
		CalcHortonOrder(system,*N,order,complete);

		for(i=0;i<*N;i++)
		{
			for(j=0;j<system[i]->numparents;j++)
			{
				drop = 0;
//printf("checking %i (%i)\n",system[i]->ID,system[i]->parents[j]->ID);
				for(k=0;k<system[i]->numupstream[j];k++)
				{
					if(order[system[i]->upstream[j][k]] <= cut_off)	//!!!! Once this happens, the rest of the links can be dropped !!!!
					{
						drop++;
					}
					else
					{
						system[i]->upstream[j][k-drop] = system[i]->upstream[j][k];
					}
				}
//printf("%i going to %i from %i\n",system[i]->ID,system[i]->numupstream[j]-drop,system[i]->numupstream[j]);
				system[i]->numupstream[j] -= drop;
			}
		}

		free(order);
		free(complete);
*/
	}

	MPI_Barrier(MPI_COMM_WORLD);
	stop = time(NULL);
	if(my_rank == 0)	printf("Time to calculate distances: %f\n",difftime(stop,start));


	//Partition the system and assign the links
	if(my_rank == 0)	printf("Partitioning system...\n");
	MPI_Barrier(MPI_COMM_WORLD);
	start = time(NULL);

	*my_data = Initialize_TransData();
	*getting = (short int*) malloc(*N*sizeof(short int));
	*assignments = Partition_System_By_Leaves(system,*N,leaves,leaves_size,my_sys,my_N,my_max_nodes,*my_data,*getting);
	//*assignments = Partition_System_By_Leaves_2(system,*N,leaves,leaves_size,my_sys,my_N,my_max_nodes,*my_data,*getting);
	//*assignments = Partition_METIS_Traditional(system,*N,leaves,leaves_size,my_sys,my_N,my_max_nodes,*my_data,*getting,GlobalVars);
	//*assignments = Partition_METIS_RainChanges(system,*N,leaves,leaves_size,my_sys,my_N,my_max_nodes,*my_data,*getting,GlobalVars);
	//*assignments = Partition_METIS_RainVolume(system,*N,leaves,leaves_size,my_sys,my_N,my_max_nodes,*my_data,*getting,GlobalVars);	//!!!! Requires params for all links !!!!

	free(leaves);

	//Allocate space for communication related data
	//Need space for nodes, number of iterations, discontinuities
	//Data: ( size(double)*(max_s*max_dim + max_dim + time)*# steps to transfer + size(int)(for stage)*# steps to transfer + size(int)*(location + # steps to transfer) ) * # of sending links
	//Upstream: + size(int) * (location + # of iterations) * # of receiving links
	//Discontinuities: + (size(int) + size(double)*discont_size + size(int)*discont_size) * # of sending links
	unsigned int bytes1 = ( (sizeof(double)*(GlobalVars->max_s*GlobalVars->num_dense + GlobalVars->dim + 1) + sizeof(int) )*GlobalVars->max_transfer_steps + sizeof(int)*2);
	unsigned int bytes2 = 2*sizeof(int);
	unsigned int bytes3 = sizeof(int) + (sizeof(int) + sizeof(double))*GlobalVars->discont_size;
	for(ii=0;ii<np;ii++)
	{
		(*my_data)->send_buffer_size[ii] = bytes1 * (*my_data)->send_size[ii] + bytes2 * (*my_data)->receive_size[ii] + bytes3 * (*my_data)->send_size[ii];
		(*my_data)->receive_buffer_size[ii] = bytes1 * (*my_data)->receive_size[ii] + bytes2*(*my_data)->send_size[ii] + bytes3 * (*my_data)->receive_size[ii];

		if((*my_data)->send_buffer_size[ii])	(*my_data)->send_buffer[ii] = malloc((*my_data)->send_buffer_size[ii]*sizeof(char));
		else					(*my_data)->send_buffer[ii] = NULL;

		if((*my_data)->receive_buffer_size[ii])	(*my_data)->receive_buffer[ii] = malloc((*my_data)->receive_buffer_size[ii]*sizeof(char));
		else					(*my_data)->receive_buffer[ii] = NULL;
	}

	MPI_Barrier(MPI_COMM_WORLD);
	stop = time(NULL);
	if(my_rank == 0) 	printf("Time for partitioning: %f\n",difftime(stop,start));

	//Read uniform initial data
	if(GlobalVars->init_flag == 1)
	{
		uniform_y_0 = v_get(dim);
		for(i=GlobalVars->diff_start;i<GlobalVars->no_ini_start;i++)
		//for(i=GlobalVars->diff_start;i<initdim;i++)
		{
			if( fscanf(initdata,"%lf",&(uniform_y_0->ve[i])) == 0 )
			{
				printf("[%i]: Error reading .uini file: Not enough initial states.\n",my_rank);
				MPI_Abort(MPI_COMM_WORLD,1);
			}
		}
	}
	else if(GlobalVars->init_flag == 3)	//Initial data from database
	{
		//Allocate space 
		db_init_buffer = (VEC**) calloc(*N,sizeof(VEC*));
		double buffer[dim];
		unsigned int loc;

		if(my_rank == 0)
		{
			if(ConnectPGDB(db_connections[ASYNCH_DB_LOC_INIT]))
			{
				printf("[%i]: Error connecting to database for init conditions.\n",my_rank);
				MPI_Abort(MPI_COMM_WORLD,1);
			}
			sprintf(db_connections[ASYNCH_DB_LOC_INIT]->query,db_connections[ASYNCH_DB_LOC_INIT]->queries[0],GlobalVars->init_timestamp);
			res = PQexec(db_connections[ASYNCH_DB_LOC_INIT]->conn,db_connections[ASYNCH_DB_LOC_INIT]->query);
			CheckResError(res,"downloading init data");

			if(PQntuples(res) != *N)
			{
				printf("[%i]: Error downloading init data. Got %i conditions, expected %u.\n",my_rank,PQntuples(res),*N);
				MPI_Abort(MPI_COMM_WORLD,1);
			}
			if(PQnfields(res) != dim + 1)
			{
				printf("[%i]: Error downloading init data. Expected link id and %u states from initial condition data, but the query only gives %u columns.\n",my_rank,dim,PQnfields(res));
				MPI_Abort(MPI_COMM_WORLD,1);
			}

			for(i=0;i<*N;i++)
			{
				loc = find_link_by_idtoloc(atoi(PQgetvalue(res,i,0)),*id_to_loc,*N);
				for(j=0;j<dim;j++)	buffer[j] = atof(PQgetvalue(res,i,j+1));
				MPI_Bcast(&loc,1,MPI_INT,0,MPI_COMM_WORLD);
				MPI_Bcast(buffer,dim,MPI_DOUBLE,0,MPI_COMM_WORLD);
				if((*assignments)[loc] == my_rank || (*getting)[loc])
				{
					db_init_buffer[loc] = v_get(dim);
					for(j=0;j<dim;j++)
						db_init_buffer[loc]->ve[j] = buffer[j];
				}
			}

			PQclear(res);
			DisconnectPGDB(db_connections[ASYNCH_DB_LOC_INIT]);
		}
		else
		{
			for(i=0;i<*N;i++)
			{
				MPI_Bcast(&loc,1,MPI_INT,0,MPI_COMM_WORLD);
				MPI_Bcast(buffer,dim,MPI_DOUBLE,0,MPI_COMM_WORLD);

				if((*assignments)[loc] == my_rank || (*getting)[loc])
				{
					db_init_buffer[loc] = v_get(dim);
					for(j=0;j<dim;j++)
						db_init_buffer[loc]->ve[j] = buffer[j];
				}
			}
		}
	}

	//Setup forcings. Read uniform forcing data and open .str files. Also initialize rainfall from database.
	FILE** forcingfiles = (FILE**) malloc(GlobalVars->num_forcings*sizeof(FILE*));
	for(i=0;i<GlobalVars->num_forcings;i++)
	{
		forcings[i]->maxtime = GlobalVars->t_0;
		forcings[i]->iteration = 0;
		//forcings[i]->passes = 0;
		forcings[i]->active = 1;

		switch(forcings[i]->flag)
		{
			case 0:
				forcings[i]->GetPasses = &PassesOther;
				forcings[i]->GetNextForcing = &NextForcingOther;
				break;
			case 1:
				forcings[i]->GetPasses = &PassesOther;
				forcings[i]->GetNextForcing = &NextForcingOther;
				break;
			case 2:
				forcings[i]->GetPasses = &PassesBinaryFiles;
				forcings[i]->GetNextForcing = &NextForcingBinaryFiles;
				break;
			case 3:
				forcings[i]->GetPasses = &PassesDatabase;
				forcings[i]->GetNextForcing = &NextForcingDatabase;
				break;
			case 4:
				forcings[i]->GetPasses = &PassesOther;
				forcings[i]->GetNextForcing = &NextForcingOther;
				break;
			case 5:
				//!!!! I think these are ok. But option 5 should die. !!!!
				forcings[i]->GetPasses = &PassesDatabase;
				forcings[i]->GetNextForcing = &NextForcingDatabase;
				break;
			case 6:
				forcings[i]->GetPasses = &PassesBinaryFiles;
				forcings[i]->GetNextForcing = &NextForcingGZBinaryFiles;
				break;
			case 7:
				forcings[i]->GetPasses = &PassesRecurring;
				forcings[i]->GetNextForcing = &NextForcingRecurring;
				break;
			default:
				printf("[%i]: Error: Bad forcing flag (%hu) in forcing %u.\n",my_rank,forcings[i]->flag,i);
				MPI_Abort(MPI_COMM_WORLD,1);
		}

		//Read in files for uniform in space data
		if(forcings[i]->flag == 4)
		{
			forcingfiles[i] = fopen(forcings[i]->filename,"r");
			if(!forcingfiles[i])
			{
				printf("[%i]: Error: cannot open uniform forcing file %s.\n",my_rank,forcings[i]->filename);
				MPI_Abort(MPI_COMM_WORLD,1);
			}

			//Read the number of times for the rain data for this link
			fscanf(forcingfiles[i],"%i",&m);
			m++;		//Increase this by one to add a "ceiling" term

			forcings[i]->GlobalForcing = (ForcingData*) malloc(sizeof(ForcingData));
			forcings[i]->GlobalForcing->rainfall = (double**) malloc(m*sizeof(double*));
			forcings[i]->GlobalForcing->n_times = m;
			for(j=0;j<m;j++)	forcings[i]->GlobalForcing->rainfall[j] = (double*) malloc(2*sizeof(double));
			m--;

			//Read in the storm data for this link
			for(j=0;j<m;j++)
			{
				double tempy;
				fscanf(forcingfiles[i],"%lf",&tempy);
				(forcings[i]->GlobalForcing->rainfall[j][0]) = tempy;
				fscanf(forcingfiles[i],"%lf",&(forcings[i]->GlobalForcing->rainfall[j][1]));
			}
			forcings[i]->GlobalForcing->rainfall[m][0] = GlobalVars->maxtime + 3.0;
			forcings[i]->GlobalForcing->rainfall[m][1] = -1;

			double rainfall_buffer = forcings[i]->GlobalForcing->rainfall[0][1];
			for(j=1;j<forcings[i]->GlobalForcing->n_times;j++)
			{
				if(rainfall_buffer != forcings[i]->GlobalForcing->rainfall[j][1])
				{
					univ_forcing_change_time[i] = forcings[i]->GlobalForcing->rainfall[j][0];
					break;
				}
			}
			if(j == forcings[i]->GlobalForcing->n_times)
				univ_forcing_change_time[i] = forcings[i]->GlobalForcing->rainfall[j-1][0];

			fclose(forcingfiles[i]);
		}
		else if(forcings[i]->flag == 7)	//Read in a monthly file
		{
			forcingfiles[i] = fopen(forcings[i]->filename,"r");
			if(!forcingfiles[i])
			{
				printf("[%i]: Error: cannot open uniform forcing file %s.\n",my_rank,forcings[i]->filename);
				MPI_Abort(MPI_COMM_WORLD,1);
			}

			//I'm assuming there are 12 months in a year. I hope that assumption isn't too strong...
			int num_read,num_months = 12;
			forcings[i]->GlobalForcing = (ForcingData*) malloc(sizeof(ForcingData));
			forcings[i]->GlobalForcing->rainfall = (double**) malloc((num_months+1)*sizeof(double*));
			forcings[i]->GlobalForcing->n_times = num_months + 1;
			for(j=0;j<num_months+1;j++)	forcings[i]->GlobalForcing->rainfall[j] = (double*) malloc(2*sizeof(double));

			for(j=0;j<num_months && !feof(forcingfiles[i]);j++)
				num_read = fscanf(forcingfiles[i],"%lf",&(forcings[i]->GlobalForcing->rainfall[j][1]));

			if(feof(forcingfiles[i]) || !num_read)
			{
				printf("[%i]: Error reading from monthly forcing file %s.\n",my_rank,forcings[i]->filename);
				MPI_Abort(MPI_COMM_WORLD,1);
			}

			//Find the starting month, and use the next month for the change time
			time_t start_time_t = forcings[i]->first_file;
			struct tm *start_time = gmtime(&start_time_t);

			(start_time->tm_mon)++;
			start_time->tm_mday = 1;
			start_time->tm_hour = 0;
			start_time->tm_min = 0;
			start_time->tm_sec = 0;
			time_t next_month = mktime(start_time);
			univ_forcing_change_time[i] = difftime(next_month,forcings[i]->first_file)/60.0;

			fclose(forcingfiles[i]);
		}
		else if(forcings[i]->flag == 1)	//Open ascii files
		{
			forcingfiles[i] = fopen(forcings[i]->filename,"r");
			if(!forcingfiles[i])
			{
				printf("[%i]: Error: cannot open uniform forcing file %s.\n",my_rank,forcings[i]->filename);
				MPI_Abort(MPI_COMM_WORLD,1);
			}

			fscanf(forcingfiles[i],"%u",&m);
			if(m != *N && (!(GlobalVars->res_flag) || i != GlobalVars->res_forcing_idx))
			{
				printf("[%i]: Error: Number of links in .str file differs from number of links in network (%u vs %u).\n",my_rank,m,*N);
				MPI_Abort(MPI_COMM_WORLD,1);
			}
		}
		else if(forcings[i]->flag == 3 || forcings[i]->flag == 5)
		{
			unsigned int good_time,is_null;
			if(my_rank == 0)
			{
				ConnectPGDB(db_connections[ASYNCH_DB_LOC_FORCING_START+i]);
				res = PQexec(db_connections[ASYNCH_DB_LOC_FORCING_START+i]->conn,db_connections[ASYNCH_DB_LOC_FORCING_START+i]->queries[2]);
				CheckResError(res,"querying a valid forcing time");
				is_null = PQgetisnull(res,0,0);
				if(is_null)	printf("[%i]: Warning: forcing %u has no data...\n",my_rank,i);
				else		good_time = atoi(PQgetvalue(res,0,0));
				PQclear(res);
				DisconnectPGDB(db_connections[ASYNCH_DB_LOC_FORCING_START+i]);
			}

			MPI_Bcast(&is_null,1,MPI_INT,0,MPI_COMM_WORLD);
			if(is_null)
				forcings[i]->good_timestamp = forcings[i]->raindb_start_time;
			else
			{
				MPI_Bcast(&good_time,1,MPI_INT,0,MPI_COMM_WORLD);
				forcings[i]->good_timestamp = good_time;
			}	
		}
	}

	//Read dam file
	if(GlobalVars->uses_dam && GlobalVars->dam_flag == 1)	//.dam file
	{
		unsigned int num_dams;
		VEC* newparams;
		FILE* damfile = fopen(GlobalVars->dam_filename,"r");
		if(damfile == NULL)
		{
			printf("[%i]: Dam file %s not found.\n",my_rank,GlobalVars->dam_filename);
			return NULL;
		}
		
		fscanf(damfile,"%u",&num_dams);

		for(i=0;i<num_dams;i++)
		{
			fscanf(damfile,"%u",&id);
			m = find_link_by_idtoloc(id,*id_to_loc,*N);
			if(m > *N)
			{
				printf("[%i]: Error: Link id %u from .dam file not present in network.\n",my_rank,id);
				MPI_Abort(MPI_COMM_WORLD,1);
			}

			if(my_rank == (*assignments)[m] || (*getting)[m] == 1)
			{
				newparams = v_get(GlobalVars->dam_params_size);

				for(j=GlobalVars->params_size;j<GlobalVars->dam_params_size;j++)
					fscanf(damfile,"%lf",&(newparams->ve[j]));

				for(j=0;j<GlobalVars->params_size;j++)		//!!!! Use realloc !!!!
					newparams->ve[j] = system[m]->params->ve[j];
				v_free(system[m]->params);
				system[m]->params = newparams;
				system[m]->dam = 1;
			}
			else
			{
				for(j=GlobalVars->params_size;j<GlobalVars->dam_params_size;j++)
					fscanf(damfile,"%*f");
			}
		}

		//Set error tolerance
		if(rk_filename[0] == '\0')
		{
			GlobalVars->discont_tol = GlobalErrors->abstol->ve[0];
			if(GlobalVars->discont_tol < 1e-12 && my_rank == 0)
				printf("Warning: Discontinuity tolerance has been set to %e.\n",GlobalVars->discont_tol);
		}
		else
		{
			GlobalVars->discont_tol = 1e-8;
			if(my_rank == 0)	printf("Notice: Discontinuity tolerance has been set to %e.\n",GlobalVars->discont_tol);
		}

		fclose(damfile);
	}
	else if(GlobalVars->uses_dam && GlobalVars->dam_flag == 2)	//.qvs file
	{
		unsigned int num_dams,num_values;
		FILE* damfile = fopen(GlobalVars->dam_filename,"r");
		if(damfile == NULL)
		{
			printf("Dam file %s not found.\n",GlobalVars->dam_filename);
			return NULL;
		}
		
		fscanf(damfile,"%u",&num_dams);
		
		for(i=0;i<num_dams;i++)
		{
			fscanf(damfile,"%u %u",&id,&num_values);

			m = find_link_by_idtoloc(id,*id_to_loc,*N);
			if(m > *N)
			{
				printf("[%i]: Error: Link id %u from .qvs file not present in network.\n",my_rank,id);
				MPI_Abort(MPI_COMM_WORLD,1);
			}

			if(my_rank == (*assignments)[m] || (*getting)[m] == 1)
			{
				current = system[m];
				current->dam = 1;
				current->qvs = (QVSData*) malloc(sizeof(QVSData));
				current->qvs->n_values = num_values;
				current->qvs->points_array = (double*) malloc(2*num_values*sizeof(double));
				current->qvs->points = (double**) malloc(num_values*sizeof(double*));

				for(j=0;j<num_values;j++)	current->qvs->points[j] = &(current->qvs->points_array[2*j]);
				//for(j=0;j<num_values;j++)	current->qvs->points[j] = (double*) malloc(2*sizeof(double));

				for(j=0;j<num_values;j++)
					fscanf(damfile,"%lf %lf",&(current->qvs->points[j][0]),&(current->qvs->points[j][1]));
			}
			else	//Just read through to the next dam
			{
				for(j=0;j<num_values;j++)	fscanf(damfile,"%*f %*f");
			}
		}

		//Set error tolerance
		if(rk_filename[0] == '\0')
		{
			GlobalVars->discont_tol = GlobalErrors->abstol->ve[0];
			if(GlobalVars->discont_tol < 1e-12 && my_rank == 0)
				printf("Warning: Discontinuity tolerance has been set to %e.\n",GlobalVars->discont_tol);
		}
		else
		{
			GlobalVars->discont_tol = 1e-8;
			if(my_rank == 0)	printf("Notice: Discontinuity tolerance has been set to %e.\n",GlobalVars->discont_tol);
		}

		fclose(damfile);
	}
	else if(GlobalVars->uses_dam && GlobalVars->dam_flag == 3)	//database connection
	{
		unsigned int num_dams = 0,num_pts,num_values;
		short int procs_sending_to[np],mine;
		double* array_holder;
		MPI_Status status;

		//Connect to the database and download data
		if(my_rank == 0)
		{
			ConnectPGDB(db_connections[ASYNCH_DB_LOC_QVS]);
			res = PQexec(db_connections[ASYNCH_DB_LOC_QVS]->conn,db_connections[ASYNCH_DB_LOC_QVS]->queries[0]);
			if(CheckResError(res,"querying qvs relations"))	return NULL;
			num_pts = PQntuples(res);

			i = 0;
			while(i < num_pts)
			{
				//Get link id and location
				id = atoi(PQgetvalue(res,i,0));
				curr_loc = find_link_by_idtoloc(id,*id_to_loc,*N);

				//Count number of qvs values for current
				for(j=i;j<num_pts && id == atoi(PQgetvalue(res,j,0));j++);
				num_values = j-i;

				//Extract the data points
				array_holder = (double*) malloc(2*num_values*sizeof(double));
				for(j=0;j<num_values;j++)
				{
					array_holder[2*j] = atof(PQgetvalue(res,i+j,1));
					array_holder[2*j+1] = atof(PQgetvalue(res,i+j,2));
				}

				//Check for an error real quick
				for(j=1;j<num_values;j++)
				{
					if(array_holder[2*(j-1)] > array_holder[2*j] || array_holder[2*(j-1)+1] > array_holder[2*j+1])
					{
						printf("[%i]: Bad storage or discharge values found at link id %u. Check that the data is sorted correctly. (%u)\n",my_rank,id,j);
						//break;
						free(array_holder);
						return NULL;
					}
				}

				if(curr_loc < *N)
				{
					//Tell everyone what link has the current dam
					MPI_Bcast(&curr_loc,1,MPI_INT,0,MPI_COMM_WORLD);
					mine = (my_rank == (*assignments)[curr_loc] || (*getting)[curr_loc]);
					MPI_Gather(&mine,1,MPI_SHORT,procs_sending_to,1,MPI_SHORT,0,MPI_COMM_WORLD);

					//Send the data to whoever needs it
					for(j=1;j<np;j++)
					{
						if(procs_sending_to[j])
						{
							MPI_Send(&num_values,1,MPI_UNSIGNED,j,0,MPI_COMM_WORLD);
							MPI_Send(array_holder,2*num_values,MPI_DOUBLE,j,0,MPI_COMM_WORLD);
						}
					}

					//Check if proc 0 needs the data
					if(mine)
					{
						current = system[curr_loc];
						current->dam = 1;
						current->qvs = (QVSData*) malloc(sizeof(QVSData));
						current->qvs->n_values = num_values;
						current->qvs->points_array = array_holder;
						current->qvs->points = (double**) malloc(num_values*sizeof(double*));
						for(j=0;j<num_values;j++)	current->qvs->points[j] = &(current->qvs->points_array[2*j]);
					}
					else
						free(array_holder);
				}

				i += num_values;
			}

			//Finish up
			PQclear(res);
			curr_loc = -1;
			MPI_Bcast(&curr_loc,1,MPI_INT,0,MPI_COMM_WORLD);
		}
		else	//Receive dam data
		{
			curr_loc = 0;
			MPI_Bcast(&curr_loc,1,MPI_INT,0,MPI_COMM_WORLD);

			while((int)curr_loc != -1)
			{
				//Check if I need the data for this link
				mine = (my_rank == (*assignments)[curr_loc] || (*getting)[curr_loc]);
				MPI_Gather(&mine,1,MPI_SHORT,procs_sending_to,1,MPI_SHORT,0,MPI_COMM_WORLD);

				//Receive data
				if(mine)
				{
					MPI_Recv(&num_values,1,MPI_UNSIGNED,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
					array_holder = (double*) malloc(2*num_values*sizeof(double));
					MPI_Recv(array_holder,2*num_values,MPI_DOUBLE,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
					current = system[curr_loc];
					current->dam = 1;
					current->qvs = (QVSData*) malloc(sizeof(QVSData));
					current->qvs->n_values = num_values;
					current->qvs->points_array = array_holder;
					current->qvs->points = (double**) malloc(num_values*sizeof(double*));
					for(j=0;j<num_values;j++)	current->qvs->points[j] = &(current->qvs->points_array[2*j]);
				}

				//Check next signal
				MPI_Bcast(&curr_loc,1,MPI_INT,0,MPI_COMM_WORLD);
			}
		}

		//Set error tolerance
		if(rk_filename[0] == '\0')
		{
			GlobalVars->discont_tol = GlobalErrors->abstol->ve[0];
			if(GlobalVars->discont_tol < 1e-12 && my_rank == 0)
				printf("Warning: Discontinuity tolerance has been set to %e.\n",GlobalVars->discont_tol);
		}
		else
		{
			GlobalVars->discont_tol = 1e-8;
			if(my_rank == 0)	printf("Notice: Discontinuity tolerance has been set to %e.\n",GlobalVars->discont_tol);
		}

	}
	else	//Some other type of discontinuity (or none at all)
	{
		//Set error tolerance
		if(rk_filename[0] == '\0')
		{
			GlobalVars->discont_tol = GlobalErrors->abstol->ve[0];
			if(GlobalVars->discont_tol < 1e-12 && my_rank == 0)
				printf("Warning: Discontinuity tolerance has been set to %e.\n",GlobalVars->discont_tol);
		}
		else
		{
			GlobalVars->discont_tol = 1e-8;
			//if(my_rank == 0)	printf("Notice: Discontinuity tolerance has been set to %e.\n",GlobalVars->discont_tol);
		}
	}

	//Find links with state forcing
	unsigned int res_size,loc,*res_list = NULL;
	if(GlobalVars->res_flag)
	{
		res_list = Create_SAV_Data(GlobalVars->rsv_filename,system,*N,&res_size,db_connections[ASYNCH_DB_LOC_RSV],GlobalVars->res_flag);

		//Setup links with forcing
		for(j=0;j<res_size;j++)
		{
			loc = find_link_by_idtoloc(res_list[j],*id_to_loc,*N);
			if(loc < *N)
			{
				system[loc]->res = 1;
			}
			else
			{
				if(my_rank == 0)
					printf("[%i]: Warning: ignoring reservoir at ID %u. ID not found in network.\n",my_rank,res_list[j]);
			}
		}
	}

	//Second pass: Read in the data for the necessary links.
	if(my_rank == 0)	printf("Second read...\n");
	MPI_Barrier(MPI_COMM_WORLD);
	start = time(NULL);

	for(i=0;i<*N;i++)
	{
		if(rk_filename[0] != '\0')	fscanf(rkdata,"%i",&j);
		else 				j = system[i]->ID;
		if(GlobalVars->init_flag != 1 && initdata)	fscanf(initdata,"%i",&k);
		else				k = system[i]->ID;
		//if(GlobalVars->rain_flag == 1)	fscanf(stormdata,"%i",&m);

		if(!( (system[i]->ID == j) && (j == k) ))
		{
			if(my_rank == 0)
			{
				printf("[%i]: Error: Links must be listed in the same order in every data file. %u %u %u\n",my_rank,system[i]->ID,j,k);

				//Give further hints about the error
				if(system[i]->ID == j && j != k)
					printf("[%i]: Hint: the problem seems to be in the initial condition source.\n",my_rank);
				else if(system[i]->ID != j && j == k)
					printf("[%i]: Hint: the problem seems to be in the topology data source.\n",my_rank);
				else if(system[i]->ID == k && j != k)
					printf("[%i]: Hint: the problem seems to be in the Runge-Kutta data source.\n",my_rank);
				else
					printf("[%i]: Hint: check the topology data source, initial condition source, and Runge-Kutta data source.\n",my_rank);
			}

			*N=0;
			return NULL;
		}

		//If assigned to this process or if data is needed from it
		if(my_rank == (*assignments)[i] || (*getting)[i])
		{
			if(rk_filename[0] != '\0')
			{
				system[i]->errorinfo = malloc(sizeof(ErrorData));
				system[i]->errorinfo->abstol = v_get(dim);
				system[i]->errorinfo->reltol = v_get(dim);
				system[i]->errorinfo->abstol_dense = v_get(dim);
				system[i]->errorinfo->reltol_dense = v_get(dim);
				system[i]->errorinfo->facmax = GlobalErrors->facmax;
				system[i]->errorinfo->facmin = GlobalErrors->facmin;
				system[i]->errorinfo->fac = GlobalErrors->fac;

				for(j=0;j<dim;j++)	fscanf(rkdata,"%lf",&(system[i]->errorinfo->abstol->ve[j]));
				//for(j=dim;j<rkdim;j++)	fscanf(rkdata,"%lf",&trash);	//Trash any unneeded data
				for(j=0;j<dim;j++)	fscanf(rkdata,"%lf",&(system[i]->errorinfo->reltol->ve[j]));
				//for(j=dim;j<rkdim;j++)	fscanf(rkdata,"%lf",&trash);
				for(j=0;j<dim;j++)	fscanf(rkdata,"%lf",&(system[i]->errorinfo->abstol_dense->ve[j]));
				//for(j=dim;j<rkdim;j++)	fscanf(rkdata,"%lf",&trash);
				for(j=0;j<dim;j++)	fscanf(rkdata,"%lf",&(system[i]->errorinfo->reltol_dense->ve[j]));
				//for(j=dim;j<rkdim;j++)	fscanf(rkdata,"%lf",&trash);

				fscanf(rkdata,"%i",&data1);
				system[i]->method = (*AllMethods)[data1];
			}
			else
			{
				system[i]->errorinfo = GlobalErrors;
				system[i]->method = (*AllMethods)[GlobalVars->method];
			}

			//GlobalVars->max_s = (GlobalVars->max_s > system[i]->method->s + 1) ? GlobalVars->max_s : system[i]->method->s + 1;
			//!!!! Note: added +1 for Radau solver !!!!
			//!!!! This should be done outside this loop !!!!

			//Set solver and ODE
			if(custom_model)	custom_model->Routines(system[i],type,system[i]->method->exp_imp,system[i]->dam,external);
			else			InitRoutines(system[i],type,system[i]->method->exp_imp,system[i]->dam,external);

			//Setup the initial stepsize
			//!!!! This should be removed, but requires changing .rkd file structure !!!!
			//if(rk_filename[0] != '\0')	fscanf(rkdata,"%lf",&(system[i]->h));

			//Make some precalculations
			if(custom_model)	custom_model->Precalculations(system[i],GlobalVars->global_params,system[i]->params,system[i]->iparams,GlobalVars->disk_params,GlobalVars->params_size,system[i]->dam,type,external);
			else			Precalculations(system[i],GlobalVars->global_params,system[i]->params,system[i]->iparams,GlobalVars->disk_params,GlobalVars->params_size,system[i]->dam,type,external);

			//Setup the initial values
			if(GlobalVars->init_flag == 0)
			{
				for(j=GlobalVars->diff_start;j<GlobalVars->no_ini_start;j++)	fscanf(initdata,"%lf",&(y_0->ve[j]));
				//for(j=GlobalVars->no_ini_start;j<dim;j++)		fscanf(initdata,"%*f");	//Trash unneeded data	//!!!! This was uncommented. I don't think it should exist... !!!!
				if(custom_model)	system[i]->state = custom_model->InitializeEqs(GlobalVars->global_params,system[i]->params,system[i]->iparams,system[i]->qvs,system[i]->dam,y_0,type,GlobalVars->diff_start,GlobalVars->no_ini_start,external);
				else			system[i]->state = ReadInitData(GlobalVars->global_params,system[i]->params,system[i]->iparams,system[i]->qvs,system[i]->dam,y_0,type,GlobalVars->diff_start,GlobalVars->no_ini_start,external);
			}
			else if(GlobalVars->init_flag == 1)
			{
				v_copy(uniform_y_0,y_0);
				if(custom_model)	system[i]->state = custom_model->InitializeEqs(GlobalVars->global_params,system[i]->params,system[i]->iparams,system[i]->qvs,system[i]->dam,y_0,type,GlobalVars->diff_start,GlobalVars->no_ini_start,external);
				else			system[i]->state = ReadInitData(GlobalVars->global_params,system[i]->params,system[i]->iparams,system[i]->qvs,system[i]->dam,y_0,type,GlobalVars->diff_start,GlobalVars->no_ini_start,external);
			}
			else if(GlobalVars->init_flag == 2)	//Reading a .rec file
			{
				for(j=0;j<dim;j++)	fscanf(initdata,"%lf",&(y_0->ve[j]));
/*
				if(system[i]->alg != NULL)
				{
					//if(GlobalVars->type < 21)	printf("Warning: Type should be 21 for dam_check (in riversys.c).\n");
					//system[i]->state = dam_check(y_0,GlobalVars->global_params,system[i]->params,system[i]->qvs,system[i]->dam);
					system[i]->state = system[i]->state_check(y_0,GlobalVars->global_params,system[i]->params,system[i]->qvs,system[i]->dam);
				}
*/
				if(system[i]->state_check != NULL)
					system[i]->state = system[i]->state_check(y_0,GlobalVars->global_params,system[i]->params,system[i]->qvs,system[i]->dam);
			}
			else if(GlobalVars->init_flag == 3)
			{
				v_copy(db_init_buffer[i],y_0);

				if(system[i]->state_check != NULL)
					system[i]->state = system[i]->state_check(y_0,GlobalVars->global_params,system[i]->params,system[i]->qvs,system[i]->dam);
			}
			else
			{
				printf("[%i]: Error: Bad initial condition flag %i.\n",my_rank,GlobalVars->init_flag);
				MPI_Abort(MPI_COMM_WORLD,1);
			}

			system[i]->list = Create_List(y_0,t_0,dim,GlobalVars->num_dense,system[i]->method->s,GlobalVars->iter_limit);
			system[i]->list->head->state = system[i]->state;
			//if(GlobalVars->assim_flag == 0)	system[i]->disk_iterations = 1;
			//else				system[i]->disk_iterations = 0;
			system[i]->disk_iterations = 0;

			//Discontinuity information
			if(system[i]->numparents > 0)
				system[i]->discont = (double*) malloc(GlobalVars->discont_size*sizeof(double));
			if(system[i]->c != NULL && my_rank != (*assignments)[system[i]->c->location])
			{
				system[i]->discont_send = (double*) malloc(GlobalVars->discont_size*sizeof(double));
				system[i]->discont_order_send = (unsigned int*) malloc(GlobalVars->discont_size*sizeof(int));
				system[i]->discont_send_count = 0;
			}

			//Setup most of the remaining data
			system[i]->last_t = t_0;
			system[i]->print_time = t_0;
			system[i]->current_iterations = 1;
			system[i]->next_save = t_0 - 1.0;
			system[i]->peak_time = t_0;
			system[i]->save_flag = 0;
			system[i]->peak_flag = 0;
			system[i]->peak_value = v_get(dim);
			v_copy(y_0,system[i]->peak_value);
			if(system[i]->numparents == 0)
				system[i]->ready = 1;
			else
				system[i]->ready = 0;
			//system[i]->rain = NULL;
			system[i]->iters_removed = 0;
			system[i]->steps_on_diff_proc = 1;	//Note: This won't be used if link isn't stored on another proc
			system[i]->forcing_values = (double*) malloc(GlobalVars->num_forcings*sizeof(double));
			system[i]->forcing_change_times = (double*) malloc(GlobalVars->num_forcings*sizeof(double));
			system[i]->forcing_indices = (unsigned int*) malloc(GlobalVars->num_forcings*sizeof(double));
			system[i]->rejected = 0;
			system[i]->last_eta = 1e10;
			system[i]->compute_J = 1;
			system[i]->compute_LU = 1;

			if(system[i]->method->exp_imp == 1)
			{
				system[i]->JMatrix = m_get(GlobalVars->dim,GlobalVars->dim);
				system[i]->CoefMat = m_get(GlobalVars->max_s*GlobalVars->dim,GlobalVars->max_s*GlobalVars->dim);
				system[i]->Z_i = malloc(GlobalVars->max_s*sizeof(VEC*));
				for(j=0;j<GlobalVars->max_s;j++)
					system[i]->Z_i[j] = v_get(GlobalVars->dim);
				system[i]->sol_diff = v_get(GlobalVars->dim);
				system[i]->h_old = -1.0;
				system[i]->value_old = -1.0;
			}
			else
			{
				system[i]->JMatrix = NULL;
				system[i]->CoefMat = NULL;
				system[i]->Z_i = NULL;
				system[i]->sol_diff = NULL;
			}

			//Set up the rain fall data if using .str or .ustr file
			system[i]->forcing_buff = (ForcingData**) malloc(GlobalVars->num_forcings*sizeof(ForcingData*));
			for(l=0;l<GlobalVars->num_forcings;l++)
			{
				if(forcings[l]->flag == 1)
				{
					if(!(GlobalVars->res_flag) || !(l == GlobalVars->res_forcing_idx) || system[i]->res)
					{
						//Read id
						fscanf(forcingfiles[l],"%i",&m);
						if(m != system[i]->ID)
						{
							printf("[%i]: Error: bad link id. Got %i, expected %u for forcing %u.\n",my_rank,m,system[i]->ID,l);
							MPI_Abort(MPI_COMM_WORLD,1);
						}

						//Read the number of times for the rain data for this link
						fscanf(forcingfiles[l],"%i",&m);
						m++;		//Increase this by one to add a "ceiling" term

						system[i]->forcing_buff[l] = (ForcingData*) malloc(sizeof(ForcingData));
						system[i]->forcing_buff[l]->rainfall = (double**) malloc(m*sizeof(double*));
						system[i]->forcing_buff[l]->n_times = m;
						for(j=0;j<m;j++)	system[i]->forcing_buff[l]->rainfall[j] = (double*) malloc(2*sizeof(double));
						m--;

						//Read in the storm data for this link
						for(j=0;j<m;j++)
						{
							double tempy;
							fscanf(forcingfiles[l],"%lf",&tempy);
							(system[i]->forcing_buff[l]->rainfall[j][0]) = tempy;
							fscanf(forcingfiles[l],"%lf",&(system[i]->forcing_buff[l]->rainfall[j][1]));
						}
						system[i]->forcing_buff[l]->rainfall[m][0] = GlobalVars->maxtime + 3.0;
						system[i]->forcing_buff[l]->rainfall[m][1] = -1;

						double rainfall_buffer = system[i]->forcing_buff[l]->rainfall[0][1];
						system[i]->forcing_values[l] = rainfall_buffer;
						system[i]->forcing_indices[l] = 0;
						for(j=1;j<system[i]->forcing_buff[l]->n_times;j++)
						{
							//if(rainfall_buffer != system[i]->forcing_buff[l]->rainfall[j][1])
							if(fabs(system[i]->forcing_buff[l]->rainfall[j][1] - rainfall_buffer) > 1e-8)
							{
								system[i]->forcing_change_times[l] = system[i]->forcing_buff[l]->rainfall[j][0];
								break;
							}
						}
						if(j == system[i]->forcing_buff[l]->n_times)
						{
							system[i]->forcing_change_times[l] = system[i]->forcing_buff[l]->rainfall[j-1][0];
							system[i]->forcing_indices[l] = j-1;
						}
					}
					else	//Reservoir, so allocate only a little memory
					{
						m = 2;	//Init value (assumed 0.0)
						
						system[i]->forcing_buff[l] = (ForcingData*) malloc(sizeof(ForcingData));
						system[i]->forcing_buff[l]->rainfall = (double**) malloc(m*sizeof(double*));
						system[i]->forcing_buff[l]->n_times = m;
						for(j=0;j<m;j++)	system[i]->forcing_buff[l]->rainfall[j] = (double*) malloc(2*sizeof(double));

						system[i]->forcing_buff[l]->rainfall[0][0] = t_0;
						system[i]->forcing_buff[l]->rainfall[0][1] = 0.0;
						system[i]->forcing_buff[l]->rainfall[1][0] = GlobalVars->maxtime + 3.0;
						system[i]->forcing_buff[l]->rainfall[1][1] = -1.0;

						double rainfall_buffer = system[i]->forcing_buff[l]->rainfall[0][1];
						system[i]->forcing_values[l] = rainfall_buffer;
						system[i]->forcing_indices[l] = 0;
						for(j=1;j<system[i]->forcing_buff[l]->n_times;j++)
						{
							//if(rainfall_buffer != system[i]->forcing_buff[l]->rainfall[j][1])
							if(fabs(system[i]->forcing_buff[l]->rainfall[j][1] - rainfall_buffer) > 1e-8)
							{
								system[i]->forcing_change_times[l] = system[i]->forcing_buff[l]->rainfall[j][0];
								break;
							}
						}
						if(j == system[i]->forcing_buff[l]->n_times)
							system[i]->forcing_change_times[l] = system[i]->forcing_buff[l]->rainfall[j-1][0];
					}
				}
				else if(forcings[l]->flag == 4 || forcings[l]->flag == 7)	//Copy all the data from GlobalForcing
				{
					system[i]->forcing_buff[l] = forcings[l]->GlobalForcing;
					system[i]->forcing_values[l] = system[i]->forcing_buff[l]->rainfall[0][1];
					system[i]->forcing_change_times[l] = univ_forcing_change_time[l];
					system[i]->forcing_indices[l] = 0;
				}
				else if(forcings[l]->flag == 3 || forcings[l]->flag == 5)
				{
					if(!(GlobalVars->res_flag) || !(l == GlobalVars->res_forcing_idx) || system[i]->res)
					{
						m = forcings[l]->increment + 4;	//+1 for init, +1 for ceiling, +2 for when init time doesn't line up with file_time
						system[i]->forcing_buff[l] = (ForcingData*) malloc(sizeof(ForcingData));
						system[i]->forcing_buff[l]->rainfall = (double**) malloc(m*sizeof(double*));
						system[i]->forcing_buff[l]->n_times = m;
						for(j=0;j<m;j++)	system[i]->forcing_buff[l]->rainfall[j] = (double*) calloc(2,sizeof(double));
							system[i]->forcing_buff[l]->rainfall[0][0] = t_0;
						system[i]->forcing_values[l] = 0.0;
						system[i]->forcing_change_times[l] = fabs(t_0 + GlobalVars->maxtime) + 10.0;	//Just pick something away from t_0, and positive
					}
					else	//Reservoir, so allocate only a little memory
					{
						m = 4;	//+1 for init, +1 for ceiling, +2 for when init time doesn't line up with file_time
						system[i]->forcing_buff[l] = (ForcingData*) malloc(sizeof(ForcingData));
						system[i]->forcing_buff[l]->rainfall = (double**) malloc(m*sizeof(double*));
						system[i]->forcing_buff[l]->n_times = m;
						for(j=0;j<m;j++)	system[i]->forcing_buff[l]->rainfall[j] = (double*) calloc(2,sizeof(double));
							system[i]->forcing_buff[l]->rainfall[0][0] = t_0;
						system[i]->forcing_values[l] = 0.0;
						system[i]->forcing_change_times[l] = fabs(t_0 + GlobalVars->maxtime) + 10.0;	//Just pick something away from t_0, and positive
					}
				}
				else if(forcings[l]->flag == 0)
				{
					system[i]->forcing_buff[l] = NULL;
					system[i]->forcing_values[l] = 0.0;
					system[i]->forcing_change_times[l] = GlobalVars->maxtime + 1.0;
				}
				else
				{
					system[i]->forcing_buff[l] = NULL;
				}
			}
		}
		else	//If link not assigned to this process, and not receiving any information from this link.
		{
			if(rk_filename[0] != '\0')
			{
				for(j=0;j<dim;j++)	fscanf(rkdata,"%*f");	//Errors
				for(j=0;j<dim;j++)	fscanf(rkdata,"%*f");
				for(j=0;j<dim;j++)	fscanf(rkdata,"%*f");
				for(j=0;j<dim;j++)	fscanf(rkdata,"%*f");
/*
				for(j=0;j<rkdim;j++)	fscanf(rkdata,"%lf",&crap1);	//Errors
				for(j=0;j<rkdim;j++)	fscanf(rkdata,"%lf",&crap1);
				for(j=0;j<rkdim;j++)	fscanf(rkdata,"%lf",&crap1);
				for(j=0;j<rkdim;j++)	fscanf(rkdata,"%lf",&crap1);
*/

				fscanf(rkdata,"%i",&crap2);	//NA method
				fscanf(rkdata,"%lf",&crap1);	//Stepsize
			}

			//if(GlobalVars->init_flag == 0)	for(j=0;j<initdim;j++)	fscanf(initdata,"%lf",&crap1);	//Initial data
			if(GlobalVars->init_flag == 0)	for(j=GlobalVars->diff_start;j<GlobalVars->no_ini_start;j++)	fscanf(initdata,"%*f");	//Initial data
			else if(GlobalVars->init_flag == 2)	for(j=0;j<dim;j++)	fscanf(initdata,"%lf",&crap1);
			system[i]->method = NULL;
			system[i]->list = NULL;
			system[i]->peak_value = NULL;
			system[i]->forcing_buff = NULL;

			v_free(system[i]->params);
			system[i]->params = NULL;

			for(l=0;l<GlobalVars->num_forcings;l++)
			{
				if(forcings[l]->flag == 1 && (!(GlobalVars->res_flag) || !(l == GlobalVars->res_forcing_idx) || system[i]->res)) //Rainfall with .str
				{
					fscanf(forcingfiles[l],"%*i %i",&m);
					for(j=0;j<m;j++)
						fscanf(forcingfiles[l],"%*f %*f");
				}
			}

			//NULL out the unneeded data
			system[i]->f = NULL;
			system[i]->Jacobian = NULL;
			system[i]->JMatrix = NULL;
			system[i]->CoefMat = NULL;
			system[i]->Z_i = NULL;
			system[i]->sol_diff = NULL;
		}
	}

	v_free(y_0);

	MPI_Barrier(MPI_COMM_WORLD);
	stop = time(NULL);
	if(my_rank == 0)	printf("Time for second read: %f\n",difftime(stop,start));

	//Reserve workspace
	*workspace = Create_Workspace(GlobalVars->dim,GlobalVars->max_s,GlobalVars->max_parents);
/*
	//Prepare the equations for the standardized template models
	if(GlobalVars->template_flag)
	{
		PrepareParser(system,*N,*my_sys,*my_N,GlobalVars);
	}
*/
	//Calculate initial step sizes
	//if(GlobalVars->rain_flag != 2 && GlobalVars->rain_flag != 3 && GlobalVars->rain_flag != 5 && GlobalVars->type != 200)
	int calc_h = 1;
	for(i=0;i<GlobalVars->num_forcings;i++)	calc_h &= (system[(*my_sys)[0]]->forcing_buff[i] != NULL);
	if(calc_h)
	{
		if(my_rank == 0)	printf("Calculating step sizes...\n");
		MPI_Barrier(MPI_COMM_WORLD);
		start = time(NULL);

		for(i=0;i<*my_N;i++)
			system[(*my_sys)[i]]->h = InitialStepSize(t_0,system[(*my_sys)[i]],GlobalVars,*workspace);
		for(i=0;i<*my_N;i++)
		{
			for(l=0;l<GlobalVars->num_forcings;l++)
				if(system[(*my_sys)[i]]->forcing_buff[l])
					system[(*my_sys)[i]]->h = min(system[(*my_sys)[i]]->h,system[(*my_sys)[i]]->forcing_change_times[l] - system[(*my_sys)[i]]->last_t);
		}


		MPI_Barrier(MPI_COMM_WORLD);
		stop = time(NULL);
		if(my_rank == 0)	printf("Time for step sizes: %f\n",difftime(stop,start));
	}

	//Read in the list of links for which data will be saved and store it
	if(my_rank == 0)	printf("Building save lists...\n");
	MPI_Barrier(MPI_COMM_WORLD);
	start = time(NULL);

	*my_save_size = 0;
	if(GlobalVars->hydros_loc_flag)
		*save_list = Create_SAV_Data(GlobalVars->hydrosave_filename,system,*N,save_size,db_connections[ASYNCH_DB_LOC_HYDROSAVE],GlobalVars->hydrosave_flag);
	else
	{
		*save_list = NULL;
		*save_size = 0;
	}

	for(j=0;j<*save_size;j++)
	{
		k = find_link_by_idtoloc((*save_list)[j],*id_to_loc,*N);

		if(k == *N+1)
		{
			if(my_rank == 0)
				printf("[%i]: Warning: Time series output for requested for a link NOT in the network (link id = %u).\n",my_rank,(*save_list)[j]);

			//Shift ids in save list to spot k, and decrement save_size
			for(k=j+1;k<*save_size;k++)
				(*save_list)[k-1] = (*save_list)[k];
			j--;
			(*save_size)--;

			//(*save_list)[j] = (*save_list)[--(*save_size)];
			//j--;
		}
		//else if((*assignments)[(*id_to_loc)[k][1]] == my_rank)
		else if((*assignments)[k] == my_rank)
		{
			(*my_save_size)++;
			//current = system[(*id_to_loc)[k][1]];
			current = system[k];

			if(GlobalVars->print_time > 0.0)
				current->print_time = GlobalVars->print_time;
			else
				current->print_time = pow(current->params->ve[GlobalVars->area_idx]*.1,.5);

			//current->next_save = t_0 + current->print_time;
			current->next_save = t_0;
			current->save_flag = 1;
		}
	}


	*peaksave_size = 0;
	unsigned int* peaksave_list;
	if(GlobalVars->peaks_loc_flag)
	{
		peaksave_list = Create_SAV_Data(GlobalVars->peaksave_filename,system,*N,peaksave_size,db_connections[ASYNCH_DB_LOC_PEAKSAVE],GlobalVars->peaksave_flag);

		for(j=0;j<*peaksave_size;j++)		//!!!! Can we loop over my_peaksave_size? !!!!
		{
			k = find_link_by_idtoloc(peaksave_list[j],*id_to_loc,*N);

			if(k == *N+1)
			{
				if(my_rank == 0)
					printf("[%i]: Warning: Peakflow output for requested for a link NOT in the network (link id = %u).\n",my_rank,(peaksave_list)[j]);

				//Shift ids in save list to spot k, and decrement save_size
				for(k=j+1;k<*peaksave_size;k++)
					peaksave_list[k-1] = peaksave_list[k];
				j--;
				(*peaksave_size)--;

				//peaksave_list[j] = peaksave_list[--(*peaksave_size)];
				//j--;
			}
			else if((*assignments)[k] == my_rank)
			{
				//my_peaksave_size++;
				current = system[(*id_to_loc)[k][1]];
				current->peak_flag = 1;
			}
		}
	}
	else
	{
		peaksave_list = NULL;
		*peaksave_size = 0;
	}

	if(peaksave_list != NULL)	free(peaksave_list);


	//Find init values at reservoirs
	if(GlobalVars->res_flag)
	{
		//Download forcing data	!!!! Not sure if this is the way to go. Maybe separate function? !!!!
		forcings[GlobalVars->res_forcing_idx]->passes = 1;
		unsigned int start_iteration = forcings[GlobalVars->res_forcing_idx]->iteration;
		forcings[GlobalVars->res_forcing_idx]->GetNextForcing(system,*N,*my_sys,*my_N,*assignments,GlobalVars,forcings[GlobalVars->res_forcing_idx],db_connections,*id_to_loc,GlobalVars->res_forcing_idx);
		forcings[GlobalVars->res_forcing_idx]->iteration = start_iteration;	//Keep this the same

		//Setup links with forcing
		for(j=0;j<res_size;j++)
		{
			loc = find_link_by_idtoloc(res_list[j],*id_to_loc,*N);
			if(loc < *N && (*assignments)[loc] == my_rank)
			{
				//system[loc]->res = 1;

				//!!!! Not sure if this the way to go... !!!!
				system[loc]->f(system[loc]->last_t,system[loc]->list->tail->y_approx,NULL,system[loc]->numparents,GlobalVars->global_params,system[loc]->forcing_values,system[loc]->qvs,system[loc]->params,system[loc]->iparams,system[loc]->state,system[loc]->upstream,system[loc]->numupstream,system[loc]->list->tail->y_approx);
			}
			//else
			//{
			//	if(my_rank == 0)
			//		printf("[%i]: Warning: ignoring reservoir at ID %u. ID not found in network.\n",my_rank,res_list[j]);
			//}
		}

		free(res_list);
	}

	MPI_Barrier(MPI_COMM_WORLD);
	stop = time(NULL);
	if(my_rank == 0)	printf("Time for save lists: %f\n",difftime(stop,start));

	//Cleanup
	if(db_init_buffer)
	{
		for(i=0;i<*N;i++)
			if(db_init_buffer[i])	free(db_init_buffer[i]);
		free(db_init_buffer);
	}
	free(stack);
	if(riverdata != NULL)		fclose(riverdata);
	if(rk_filename[0] != '\0')	fclose(rkdata);
	if(initdata)	fclose(initdata);
	if(paramdata != NULL)		fclose(paramdata);
	//if(GlobalInit != NULL)		free(GlobalInit);
	for(l=0;l<GlobalVars->num_forcings;l++)
		if(forcings[l]->flag == 1)	fclose(forcingfiles[l]);
	free(forcingfiles);
	if(uniform_y_0 != NULL)	v_free(uniform_y_0);
	//free(leaves);
	//free(id_tracker);
	//if(peaksave_list != NULL)	free(peaksave_list);
	//free(upstream_order);
	//v_free(y_0);

	return system;
}


//Reads in the data from a .gbl file. All data will be global data for the entire river system.
//char globalfilename[]: String with the filename of the .gbl file.
//ErrorData** GlobalErrors (set by this method): Will contain the error data for the entire river system, if the error data is global.
//PGconn** conn:	NULL pointer that will be set to an SQL database, if needed.
//char* rkdfilename (set by this method): Will be the filename of the .rkd file, if the error data is not global.
//Returns a UnivVars that contains all the global data read in from the file globalfilename.
UnivVars* Read_Global_Data(char globalfilename[],ErrorData** GlobalErrors,Forcing** forcings,ConnData** db_connections,char* rkdfilename,model* custom_model,void* external)
{
	unsigned int i,j,k,total,written;
	int flag,valsread;
	char endmark;
	unsigned int string_size = 256;	//This is the max size of all strings for filenames and location
	unsigned int buff_size = string_size + 20;
	char* linebuffer = (char*) malloc(buff_size*sizeof(char));
	UnivVars* GlobalVars = (UnivVars*) malloc(sizeof(UnivVars));
	GlobalVars->string_size = string_size;
	char* db_filename = (char*) malloc(string_size*sizeof(char));
	GlobalVars->query_size = 1024;

	*GlobalErrors = malloc(sizeof(ErrorData));
	FILE* globalfile = fopen(globalfilename,"r");
	if(globalfile == NULL)
	{
		printf("Error: Global file %s was not found.\n",globalfilename);
		return NULL;
	}
	GlobalVars->rain_filename = NULL;

	//Grab the type and maxtime
	ReadLineFromTextFile(globalfile,linebuffer,buff_size,string_size);
	valsread = sscanf(linebuffer,"%hu %lf",&(GlobalVars->type),&(GlobalVars->maxtime));
	if(ReadLineError(valsread,2,"type and maxtime"))	return NULL;

	//Grab the output filename info
	ReadLineFromTextFile(globalfile,linebuffer,buff_size,string_size);
	valsread = sscanf(linebuffer,"%hu",&(GlobalVars->print_par_flag));
	if(ReadLineError(valsread,1,"to print filename parameters"))	return NULL;

	//Grab components to print
	ReadLineFromTextFile(globalfile,linebuffer,buff_size,string_size);
	valsread = sscanf(linebuffer,"%u",&(GlobalVars->num_print));
	if(ReadLineError(valsread,1,"number of indices to print"))	return NULL;
	GlobalVars->output_names = (char**) malloc(GlobalVars->num_print*sizeof(char*));
	for(i=0;i<GlobalVars->num_print;i++)
	{
		GlobalVars->output_names[i] = (char*) malloc(string_size*sizeof(char));
		ReadLineFromTextFile(globalfile,linebuffer,buff_size,string_size);
		valsread = sscanf(linebuffer,"%s",GlobalVars->output_names[i]);
		if(ReadLineError(valsread,1,"a component to print"))	return NULL;
	}

	//Peakflow function
	GlobalVars->peakflow_function_name = (char*) malloc(string_size*sizeof(char));
	ReadLineFromTextFile(globalfile,linebuffer,buff_size,string_size);
	valsread = sscanf(linebuffer,"%s",GlobalVars->peakflow_function_name);
	if(ReadLineError(valsread,1,"peakflow function name"))	return NULL;
	SetPeakflowOutputFunctions(GlobalVars->peakflow_function_name,&(GlobalVars->peakflow_output));

	//Grab the parameters
	ReadLineFromTextFile(globalfile,linebuffer,buff_size,string_size);
	valsread = sscanf(linebuffer,"%u%n",&i,&total);
	if(ReadLineError(valsread,1,"number of global parameters"))	return NULL;
	GlobalVars->global_params = v_get(i);
	for(i=0;i<GlobalVars->global_params->dim;i++)
	{
		valsread = sscanf(&(linebuffer[total]),"%lf%n",&(GlobalVars->global_params->ve[i]),&written);
		if(ReadLineError(valsread,1,"a global parameter"))	return NULL;
		total += written;
	}

	//Set dim and other sizes
	if(custom_model)	custom_model->SetParamSizes(GlobalVars,external);
	else			SetParamSizes(GlobalVars,external);

	//Find the states needed for printing
	GlobalVars->print_indices = (unsigned int*) calloc(GlobalVars->dim,sizeof(unsigned int));
	GlobalVars->outputs_i = (int (**)(double,VEC*,VEC*,VEC*,IVEC*,int,void*)) calloc( GlobalVars->num_print, sizeof(int (*)(double,VEC*,VEC*,VEC*,IVEC*,int,void*)) );
	GlobalVars->outputs_d = (double (**)(double,VEC*,VEC*,VEC*,IVEC*,int,void*)) calloc( GlobalVars->num_print, sizeof(double (*)(double,VEC*,VEC*,VEC*,IVEC*,int,void*)) );
	GlobalVars->output_types = (short int*) malloc(GlobalVars->num_print*sizeof(short int));
	//GlobalVars->output_names = (char**) malloc(GlobalVars->num_print*sizeof(char*));
	GlobalVars->output_sizes = (short int*) malloc(GlobalVars->num_print*sizeof(short int));
	GlobalVars->output_specifiers = (char**) malloc(GlobalVars->num_print*sizeof(char*));
	for(i=0;i<GlobalVars->num_print;i++)
	{
		GlobalVars->output_specifiers[i] = (char*) malloc(16*sizeof(char));
		SetOutputFunctions(GlobalVars->output_names[i],GlobalVars->output_specifiers[i],GlobalVars->num_print,GlobalVars->print_indices,&(GlobalVars->output_sizes[i]),&(GlobalVars->output_types[i]),&(GlobalVars->outputs_i[i]),&(GlobalVars->outputs_d[i]));
	}

	//Adjust dense_indicies to include print indices
	//Check first for redundancy
	for(i=0;i<GlobalVars->dim;i++)
	{
		if(GlobalVars->print_indices[i])
		{
			for(j=0;j<GlobalVars->num_dense;j++)
			{
				if(GlobalVars->dense_indices[j] == i)
				{
					for(k=j+1;k<GlobalVars->num_dense;k++)
						GlobalVars->dense_indices[k-1] = GlobalVars->dense_indices[k];
					(GlobalVars->num_dense)--;
				}
			}
		}
	}

	//Find how many states are to be printed
	unsigned int total_states_to_print = 0;
	for(i=0;i<GlobalVars->dim;i++)	total_states_to_print += GlobalVars->print_indices[i];

	//Add print_indices to the dense list
	GlobalVars->dense_indices = realloc(GlobalVars->dense_indices,(GlobalVars->num_dense + GlobalVars->num_print)*sizeof(unsigned int));
	j = 0;
	for(i=GlobalVars->num_dense;i<GlobalVars->num_dense + GlobalVars->num_print;i++)
	{
		while(j < GlobalVars->dim && !(GlobalVars->print_indices[j]))	j++;
		GlobalVars->dense_indices[i] = j;
		j++;
	}
	//GlobalVars->num_dense += GlobalVars->num_print;
	GlobalVars->num_dense += total_states_to_print;
	merge_sort_1D(GlobalVars->dense_indices,GlobalVars->num_dense);

	//Grab the stored steps limits
	ReadLineFromTextFile(globalfile,linebuffer,buff_size,string_size);
	valsread = sscanf(linebuffer,"%u %i %u",&(GlobalVars->iter_limit),&(GlobalVars->max_transfer_steps),&(GlobalVars->discont_size));
	if(ReadLineError(valsread,3,"steps stored, steps transfered, and discontinuity buffer size"))	return NULL;

	//Grab the topology data filename
	GlobalVars->outletlink = 0;
	ReadLineFromTextFile(globalfile,linebuffer,buff_size,string_size);
	valsread = sscanf(linebuffer,"%hu",&(GlobalVars->rvr_flag));
	if(ReadLineError(valsread,1,"topology data flag"))	return NULL;
	if(GlobalVars->rvr_flag == 0)
	{
		GlobalVars->rvr_filename = (char*) malloc(string_size*sizeof(char));
		valsread = sscanf(linebuffer,"%*u %s",GlobalVars->rvr_filename);
		if(ReadLineError(valsread,1,"filename for topology data"))	return NULL;
	}
	else	//Reading from database
	{
		valsread = sscanf(linebuffer,"%*u %u %s",&(GlobalVars->outletlink),db_filename);
		if(ReadLineError(valsread,2,"link id of downstream link for topology data or .dbc for topology"))	return NULL;
		GlobalVars->rvr_filename = NULL;
		db_connections[ASYNCH_DB_LOC_TOPO] = ReadDBC(db_filename,string_size);
		if(!db_connections[ASYNCH_DB_LOC_TOPO])	return NULL;
	}

	//Grab the parameter data filename
	ReadLineFromTextFile(globalfile,linebuffer,buff_size,string_size);
	valsread = sscanf(linebuffer,"%hu",&(GlobalVars->prm_flag));
	if(ReadLineError(valsread,1,"parameter flag"))	return NULL;
	if(GlobalVars->prm_flag == 0)
	{
		GlobalVars->prm_filename = (char*) malloc(string_size*sizeof(char));
		valsread = sscanf(linebuffer,"%*u %s",GlobalVars->prm_filename);
		if(ReadLineError(valsread,1,".prm filename"))	return NULL;
	}
	else
	{
		valsread = sscanf(linebuffer,"%*u %s",db_filename);
		if(ReadLineError(valsread,1,".dbc for parameters"))	return NULL;
		GlobalVars->prm_filename = NULL;
		db_connections[ASYNCH_DB_LOC_PARAMS] = ReadDBC(db_filename,string_size);
		if(!db_connections[ASYNCH_DB_LOC_PARAMS])	return NULL;
	}

	//Grab the initial data file
	GlobalVars->init_filename = NULL;
	ReadLineFromTextFile(globalfile,linebuffer,buff_size,string_size);
	valsread = sscanf(linebuffer,"%hu",&(GlobalVars->init_flag));
	if(ReadLineError(valsread,1,"initial data flag"))	return NULL;
	if(GlobalVars->init_flag == 0 || GlobalVars->init_flag == 1 || GlobalVars->init_flag == 2)
	{
		GlobalVars->init_filename = (char*) malloc(string_size*sizeof(char));
		valsread = sscanf(linebuffer,"%*u %s",GlobalVars->init_filename);
		if(ReadLineError(valsread,1,"initial data flag"))	return NULL;
	}
	else if(GlobalVars->init_flag == 3)
	{
		valsread = sscanf(linebuffer,"%*u %s %u",db_filename,&(GlobalVars->init_timestamp));
		if(ReadLineError(valsread,1,".dbc for parameters"))	return NULL;
		db_connections[ASYNCH_DB_LOC_INIT] = ReadDBC(db_filename,string_size);
		if(!db_connections[ASYNCH_DB_LOC_INIT])	return NULL;
	}

	//Grab number of forcings
	unsigned int got_forcings;
	ReadLineFromTextFile(globalfile,linebuffer,buff_size,string_size);
	valsread = sscanf(linebuffer,"%u",&got_forcings);
	if(ReadLineError(valsread,1,"rainfall flag"))	return NULL;
	if(got_forcings < GlobalVars->num_forcings && my_rank == 0)
	{
		printf("[%i]: Error: Got %u forcings in the .gbl file. Expected %u for model %u.\n",my_rank,got_forcings,GlobalVars->num_forcings,GlobalVars->type);
		MPI_Abort(MPI_COMM_WORLD,1);
	}
	if(got_forcings > GlobalVars->num_forcings && my_rank == 0)
	{
		printf("[%i]: Warning: Got %u forcings in the .gbl file. Expected %u for model %u.\n",my_rank,got_forcings,GlobalVars->num_forcings,GlobalVars->type);
		GlobalVars->num_forcings = got_forcings;
	}

	//Grab the forcing parameters
	//0 for no rain, 1 for .str file, 2 for binary files, 3 for database, 4 for uniform rain (.ustr)
	GlobalVars->hydro_table = GlobalVars->peak_table = NULL;
	for(i=0;i<GlobalVars->num_forcings;i++)
	{
		ReadLineFromTextFile(globalfile,linebuffer,buff_size,string_size);
		valsread = sscanf(linebuffer,"%hi",&(forcings[i]->flag));
		if(ReadLineError(valsread,1,"forcings flag"))	return NULL;

		if(forcings[i]->flag == 1 || forcings[i]->flag == 2 || forcings[i]->flag == 4 || forcings[i]->flag == 6)
		{
			forcings[i]->filename = (char*) malloc(string_size*sizeof(char));
			valsread = sscanf(linebuffer,"%*i %s",forcings[i]->filename);
			if(ReadLineError(valsread,1,"forcing data filename"))	return NULL;
			if(forcings[i]->flag == 2 || forcings[i]->flag == 6)
			{
				ReadLineFromTextFile(globalfile,linebuffer,buff_size,string_size);
				valsread = sscanf(linebuffer,"%u %lf %u %u",&(forcings[i]->increment),&(forcings[i]->file_time),&(forcings[i]->first_file),&(forcings[i]->last_file));
				if(ReadLineError(valsread,4,"time increment, file time, first file, and last file"))	return NULL;
			}
		}
		//else if(forcings[i]->flag == 3 || forcings[i]->flag == 5)	//Database
		else if(forcings[i]->flag == 3)	//Database
		{
			valsread = sscanf(linebuffer,"%*i %s",db_filename);
			if(ReadLineError(valsread,1,".dbc for rainfall"))	return NULL;
			db_connections[ASYNCH_DB_LOC_FORCING_START+i] = ReadDBC(db_filename,string_size);

			ReadLineFromTextFile(globalfile,linebuffer,buff_size,string_size);
			valsread = sscanf(linebuffer,"%u %lf %u %u",&(forcings[i]->increment),&(forcings[i]->file_time),&(forcings[i]->first_file),&(forcings[i]->last_file));
			if(ReadLineError(valsread,4,"time increment, file time, first file, and last file"))	return NULL;
			forcings[i]->raindb_start_time = forcings[i]->first_file;
/*
			if(forcings[i]->flag == 5)	//Flood forecaster	!!!! This should not exist !!!!
			{
				if(GlobalVars->hydro_table)	printf("[%i]: Warning: Only one forcing of type 5 (forecasting) should be selected.\n",my_rank);
				GlobalVars->halt_filename = (char*) malloc(string_size*sizeof(char));

				ReadLineFromTextFile(globalfile,linebuffer,buff_size,string_size);
				valsread = sscanf(linebuffer,"%u",&(forcings[i]->num_rainsteps));
				if(ReadLineError(valsread,1,"number of rainfall steps"))	return NULL;
				ReadLineFromTextFile(globalfile,linebuffer,buff_size,string_size);
				valsread = sscanf(linebuffer,"%s",GlobalVars->halt_filename);
				if(ReadLineError(valsread,1,"halt filename"))	return NULL;
			}
*/
		}
		else if(forcings[i]->flag == 7)	//Recurring
		{
			forcings[i]->filename = (char*) malloc(string_size*sizeof(char));
			valsread = sscanf(linebuffer,"%*i %s",forcings[i]->filename);
			if(ReadLineError(valsread,1,"recurring rainfall filename"))	return NULL;

			ReadLineFromTextFile(globalfile,linebuffer,buff_size,string_size);
			valsread = sscanf(linebuffer,"%u %u",&(forcings[i]->first_file),&(forcings[i]->last_file));
			if(ReadLineError(valsread,2,"first time, and last time"))	return NULL;
			forcings[i]->raindb_start_time = forcings[i]->first_file;
		}
		else if(forcings[i]->flag == 0) //No forcing
		{
			forcings[i]->filename = NULL;
		}
		else
		{
			printf("[%i]: Error reading %s: Invalid forcing flag %i.\n",my_rank,globalfilename,forcings[i]->flag);
			return NULL;
		}
	}

	//Grab the .dam filename
	ReadLineFromTextFile(globalfile,linebuffer,buff_size,string_size);
	valsread = sscanf(linebuffer,"%hu",&(GlobalVars->dam_flag));
	if(ReadLineError(valsread,1,"dam flag"))	return NULL;

	if(GlobalVars->dam_flag)
	{
		GlobalVars->dam_filename = (char*) malloc(string_size*sizeof(char));
		valsread = sscanf(linebuffer,"%*u %s",GlobalVars->dam_filename);
		if(ReadLineError(valsread,1,"filename for dam info"))	return NULL;
		if(GlobalVars->dam_flag == 3)
			db_connections[ASYNCH_DB_LOC_QVS] = ReadDBC(GlobalVars->dam_filename,string_size);
	}
	else
		GlobalVars->dam_filename = NULL;

	//Get the link ids where reservoirs exist
	ReadLineFromTextFile(globalfile,linebuffer,buff_size,string_size);
	valsread = sscanf(linebuffer,"%hu",&(GlobalVars->res_flag));
	if(ReadLineError(valsread,1,"res flag"))	return NULL;

	if(GlobalVars->res_flag)
	{
		if(GlobalVars->res_flag == 1)
		{
			GlobalVars->rsv_filename = (char*) malloc(string_size*sizeof(char));
			valsread = sscanf(linebuffer,"%*u %s %hi",GlobalVars->rsv_filename,&(GlobalVars->res_forcing_idx));
			if(ReadLineError(valsread,2,".rsv filename"))	return NULL;
		}
		else	//Flag is 2
		{
			GlobalVars->rsv_filename = (char*) malloc(string_size*sizeof(char));
			valsread = sscanf(linebuffer,"%*u %s %hi",GlobalVars->rsv_filename,&(GlobalVars->res_forcing_idx));
			if(ReadLineError(valsread,2,".dbc for reservoirs"))	return NULL;
			db_connections[ASYNCH_DB_LOC_RSV] = ReadDBC(GlobalVars->rsv_filename,string_size);
		}

		if(GlobalVars->res_forcing_idx >= GlobalVars->num_forcings)
		{
			printf("[%i]: Bad forcing index for a reservoir feed (%hi). Only %i forcings available.\n",my_rank,GlobalVars->res_forcing_idx,GlobalVars->num_forcings);
			return NULL;
		}
	}
	else
	{
		GlobalVars->rsv_filename = NULL;
		GlobalVars->res_forcing_idx = -1;
	}

	//Grab where to write the hydrographs
	ReadLineFromTextFile(globalfile,linebuffer,buff_size,string_size);
	valsread = sscanf(linebuffer,"%hu",&(GlobalVars->hydros_loc_flag));
	if(ReadLineError(valsread,1,"hydrographs location"))	return NULL;

	//GlobalVars->results_folder = NULL;
	GlobalVars->hydros_loc_filename = NULL;
	GlobalVars->hydro_table = NULL;

	if(GlobalVars->hydros_loc_flag == 1 || GlobalVars->hydros_loc_flag == 2)
	{
		GlobalVars->hydros_loc_filename = (char*) malloc(string_size*sizeof(char));
		valsread = sscanf(linebuffer,"%*u %lf %s",&(GlobalVars->print_time),GlobalVars->hydros_loc_filename);
		if(ReadLineError(valsread,2,"hydrographs location"))	return NULL;
		GlobalVars->output_flag = (GlobalVars->hydros_loc_flag == 1) ? 0 : 1;

		if(GlobalVars->hydros_loc_flag == 1)	RemoveSuffix(GlobalVars->hydros_loc_filename,".dat");
		else if(GlobalVars->hydros_loc_flag == 2)	RemoveSuffix(GlobalVars->hydros_loc_filename,".csv");
	}
	else if(GlobalVars->hydros_loc_flag == 3)
	{
		GlobalVars->hydros_loc_filename = (char*) malloc(string_size*sizeof(char));
		GlobalVars->hydro_table = (char*) malloc(string_size*sizeof(char));
		valsread = sscanf(linebuffer,"%*u %lf %s %s",&(GlobalVars->print_time),GlobalVars->hydros_loc_filename,GlobalVars->hydro_table);
		if(ReadLineError(valsread,3,"hydrographs location"))	return NULL;
		db_connections[ASYNCH_DB_LOC_HYDRO_OUTPUT] = ReadDBC(GlobalVars->hydros_loc_filename,string_size);
	}

	//Grab where to write the peakflow data
	ReadLineFromTextFile(globalfile,linebuffer,buff_size,string_size);
	valsread = sscanf(linebuffer,"%hu",&(GlobalVars->peaks_loc_flag));
	if(ReadLineError(valsread,1,"peakflow location"))	return NULL;

	if(GlobalVars->peaks_loc_flag == 0)
	{
		GlobalVars->peaks_loc_filename = NULL;
	}
	else if(GlobalVars->peaks_loc_flag == 1)
	{
		GlobalVars->peaks_loc_filename = (char*) malloc(string_size*sizeof(char));
		valsread = sscanf(linebuffer,"%*u %s",GlobalVars->peaks_loc_filename);
		if(ReadLineError(valsread,1,"peakflow location"))	return NULL;

		RemoveSuffix(GlobalVars->peaks_loc_filename,".pea");
	}
	else if(GlobalVars->peaks_loc_flag == 2)
	{
		GlobalVars->peaks_loc_filename = (char*) malloc(string_size*sizeof(char));
		GlobalVars->peak_table = (char*) malloc(string_size*sizeof(char));
		valsread = sscanf(linebuffer,"%*u %s %s",GlobalVars->peaks_loc_filename,GlobalVars->peak_table);
		if(ReadLineError(valsread,2,"peakflow location"))	return NULL;
		db_connections[ASYNCH_DB_LOC_PEAK_OUTPUT] = ReadDBC(GlobalVars->peaks_loc_filename,string_size);
	}

	//Grab the .sav files
	ReadLineFromTextFile(globalfile,linebuffer,buff_size,string_size);
	valsread = sscanf(linebuffer,"%hu",&(GlobalVars->hydrosave_flag));
	if(ReadLineError(valsread,1,"hydrographs save flag"))	return NULL;

	if(GlobalVars->hydrosave_flag == 1)
	{
		GlobalVars->hydrosave_filename = (char*) malloc(string_size*sizeof(char));
		valsread = sscanf(linebuffer,"%*u %s",GlobalVars->hydrosave_filename);
		if(ReadLineError(valsread,1,"hydrographs .sav filename"))	return NULL;
	}
	else if(GlobalVars->hydrosave_flag == 2)
	{
		valsread = sscanf(linebuffer,"%*u %s",db_filename);
		if(ReadLineError(valsread,1,".dbc for hydrograph save ids"))	return NULL;
		GlobalVars->hydrosave_filename = NULL;
		db_connections[ASYNCH_DB_LOC_HYDROSAVE] = ReadDBC(db_filename,string_size);
	}
	else
		GlobalVars->hydrosave_filename = NULL;

	ReadLineFromTextFile(globalfile,linebuffer,buff_size,string_size);
	valsread = sscanf(linebuffer,"%hu",&(GlobalVars->peaksave_flag));
	if(ReadLineError(valsread,1,"peakflows save flag"))	return NULL;

	if(GlobalVars->peaksave_flag == 1)
	{
		GlobalVars->peaksave_filename = (char*) malloc(string_size*sizeof(char));
		valsread = sscanf(linebuffer,"%*u %s",GlobalVars->peaksave_filename);
		if(ReadLineError(valsread,1,"peakflows .sav filename"))	return NULL;
	}
	else if(GlobalVars->peaksave_flag == 2)
	{
		valsread = sscanf(linebuffer,"%*u %s",db_filename);
		if(ReadLineError(valsread,1,".dbc for peakflow save ids"))	return NULL;
		GlobalVars->peaksave_filename = NULL;
		db_connections[ASYNCH_DB_LOC_PEAKSAVE] = ReadDBC(db_filename,string_size);
	}
	else
		GlobalVars->peaksave_filename = NULL;
	GlobalVars->peakfilename = NULL;

	//Grab data dump info
	ReadLineFromTextFile(globalfile,linebuffer,buff_size,string_size);
	valsread = sscanf(linebuffer,"%hu",&(GlobalVars->dump_loc_flag));
	if(ReadLineError(valsread,1,"snapshot save flag"))	return NULL;

	GlobalVars->dump_loc_filename = NULL;
	GlobalVars->dump_table = NULL;

	if(GlobalVars->dump_loc_flag == 1)
	{
		GlobalVars->dump_loc_filename = (char*) malloc(string_size*sizeof(char));
		valsread = sscanf(linebuffer,"%*u %s",GlobalVars->dump_loc_filename);
		if(ReadLineError(valsread,1,"snapshot filename"))	return NULL;
	}
	else if(GlobalVars->dump_loc_flag == 2)
	{
		GlobalVars->dump_table = (char*) malloc(string_size*sizeof(char));
		valsread = sscanf(linebuffer,"%*u %s %s",db_filename,GlobalVars->dump_table);
		if(ReadLineError(valsread,1,".dbc for snapshots"))	return NULL;
		GlobalVars->dump_loc_filename = NULL;
		db_connections[ASYNCH_DB_LOC_SNAPSHOT_OUTPUT] = ReadDBC(db_filename,string_size);
	}

	//Grab folder locations
	//GlobalVars->results_folder = (char*) malloc(string_size*sizeof(char));
	GlobalVars->temp_filename = (char*) malloc(string_size*sizeof(char));
	//ReadLineFromTextFile(globalfile,linebuffer,buff_size,string_size);
	//valsread = sscanf(linebuffer,"%s",GlobalVars->results_folder);
	//if(ReadLineError(valsread,1,"results folder"))	return NULL;
	ReadLineFromTextFile(globalfile,linebuffer,buff_size,string_size);
	valsread = sscanf(linebuffer,"%s",GlobalVars->temp_filename);
	if(ReadLineError(valsread,1,"scratch work folder"))	return NULL;

	if(GlobalVars->print_par_flag)	//!!!! Is this needed? Why bother? !!!!
	{
		if(AttachParameters(GlobalVars->temp_filename,string_size,GlobalVars->global_params,string_size))
		{
			printf("[%i]: Error attaching global parameters to temporary filenames.\n",my_rank);
			return NULL;
		}
	}

	sprintf(db_filename,"_%i_%i",getpid(),my_rank);
	strcat(GlobalVars->temp_filename,db_filename);

	//Grab adapative data
	ReadLineFromTextFile(globalfile,linebuffer,buff_size,string_size);
	valsread = sscanf(linebuffer,"%lf %lf %lf",&((*GlobalErrors)->facmin),&((*GlobalErrors)->facmax),&((*GlobalErrors)->fac));
	if(ReadLineError(valsread,3,"facmin, facmax, fac"))	return NULL;

	//Read in the flag for the error tolerances
	ReadLineFromTextFile(globalfile,linebuffer,buff_size,string_size);
	valsread = sscanf(linebuffer,"%i",&flag);
	if(ReadLineError(valsread,1,"error tolerance flag"))	return NULL;

	//Set some parameters
	GlobalVars->max_s = 0;
	GlobalVars->max_parents = 0;
	//GlobalVars->discont_size = 30;

	//Connect to SQL database if needed
/*
	if(GlobalVars->rvr_flag == 1 || GlobalVars->prm_flag == 1 || GlobalVars->rain_flag == 3 || GlobalVars->rain_flag == 5 || GlobalVars->hydrosave_flag == 2 || GlobalVars->peaksave_flag == 2)
	{
		if(my_rank == 0)	ConnectPGDB(conninfo);
	}
*/

	if(flag == 0)	//Error data is found in the universal file
	{
		rkdfilename[0] = '\0';
		ReadLineFromTextFile(globalfile,linebuffer,buff_size,string_size);
		valsread = sscanf(linebuffer,"%u",&flag);
		if(ReadLineError(valsread,1,"RK method index"))	return NULL;
		GlobalVars->method = flag;

		(*GlobalErrors)->abstol = v_get(GlobalVars->dim);
		(*GlobalErrors)->reltol = v_get(GlobalVars->dim);
		(*GlobalErrors)->abstol_dense = v_get(GlobalVars->dim);
		(*GlobalErrors)->reltol_dense = v_get(GlobalVars->dim);

		ReadLineFromTextFile(globalfile,linebuffer,buff_size,string_size);
		total = 0;
		for(i=0;i<GlobalVars->dim;i++)
		{
			valsread = sscanf(&(linebuffer[total]),"%lf%n",&((*GlobalErrors)->abstol->ve[i]),&written);
			if(ReadLineError(valsread,1,"an abstol component"))	return NULL;
			total += written;
		}

		ReadLineFromTextFile(globalfile,linebuffer,buff_size,string_size);
		total = 0;
		for(i=0;i<GlobalVars->dim;i++)
		{
			valsread = sscanf(&(linebuffer[total]),"%lf%n",&((*GlobalErrors)->reltol->ve[i]),&written);
			if(ReadLineError(valsread,1,"a reltol component"))	return NULL;
			total += written;
		}

		ReadLineFromTextFile(globalfile,linebuffer,buff_size,string_size);
		total = 0;
		for(i=0;i<GlobalVars->dim;i++)
		{
			valsread = sscanf((&linebuffer[total]),"%lf%n",&((*GlobalErrors)->abstol_dense->ve[i]),&written);
			if(ReadLineError(valsread,1,"an abstol dense output component"))	return NULL;
			total += written;
		}

		ReadLineFromTextFile(globalfile,linebuffer,buff_size,string_size);
		total = 0;
		for(i=0;i<GlobalVars->dim;i++)
		{
			valsread = sscanf(&(linebuffer[total]),"%lf%n",&((*GlobalErrors)->reltol_dense->ve[i]),&written);
			if(ReadLineError(valsread,1,"a reltol dense output component"))	return NULL;
			total += written;
		}
	}
	else if(flag == 1)	//Error data is in an .rkd file
	{
		GlobalVars->method = -1;
		(*GlobalErrors)->abstol = NULL;
		(*GlobalErrors)->reltol = NULL;
		(*GlobalErrors)->abstol_dense = NULL;
		(*GlobalErrors)->reltol_dense = NULL;

		//ReadLineFromTextFile(globalfile,linebuffer,buff_size,string_size);
		valsread = sscanf(linebuffer,"%*i %s",rkdfilename);
		if(ReadLineError(valsread,1,".rkd filename"))	return NULL;
	}

	//Check for end mark
	ReadLineFromTextFile(globalfile,linebuffer,buff_size,string_size);
	sscanf(linebuffer,"%c",&endmark);
	if(endmark != '#')
	{
		printf("Error: an ending # not seen in %s file. Got %c.\n",globalfilename,endmark);
		return NULL;
	}

	//Setup an io object
	GlobalVars->output_data = BuildIO(GlobalVars);

	//Clean up
	free(db_filename);
	free(linebuffer);
	fclose(globalfile);
	return GlobalVars;
}

//Reads in the contents of a .sav file.
//char filename[]: The filename of the .sav file.
//int N: The number of links in the river system.
//int* size (set by this method): Will be the number of links for which data must be written to disk (number of links in the .sav file).
//Returns a list of ids for which data is to be saved to disk.
unsigned int* Create_SAV_Data(char filename[],Link** sys,unsigned int N,unsigned int* size,ConnData *conninfo,unsigned short int flag)
{
	unsigned int i,id;
	unsigned int* save_list = NULL;
	*size = 0;

	if(flag == 1)
	{
		unsigned int* initial_list = NULL;
		FILE* save_file = fopen(filename,"r");
		if(!save_file)
		{
			printf("Error opening file %s\n",filename);
			*size = -1;
			MPI_Abort(MPI_COMM_WORLD,1);
		}
		initial_list = malloc(N*sizeof(unsigned int));

		//Read in the .sav file
		while( fscanf(save_file,"%u",&id) != EOF )
		{
			initial_list[*size] = id;
			(*size)++;
		}

		save_list = malloc(*size * sizeof(unsigned int));
		for(i=0;i<*size;i++)	save_list[i] = initial_list[i];		//!!!! Blah. Realloc. !!!!

		fclose(save_file);
		free(initial_list);

		//Sort the list (for use in outputting data)
		//merge_sort_1D(save_list,*size);
	}
	else if(flag == 2)	//Grab from database
	{
		//char* query = conninfo->query;
		PGresult *res;

		if(my_rank == 0)
		{
			ConnectPGDB(conninfo);
			sprintf(conninfo->query,conninfo->queries[0]);
			res = PQexec(conninfo->conn,conninfo->query);
			CheckResError(res,"locating links with sensors");

			*size = PQntuples(res);
			save_list = malloc(*size * sizeof(unsigned int));
			for(i=0;i<*size;i++)	save_list[i] = atoi(PQgetvalue(res,i,0));
			PQclear(res);
			DisconnectPGDB(conninfo);
		}

		MPI_Bcast(size,1,MPI_INT,0,MPI_COMM_WORLD);
		if(my_rank != 0)	save_list = malloc(*size * sizeof(unsigned int));
		MPI_Bcast(save_list,*size,MPI_INT,0,MPI_COMM_WORLD);
	}
	else if(flag == 3)	//All links
	{
		*size = N;
		save_list = malloc(*size * sizeof(unsigned int));
		for(i=0;i<N;i++)	save_list[i] = sys[i]->ID;
	}

	return save_list;
}


void ReadLineFromTextFile(FILE* globalfile,char* linebuffer,unsigned int size,unsigned int string_size)
{
	linebuffer[0] = '%';
	while(!feof(globalfile) && (linebuffer[0] == '%' || linebuffer[0] == '\n'))	fgets(linebuffer,size,globalfile);
	if(my_rank == 0 && (strlen(linebuffer) > string_size + 10) )	printf("Warning: %zu %u Line in .gbl file may be too long. Read in the long line:\n\"%s\"\n",strlen(linebuffer),string_size,linebuffer);
}

int ReadLineError(int valsread,int valswant,char message[])
{
	if(valsread < valswant)
	{
		if(my_rank == 0)	printf("Error: Did not get a value from .gbl file for %s. %i\n",message,valsread);
		return 1;
	}
	return 0;
}

//Removes a suffix from filename, if present.
//Returns 1 if suffix removed
//0 if not (not present)
int RemoveSuffix(char* filename,char suffix[])
{
	unsigned int filename_length = strlen(filename);
	unsigned int suffix_length = strlen(suffix);
	int i,j;

	if(suffix_length > filename_length)	return 0;

	i = filename_length-1;
	for(j=suffix_length-1;j>-1;j--)
		if(suffix[j] != filename[i--])	return 0;

	filename[i+1] = '\0';
	return 1;
}

//Put a vector of global_params onto the end of a filename.
//Returns 1 if filename is not long enough to support this.
//Retruns 0 if the parameters are attached.
int AttachParameters(char* filename,unsigned int max_size,VEC* v,unsigned int string_size)
{
	unsigned int i,count,total=0;
	char buffer[string_size];
	unsigned int length = strlen(filename);

	for(i=0;i<v->dim;i++)
	{
		sprintf(buffer,"_%.4e%n",v->ve[i],&count);
		total += count;
		if(count+1 > string_size)	return 1;
		if(total+1 > max_size)		return 1;
		strcat(filename,buffer);
	}

	return 0;
}





