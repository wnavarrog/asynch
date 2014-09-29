#include "processdata.h"

//Reads the results stored in temporary files and outputs them conveniently to a new file.
//Use this for parallel implementations.
//There could be problems if there is lots of data in the temporary files. Improve on this in the future.
//Link** sys: The data for the riversystem.
//int* my_sys: The index in sys of all links assigned to this process.
//int my_N: The number of links assigned to this process.
//int save_size: Number of links for which data is to be written to file. This is the grand total for all processes.
//int my_save_size: Number of links for which data is to be written to file. This is the number for just this process.
//Return value = 0 means everything is good.
//1 means a database related error.
//2 means a file system related error.
//!!!! Is additional_temp needed? !!!!
int Process_Data(Link** sys,UnivVars* GlobalVars,unsigned int N,unsigned int* save_list,unsigned int save_size,unsigned int my_save_size,unsigned int** id_to_loc,int* assignments,char* additional_temp,char* additional_out,ConnData* conninfo,FILE** my_tempfile)
{
	unsigned int blocks = 10;
	int i,proc,k;
	unsigned int j,l,m,loc,id,counter,total_spaces,max_disk,**space_counter;
	char filename[GlobalVars->string_size],filenamespace[GlobalVars->string_size],outputfilename[GlobalVars->string_size];
	double total_time;
	char data_storage[16];
	fpos_t **positions;
	unsigned int dim = GlobalVars->num_print;
	FILE *inputfile,*outputfile,**allinputfiles;
	unsigned int* sizes = NULL;
	if(my_rank == 0)	sizes = (unsigned int*) malloc(np*sizeof(unsigned int));
	unsigned int* my_steps_on_disk = (unsigned int*) malloc(2*my_save_size*sizeof(unsigned int));	//1st for id, 2nd for disk_iterations
	unsigned int** steps_on_disk = NULL;
	if(my_rank == 0)	steps_on_disk = (unsigned int**) malloc(np*sizeof(unsigned int*));

	//Close the temp file, if open
	if(my_tempfile && *my_tempfile)
	{
		fclose(*my_tempfile);
		*my_tempfile = NULL;
	}

	//Send the total number of links to proc 0
	MPI_Gather(&my_save_size,1,MPI_UNSIGNED,sizes,1,MPI_UNSIGNED,0,MPI_COMM_WORLD);
	if(my_rank == 0)
		for(i=1;i<np;i++)	steps_on_disk[i] = (unsigned int*) malloc(2*sizes[i]*sizeof(unsigned int));

	//Now send the total number of iterations on disk to proc 0
	counter = 0;
	for(j=0;j<save_size;j++)
	{
		loc = find_link_by_idtoloc(save_list[j],id_to_loc,N);
		if(assignments[loc] == my_rank)
		{
			my_steps_on_disk[counter] = save_list[j];
			my_steps_on_disk[counter+1] = sys[loc]->disk_iterations;
			counter += 2;
		}
	}

	if(my_rank == 0)	//!!!! Gather? !!!!
	{
		steps_on_disk[0] = my_steps_on_disk;
		for(i=1;i<np;i++)		MPI_Recv(steps_on_disk[i],2*sizes[i],MPI_UNSIGNED,i,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	}
	else
		MPI_Send(my_steps_on_disk,2*my_save_size,MPI_UNSIGNED,0,0,MPI_COMM_WORLD);

	//Send to proc 0 the name of all temp files
	char** tempfilenames = NULL;
	if(my_rank == 0)
	{
		int size;
		tempfilenames = (char**) malloc(np*sizeof(char*));
		for(i=0;i<np;i++)	tempfilenames[i] = (char*) malloc(GlobalVars->string_size*sizeof(char));
		strcpy(tempfilenames[0],GlobalVars->temp_filename);
		for(i=1;i<np;i++)
		{
			MPI_Recv(&size,1,MPI_INT,i,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			MPI_Recv(tempfilenames[i],size+1,MPI_CHAR,i,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		}
	}
	else
	{
		int size = strlen(GlobalVars->temp_filename);
		MPI_Send(&size,1,MPI_INT,0,0,MPI_COMM_WORLD);
		MPI_Send(GlobalVars->temp_filename,size+1,MPI_CHAR,0,0,MPI_COMM_WORLD);
	}

	//Find total size of a line in the temp files
	unsigned int line_size = CalcTotalOutputSize(GlobalVars);

	if(my_rank == 0 && GlobalVars->output_flag == 0)	//.dat
	{
		//Create output .dat file
		if(GlobalVars->print_par_flag == 1)
		{
			if(!additional_out)
				sprintf(outputfilename,"%s",GlobalVars->hydros_loc_filename);
			else
				sprintf(outputfilename,"%s_%s",GlobalVars->hydros_loc_filename,additional_out);
			for(i=0;i<GlobalVars->global_params->dim;i++)
			{
				sprintf(filenamespace,"_%.4e",GlobalVars->global_params->ve[i]);
				strcat(outputfilename,filenamespace);
			}
			sprintf(filenamespace,".dat");
			strcat(outputfilename,filenamespace);
		}
		else
			if(!additional_out)
				sprintf(outputfilename,"%s.dat",GlobalVars->hydros_loc_filename);
			else
				sprintf(outputfilename,"%s_%s.dat",GlobalVars->hydros_loc_filename,additional_out);
		outputfile = fopen(outputfilename,"w");
		if(!outputfile)
		{
			printf("[%i]: Error opening outputfile %s.\n",my_rank,outputfilename);
			return 2;
		}

		//Header
		fprintf(outputfile,"%i\n%i\n",save_size,dim);

		//Go through each temporary file
		for(i=0;i<save_size;i++)
		{
			loc = find_link_by_idtoloc(save_list[i],id_to_loc,N);
			proc = assignments[loc];

			//Open the temp file
			if(!additional_temp)
				sprintf(filename,"%s",tempfilenames[proc]);
				//sprintf(filename,"%s_%.3i",GlobalVars->temp_filename,proc);
			else
				sprintf(filename,"%s_%s",tempfilenames[proc],additional_temp);				//!!!! When is additional_temp useful? Get rid of it? Global params?? !!!!
				//sprintf(filename,"%s_%s_%.3i",GlobalVars->temp_filename,additional_temp,proc);
			inputfile = fopen(filename,"rb");
			if(!inputfile)
			{
				printf("\n[%i]: Error opening inputfile %s.\n",my_rank,filename);
				return 2;
			}

			//Find the link in the temp file inputfile
			counter = 0;	//index of link in my_save_list of its process
			fread(&id,sizeof(unsigned int),1,inputfile);
			fread(&total_spaces,sizeof(unsigned int),1,inputfile);
			while(id != save_list[i] && !feof(inputfile))
			{
				//fseek(inputfile,total_spaces * dim * sizeof(double),SEEK_CUR);
				fseek(inputfile,total_spaces * line_size,SEEK_CUR);
				counter++;
				fread(&id,sizeof(unsigned int),1,inputfile);
				fread(&total_spaces,sizeof(unsigned int),1,inputfile);
			}
			if(feof(inputfile))
			{
				printf("\n[%i]: Error: could not find id %u in temp file %s.\n",my_rank,save_list[i],filename);
				return 2;
			}

			//Write id and number of steps
			fprintf(outputfile,"\n%u %u\n",id,steps_on_disk[proc][2*counter+1]);

			//Read data in the temp file
			for(k=0;k<steps_on_disk[proc][2*counter+1];k++)
			{
				for(m=0;m<dim;m++)
				{
					fread(data_storage,GlobalVars->output_sizes[m],1,inputfile);
					WriteValue(outputfile,GlobalVars->output_specifiers[m],data_storage,GlobalVars->output_types[m]," ");
				}
				fprintf(outputfile,"\n");
			}

			//Skip over the last unused space
			for(k=0;k<total_spaces-steps_on_disk[proc][2*counter+1];k++)
				for(m=0;m<dim;m++)	fread(data_storage,GlobalVars->output_types[m],1,inputfile);

			fclose(inputfile);
		}

		//Clean up
		for(i=1;i<np;i++)	free(steps_on_disk[i]);
		free(steps_on_disk);
		fclose(outputfile);
		free(sizes);
		for(i=0;i<np;i++)	free(tempfilenames[i]);
		free(tempfilenames);
	}
	else if(my_rank == 0 && GlobalVars->output_flag == 1)	//.csv
	{
		if(GlobalVars->print_par_flag == 1)
		{
			if(!additional_out)
				sprintf(outputfilename,"%s",GlobalVars->hydros_loc_filename);
			else
				sprintf(outputfilename,"%s_%s",GlobalVars->hydros_loc_filename,additional_out);
			for(i=0;i<GlobalVars->global_params->dim;i++)
			{
				sprintf(filenamespace,"_%.4e",GlobalVars->global_params->ve[i]);
				strcat(outputfilename,filenamespace);
			}
			sprintf(filenamespace,".csv");
			strcat(outputfilename,filenamespace);
		}
		else
			if(!additional_out)
				sprintf(outputfilename,"%s.csv",GlobalVars->hydros_loc_filename);
			else
				sprintf(outputfilename,"%s_%s.csv",GlobalVars->hydros_loc_filename,additional_out);
		outputfile = fopen(outputfilename,"w");
		if(!outputfile)
		{
			printf("\n%i: Error opening outputfile %s.\n",my_rank,outputfilename);
			abort();
		}

		//Initializations
		space_counter = (unsigned int**) malloc(np*sizeof(unsigned int*));
		for(i=0;i<np;i++)	space_counter[i] = (unsigned int*) calloc(sizes[i],sizeof(unsigned int));
		max_disk = 0;
		for(i=0;i<np;i++)
		{
			for(j=0;j<sizes[i];j++)
				max_disk = (max_disk < steps_on_disk[i][2*j+1]) ? steps_on_disk[i][2*j+1] : max_disk;
		}
		positions = (fpos_t**) malloc(np*sizeof(fpos_t*));
		for(i=0;i<np;i++)	positions[i] = (fpos_t*) malloc(sizes[i]*sizeof(fpos_t));
		allinputfiles = (FILE**) malloc(np*sizeof(FILE*));
		for(i=0;i<np;i++)
		{
			if(sizes[i] > 0)
			{
				if(!additional_temp)
					sprintf(filename,"%s",tempfilenames[proc]);
					//sprintf(filename,"%s_%.3i",GlobalVars->temp_filename,i);
				else
					sprintf(filename,"%s_%s",tempfilenames[proc],additional_temp);
					//sprintf(filename,"%s_%s_%.3i",GlobalVars->temp_filename,additional_temp,i);
				allinputfiles[i] = fopen(filename,"rb");
				if(!allinputfiles[i])
				{
					printf("\n%i: Error opening inputfile %s.\n",my_rank,filename);
					return 2;
				}
			}
			else
				allinputfiles[i] = NULL;
		}

		//Find starting position in each file
		for(i=0;i<np;i++)
		{
			for(j=0;j<sizes[i];j++)
			{
				fseek(allinputfiles[i],sizeof(unsigned int),SEEK_CUR);	//Skip ID
				fread(&total_spaces,sizeof(unsigned int),1,allinputfiles[i]);	//Grab total spaces for this link
				fgetpos(allinputfiles[i],&(positions[i][j]));		//Here is the place to start for this link
				fseek(allinputfiles[i],total_spaces*line_size,SEEK_CUR);	//Skip to next link
			}
		}

		//Make the .csv header
		for(i=0;i<np;i++)
		{
			for(j=0;j<sizes[i];j++)
			{
				fprintf(outputfile,"Link %u,",steps_on_disk[i][2*j]);
				for(k=0;k<dim-1;k++)	fprintf(outputfile," ,");
			}
		}
		fprintf(outputfile,"\n");
		for(i=0;i<np;i++)
		{
			for(j=0;j<sizes[i];j++)
				for(k=0;k<dim;k++)	fprintf(outputfile,"Output_%u,",k);	//!!!! Use names !!!!
		}
		fprintf(outputfile,"\n");

		//Now make the .csv body
		for(m=0;m<max_disk;m++)
		{
			for(i=0;i<np;i++)
			{
				for(k=0;k<sizes[i];k++)
				{
					if(space_counter[i][k] > steps_on_disk[i][2*k+1])	//This link is done, leave blanks
						for(l=0;l<dim;l++)	fprintf(outputfile,",");
					else
					{
						fsetpos(allinputfiles[i],&(positions[i][k]));
						for(l=0;l<dim;l++)
						{
							fread(data_storage,GlobalVars->output_sizes[l],1,allinputfiles[i]);
							WriteValue(outputfile,GlobalVars->output_specifiers[l],data_storage,GlobalVars->output_types[l],",");
						}
						fgetpos(allinputfiles[i],&(positions[i][k]));
						(space_counter[i][k])++;
					}
				}
			}
			fprintf(outputfile,"\n");
		}

		//Clean up
		free(allinputfiles);
		for(i=1;i<np;i++)	free(steps_on_disk[i]);
		free(steps_on_disk);
		fclose(outputfile);
		free(sizes);
		for(i=0;i<np;i++)
		{
			free(positions[i]);
			free(space_counter[i]);
		}
		free(positions);
		free(space_counter);
		for(i=0;i<np;i++)	free(tempfilenames[i]);
		free(tempfilenames);
	}

	if(my_rank == 0)
		printf("\nResults written to file %s.\n",outputfilename);

	free(my_steps_on_disk);

	MPI_Barrier(MPI_COMM_WORLD);

	//Reopen the tempfile
	if(my_tempfile && my_save_size > 0)
	{
		if(additional_temp != NULL)
			sprintf(filename,"%s_%s",GlobalVars->temp_filename,additional_temp);
		else
			sprintf(filename,"%s",GlobalVars->temp_filename);
		*my_tempfile = fopen(filename,"r+b");
		if(!*my_tempfile)
			printf("[%i]: Error reopening temp file %s.\n",my_rank,filename);
	}

	return 0;
}

void PrepareDatabaseTable(UnivVars* GlobalVars,ConnData* conninfo)
{
	unsigned int num_cols,i;
	PGresult* res;
	char* query;

	if(conninfo == NULL)	return;

	if(my_rank == 0)
	{
		ConnectPGDB(conninfo);
		query = conninfo->query;

		//Create table
		if(conninfo->num_queries > 1)
		{
			sprintf(query,conninfo->queries[1],GlobalVars->hydro_table);
			res = PQexec(conninfo->conn,query);
			CheckResError(res,"deleting hydrograph table");
			PQclear(res);
		}
		sprintf(query,conninfo->queries[0],GlobalVars->hydro_table);
		res = PQexec(conninfo->conn,query);
		CheckResError(res,"creating hydrograph table");
		PQclear(res);

		//Make sure table is consistent with outputs
		sprintf(query,"SELECT data_type FROM information_schema.columns WHERE table_name='%s' ORDER BY ordinal_position;",GlobalVars->hydro_table);
		res = PQexec(conninfo->conn,query);
		CheckResError(res,"obtaining data types from hydrograph table");

		num_cols = PQntuples(res);
		if(num_cols > GlobalVars->num_print)
		{
			printf("[%i]: Error: need more outputs in .gbl file.\n",my_rank);
			MPI_Abort(MPI_COMM_WORLD,1);
		}
		else if(num_cols < GlobalVars->num_print)
		{
			printf("[%i]: Error: need more columns in table %s. Got %i, expected %i.\n",my_rank,GlobalVars->hydro_table,num_cols,GlobalVars->num_print);
			MPI_Abort(MPI_COMM_WORLD,1);
		}

		for(i=0;i<num_cols;i++)
		{
			if(GlobalVars->output_types[i] == ASYNCH_BAD_TYPE)
			{
				printf("[%i]: Error: output %i must be set to prepare database tables.\n",my_rank,i);
				MPI_Abort(MPI_COMM_WORLD,1);
			}
			else if( strcmp(PQgetvalue(res,i,0),"double precision") == 0 )
			{
				if(GlobalVars->output_types[i] != ASYNCH_DOUBLE)
				{
					printf("[%i]: Error: Output %i is of type double precision in output table, but should not be. %hi\n",my_rank,i,GlobalVars->output_types[i]);
					MPI_Abort(MPI_COMM_WORLD,1);
				}
			}
			else if( strcmp(PQgetvalue(res,i,0),"integer") == 0 )
			{
				if(GlobalVars->output_types[i] != ASYNCH_INT)
				{
					printf("[%i]: Error: Output %i is of type integer in output table, but should not be. %hi\n",my_rank,i,GlobalVars->output_types[i]);
					MPI_Abort(MPI_COMM_WORLD,1);
				}
			}
			else if( strcmp(PQgetvalue(res,i,0),"single precision") == 0 )
			{
				if(GlobalVars->output_types[i] != ASYNCH_FLOAT)
				{
					printf("[%i]: Error: Output %i is of type single precision in output table, but should not be. %hi\n",my_rank,i,GlobalVars->output_types[i]);
					MPI_Abort(MPI_COMM_WORLD,1);
				}
			}
			else if( strcmp(PQgetvalue(res,i,0),"short integer") == 0 )
			{
				if(GlobalVars->output_types[i] != ASYNCH_SHORT)
				{
					printf("[%i]: Error: Output %i is of type short integer in output table, but should not be. %hi\n",my_rank,i,GlobalVars->output_types[i]);
					MPI_Abort(MPI_COMM_WORLD,1);
				}
			}
			else if( strcmp(PQgetvalue(res,i,0),"character") == 0 )
			{
				if(GlobalVars->output_types[i] != ASYNCH_CHAR)
				{
					printf("[%i]: Error: Output %i is of type character in output table, but should not be. %hi\n",my_rank,i,GlobalVars->output_types[i]);
					MPI_Abort(MPI_COMM_WORLD,1);
				}
			}
			else
			{
				printf("[%i]: Error: Bad datatype for output %i while preparing database tables. (%s)\n",my_rank,i,PQgetvalue(res,i,0));
				MPI_Abort(MPI_COMM_WORLD,1);
			}
		}

		PQclear(res);
		DisconnectPGDB(conninfo);
	}
}


//Assumes temp files.
//Return value = 0 means everything is good.
//1 means a database related error.
//2 means a file system related error.
int UploadHydrosDB(Link** sys,UnivVars* GlobalVars,unsigned int N,unsigned int* save_list,unsigned int save_size,unsigned int my_save_size,unsigned int** id_to_loc,int* assignments,char* additional_temp,char* additional_out,ConnData* conninfo,FILE** my_tempfile)
{
	unsigned int blocks = 10;
	int i,proc,k,nbytes,result,return_val = 0;
	unsigned int j,l,m,loc,id,counter,total_spaces,max_disk,**space_counter;
	char filename[GlobalVars->string_size],filenamespace[GlobalVars->string_size],temptablename[GlobalVars->query_size];
	char* submission;
	//double *data_storage;
	char data_storage[16];
	fpos_t **positions;
	unsigned int dim = GlobalVars->num_print;
	FILE *inputfile,**allinputfiles;
	unsigned int* sizes = NULL;
	if(my_rank == 0)	sizes = (unsigned int*) malloc(np*sizeof(unsigned int));
	unsigned int* my_steps_on_disk = (unsigned int*) malloc(2*my_save_size*sizeof(unsigned int));	//1st for id, 2nd for disk_iterations
	unsigned int** steps_on_disk = NULL;
	PGresult *res;
	//unsigned int line_size = CalcTotalOutputSize(GlobalVars);

	//Close the temp file, if open
	if(my_tempfile && *my_tempfile)
	{
		fclose(*my_tempfile);
		*my_tempfile = NULL;
	}

	if(my_rank == 0)
	{
		steps_on_disk = (unsigned int**) malloc(np*sizeof(unsigned int*));
		//Assumes 24 chars max per double, plus delims, plus link id. So %.8e should work.
		submission = (char*) malloc( (dim*24 + dim) * sizeof(char));
	}

	//Send the total number of links to proc 0
	MPI_Gather(&my_save_size,1,MPI_UNSIGNED,sizes,1,MPI_UNSIGNED,0,MPI_COMM_WORLD);
	if(my_rank == 0)
		for(i=1;i<np;i++)	steps_on_disk[i] = (unsigned int*) malloc(2*sizes[i]*sizeof(unsigned int));

	//Find total size of a line in the temp files
	unsigned int line_size = CalcTotalOutputSize(GlobalVars);

	//Now send the total number of iterations on disk to proc 0
	counter = 0;
	for(j=0;j<save_size;j++)
	{
		loc = find_link_by_idtoloc(save_list[j],id_to_loc,N);
		if(assignments[loc] == my_rank)
		{
			my_steps_on_disk[counter] = save_list[j];
			my_steps_on_disk[counter+1] = sys[loc]->disk_iterations;
			counter += 2;
		}
	}

	if(my_rank == 0)	//!!!! Gather? !!!!
	{
		steps_on_disk[0] = my_steps_on_disk;
		for(i=1;i<np;i++)		MPI_Recv(steps_on_disk[i],2*sizes[i],MPI_UNSIGNED,i,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	}
	else
		MPI_Send(my_steps_on_disk,2*my_save_size,MPI_UNSIGNED,0,0,MPI_COMM_WORLD);

	//Send to proc 0 the name of all temp files
	char** tempfilenames = NULL;
	if(my_rank == 0)
	{
		int size = 0;
		tempfilenames = (char**) malloc(np*sizeof(char*));
		for(i=0;i<np;i++)	tempfilenames[i] = (char*) malloc(GlobalVars->string_size*sizeof(char));
		strcpy(tempfilenames[0],GlobalVars->temp_filename);
		for(i=1;i<np;i++)
		{
			MPI_Recv(&size,1,MPI_INT,i,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			MPI_Recv(tempfilenames[i],size+1,MPI_CHAR,i,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		}
	}
	else
	{
		int size = strlen(GlobalVars->temp_filename);
		MPI_Send(&size,1,MPI_INT,0,0,MPI_COMM_WORLD);
		MPI_Send(GlobalVars->temp_filename,size+1,MPI_CHAR,0,0,MPI_COMM_WORLD);
	}

	if(my_rank == 0)
	{
		sprintf(temptablename,"tmp_%s",GlobalVars->hydro_table);
		return_val = 0;
		ConnectPGDB(conninfo);

		//Delete temporary table, if it exists
		sprintf(conninfo->query,"DROP TABLE IF EXISTS %s;",temptablename);
		res = PQexec(conninfo->conn,conninfo->query);
		if(CheckResError(res,"dropping temporary output table"))
		{
			return_val = 1;
			goto finished;
		}
		PQclear(res);

		//Create temporary table
		sprintf(conninfo->query,conninfo->queries[0],temptablename);
		res = PQexec(conninfo->conn,conninfo->query);
		if(CheckResError(res,"creating temporary output table"))
		{
			return_val = 1;
			goto finished;
		}
		PQclear(res);

		//Tell database to prepare for copying
		sprintf(conninfo->query,"COPY %s FROM STDIN WITH DELIMITER ',';",temptablename);
		//sprintf(conninfo->query,"COPY %s FROM STDIN WITH DELIMITER ',';",GlobalVars->hydro_table);	//!!!! This should use a temp table !!!!
		res = PQexec(conninfo->conn,conninfo->query);
		if(CheckResState(res,PGRES_COPY_IN))
		{
			return_val = 1;
			goto finished;
		}
		PQclear(res);

		for(i=0;i<save_size;i++)	//!!!! Yeah, this can be done better. Go file by file. !!!!
		{
			loc = find_link_by_idtoloc(save_list[i],id_to_loc,N);
			proc = assignments[loc];

			//Open the temp file
			if(!additional_temp)
				sprintf(filename,"%s",tempfilenames[proc]);
				//sprintf(filename,"%s_%.3i",GlobalVars->temp_filename,proc);
			else
				sprintf(filename,"%s_%s",tempfilenames[proc],additional_temp);
				//sprintf(filename,"%s_%s_%.3i",GlobalVars->temp_filename,additional_temp,proc);
			inputfile = fopen(filename,"rb");
			if(!inputfile)
			{
				printf("\n[%i]: Error opening inputfile %s.\n",my_rank,filename);
				return_val = 2;
				goto finished;	//Don't bother uploading anymore
			}

			//Find the link in the temp file inputfile
			counter = 0;	//index of link in my_save_list of its process
			fread(&id,sizeof(unsigned int),1,inputfile);
			fread(&total_spaces,sizeof(unsigned int),1,inputfile);

			while(id != save_list[i] && !feof(inputfile))
			{
				//fseek(inputfile,total_spaces * dim * sizeof(double),SEEK_CUR);
				fseek(inputfile,total_spaces * line_size,SEEK_CUR);
				counter++;
				fread(&id,sizeof(unsigned int),1,inputfile);
				fread(&total_spaces,sizeof(unsigned int),1,inputfile);
			}
			if(feof(inputfile))
			{
				printf("\n[%i]: Error: could not find id %u in temp file %s.\n",my_rank,save_list[i],filename);
				return_val = 2;
				goto finished;	//Don't bother uploading anymore
			}

			//Now read in the data in the temp file, and submit them to the database
			for(k=0;k<steps_on_disk[proc][2*counter+1];k++)
			{
				nbytes = 0;
				for(j=0;j<dim-1;j++)
				{
					fread(data_storage,GlobalVars->output_sizes[j],1,inputfile);
					nbytes += CatBinaryToString(&(submission[nbytes]),GlobalVars->output_specifiers[j],data_storage,GlobalVars->output_types[j],",");
				}
				fread(data_storage,GlobalVars->output_sizes[j],1,inputfile);
				nbytes += CatBinaryToString(&(submission[nbytes]),GlobalVars->output_specifiers[j],data_storage,GlobalVars->output_types[j],"\n");
				result = PQputCopyData(conninfo->conn,submission,nbytes);
				if(result != 1)
				{
					printf("[%i]: Error: copy returned result %i.\n",my_rank,result);
					return_val = 1;
					goto finished;	//Don't bother uploading anymore
				}
			}

			//Skip over the last unused space	!!!! Is this really necesssary? !!!!
			for(k=0;k<total_spaces-steps_on_disk[proc][2*counter+1];k++)
				for(j=0;j<dim;j++)
					fread(data_storage,GlobalVars->output_types[j],1,inputfile);

			fclose(inputfile);
		}

		finished:

		//Tell database everything is uploaded
		i = PQputCopyEnd(conninfo->conn,NULL);
		if(i != 1)
		{
			printf("[%i]: Returned %i while closing copy to hydrograph database.\n",my_rank,i);
			return_val = 1;
		}

		//If the temporary table was loaded successfully, inserted it into the main table
		if(!return_val)
		{
			sprintf(conninfo->query,"INSERT INTO %s (SELECT * FROM %s);",GlobalVars->hydro_table,temptablename);
			res = PQexec(conninfo->conn,conninfo->query);
			CheckResError(res,"dropping temporary output table");
			PQclear(res);

			//Delete temporary table, if it exists
			sprintf(conninfo->query,"DROP TABLE IF EXISTS %s;",temptablename);
			res = PQexec(conninfo->conn,conninfo->query);
			CheckResError(res,"dropping temporary output table");
			PQclear(res);
		}

		//Clean up
		DisconnectPGDB(conninfo);
		for(i=1;i<np;i++)	free(steps_on_disk[i]);
		free(steps_on_disk);
		free(sizes);
		free(submission);
		for(i=0;i<np;i++)	free(tempfilenames[i]);
		free(tempfilenames);

		if(!return_val)	printf("\nResults written to table %s.\n",GlobalVars->hydro_table);
	}

	free(my_steps_on_disk);	
	MPI_Bcast(&return_val,1,MPI_INT,0,MPI_COMM_WORLD);

	//Reopen the tempfile
	if(my_save_size > 0)
	{
		if(additional_temp != NULL)
			sprintf(filename,"%s_%s",GlobalVars->temp_filename,additional_temp);
			//sprintf(filename,"%s_%s_%.3i",GlobalVars->temp_filename,additional_temp,my_rank);
		else
			sprintf(filename,"%s",GlobalVars->temp_filename);
			//sprintf(filename,"%s_%.3i",GlobalVars->temp_filename,my_rank);
		*my_tempfile = fopen(filename,"r+b");
		if(!*my_tempfile)
			printf("[%i]: Error reopening temp file %s.\n",my_rank,filename);
	}

	return return_val;
}


//Creates a peakflow file.
//Make sure GlobalVars has been initialized first.
int PreparePeakFlowFiles(UnivVars* GlobalVars,unsigned int peaksave_size)
{
	unsigned int i;

	if(!GlobalVars)
	{
		printf("[%i]: Error: A .gbl file must be read before preparing peakflow files.\n",my_rank);
		return 1;
	}

	if(peaksave_size)
	{
		//Put together the output filename string
		char outputfilename[GlobalVars->string_size],filename[GlobalVars->string_size];
		sprintf(outputfilename,"%s",GlobalVars->peaks_loc_filename);
		if(GlobalVars->print_par_flag == 1)
		{
			for(i=0;i<GlobalVars->global_params->dim;i++)
			{
				sprintf(filename,"_%.4e",GlobalVars->global_params->ve[i]);
				strcat(outputfilename,filename);
			}
		}

		//Setup peakflow file
		FILE* peakfile = NULL;
		GlobalVars->peakfilename = (char*) malloc(GlobalVars->string_size*sizeof(char));

		//Setup the output peakflows files
		sprintf(GlobalVars->peakfilename,"%s.pea",outputfilename);

		if(my_rank == 0)
		{
			peakfile = fopen(GlobalVars->peakfilename,"w");
			if(peakfile == NULL)
			{
				printf("[%i]: Error creating peak file %s.\n",my_rank,GlobalVars->peakfilename);
				return 1;
			}

			fprintf(peakfile,"%i\n%i\n\n",peaksave_size,GlobalVars->type);
			fclose(peakfile);
		}
	}

	return 0;
}


//Writes to disk the current peakflow information.
int DumpPeakFlowData(Link** sys,UnivVars* GlobalVars,unsigned int* my_sys,unsigned int my_N,unsigned int peaksave_size,ConnData* conninfo)
{
	unsigned int i,j;
	Link* current;
	FILE* peakfile;
	char buffer[256];

	if(peaksave_size)
	{
		if(!GlobalVars)
		{
			printf("[%i]: Error: A .gbl file must be read before creating peakflow files.\n",my_rank);
			return 1;
		}

		if(!(GlobalVars->peakfilename))
		{
			printf("[%i]: Error: Cannot write peakflow data to disk. Peakflow file not prepared.\n",my_rank);
			return 1;
		}

		//Write number of iterations to the peak data to peakfile
		double conversion = (GlobalVars->convertarea_flag) ? 1e-6 : 1.0;
		for(i=0;i<(unsigned int)np;i++)
		{
			if((unsigned int)my_rank == i)
			{
				peakfile = fopen(GlobalVars->peakfilename,"a");
				for(j=0;j<my_N;j++)
				{
					current = sys[my_sys[j]];
					if(current->peak_flag)
					{
						GlobalVars->peakflow_output(current->ID,current->peak_time,current->peak_value,current->params,GlobalVars->global_params,conversion,GlobalVars->area_idx,current->peakoutput_user,buffer);
						fprintf(peakfile,"%s",buffer);
					}
				}
				fclose(peakfile);
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
	}

	if(my_rank == 0)
		printf("Peakflows written to file %s.\n",GlobalVars->peakfilename);

	return 0;
}

//Uploads the current peakflow information to a database.
int UploadPeakFlowData(Link** sys,UnivVars* GlobalVars,unsigned int* my_sys,unsigned int my_N,unsigned int peaksave_size,ConnData* conninfo)
{
	unsigned int i,j,result,return_val = 0;
	char temptablename[GlobalVars->query_size];
	Link* current;
	PGresult *res;

if(conninfo->num_queries < 1)
printf("[%i]: Warning: I think you need a file to create a peakflow table...\n",my_rank);

	if(peaksave_size)
	{
		if(!GlobalVars)
		{
			printf("[%i]: Error: A .gbl file must be read before creating peakflow files.\n",my_rank);
			return 2;
		}

		//Prepare temporary table name
		sprintf(temptablename,"tmp_%s_%i",GlobalVars->peak_table,my_rank);

		if(my_rank == 0)
		{
			ConnectPGDB(conninfo);

			//Delete and create table. This should NOT be done in an init routine, in case there is an error.
			if(conninfo->num_queries > 1)
			{
				sprintf(conninfo->query,conninfo->queries[1],GlobalVars->peak_table);
				res = PQexec(conninfo->conn,conninfo->query);
				CheckResError(res,"deleting peakflow table");
				PQclear(res);		
			}
			//if(conninfo->num_queries > 0)
			{
				sprintf(conninfo->query,conninfo->queries[0],GlobalVars->peak_table);
				res = PQexec(conninfo->conn,conninfo->query);
				CheckResError(res,"creating peakflow table");
				PQclear(res);
			}

			DisconnectPGDB(conninfo);
		}

		MPI_Barrier(MPI_COMM_WORLD);

		//Write number of iterations to the peak data to peakfile	!!!! Were you smoking dope while writing this comment? What does it mean?? !!!!
		double conversion = (GlobalVars->convertarea_flag) ? 1e-6 : 1.0;
		for(i=0;i<(unsigned int)np+DB_CONNS_AT_ONCE;i+=DB_CONNS_AT_ONCE)
		{
			if(i <= (unsigned int)my_rank && (unsigned int)my_rank < i+DB_CONNS_AT_ONCE)
			{
				//Connect
				ConnectPGDB(conninfo);

				//Delete temporary table
				sprintf(conninfo->query,"DROP TABLE IF EXISTS %s;",temptablename);
				res = PQexec(conninfo->conn,conninfo->query);
				if(CheckResError(res,"dropping temporary peakflow table"))
					return_val = 1;
				PQclear(res);

				//Create temporary table
				sprintf(conninfo->query,conninfo->queries[0],temptablename);
				res = PQexec(conninfo->conn,conninfo->query);
				if(CheckResError(res,"creating temporary peakflow table"))
					return_val = 1;
				PQclear(res);

				//Tell database to prepare for copying
				sprintf(conninfo->query,"COPY %s FROM STDIN WITH DELIMITER ',';",temptablename);
				//sprintf(conninfo->query,"COPY %s FROM STDIN WITH DELIMITER ',';",GlobalVars->peak_table);
				res = PQexec(conninfo->conn,conninfo->query);
				if(CheckResState(res,PGRES_COPY_IN))
					return_val = 1;
				PQclear(res);

				//Transfer data
				if(!return_val)
				for(j=0;j<my_N;j++)
				{
					current = sys[my_sys[j]];
					if(current->peak_flag)
					{
						GlobalVars->peakflow_output(current->ID,current->peak_time,current->peak_value,current->params,GlobalVars->global_params,conversion,GlobalVars->area_idx,current->peakoutput_user,conninfo->query);
						result = PQputCopyData(conninfo->conn,conninfo->query,strlen(conninfo->query));
						if(result != 1)
						{
							printf("[%i]: Returned %i while copying to peakflow database.\n",my_rank,result);
							return_val = 1;
							break;
						}
					}
				}

				//Finish copy
				result = PQputCopyEnd(conninfo->conn,NULL);
				if(result != 1)
				{
					printf("[%i]: Returned %i while closing copy to peakflow database.\n",my_rank,result);
					return_val = 1;
				}

				//Insert temporary table, if no errors have occurred
				if(!return_val)
				{
					sprintf(conninfo->query,"INSERT INTO %s (SELECT * FROM %s);",GlobalVars->peak_table,temptablename);
					res = PQexec(conninfo->conn,conninfo->query);
					if(CheckResError(res,"inserting peakflow temporary table into permanent table"))
						return_val = 1;
					PQclear(res);
				}

				//Delete temporary table
				sprintf(conninfo->query,"DROP TABLE IF EXISTS %s;",temptablename);
				res = PQexec(conninfo->conn,conninfo->query);
				if(CheckResError(res,"dropping temporary peakflow table"))
					return_val = 1;
				PQclear(res);

				DisconnectPGDB(conninfo);
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
	}

	//Make sure everyone knows if an error occured
	MPI_Allreduce(&return_val,&i,1,MPI_INT,MPI_LOR,MPI_COMM_WORLD);

	if(my_rank == 0 && i == 0)
		printf("[%i]: Peakflows written to table %s.\n",my_rank,GlobalVars->peak_table);

	return i;
}

/*
//Writes the initial conditions of each link to disk as an .rec file.
//All data is sent to proc 0 so that the data file is stored correctly on disk.
void DataDump(Link** sys,unsigned int N,int* assignments,UnivVars* GlobalVars,unsigned int last_file)
{
	unsigned int i,j;
	FILE* output;
	char filename[256];
	unsigned int dim = GlobalVars->dim;
	double buffer[dim];

	if(my_rank == 0)	//Creating the file
	{
		sprintf(filename,"%s%u.rec",GlobalVars->dump_location,last_file);
		output = fopen(filename,"w");
		if(output == NULL)	printf("[%i]: Error opening file %s.\n",my_rank,filename);
		fprintf(output,"%hu\n%u\n0.0\n\n",GlobalVars->type,N);

		for(i=0;i<N;i++)
		{
			if(assignments[i] != 0)
				MPI_Recv(buffer,dim,MPI_DOUBLE,assignments[i],i,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			else
				for(j=0;j<dim;j++)	buffer[j] = sys[i]->list->tail->y_approx->ve[j];

			fprintf(output,"%u\n",sys[i]->ID);
			for(j=0;j<dim;j++)	fprintf(output,"%.6e ",buffer[j]);
			fprintf(output,"\n");
		}

		fclose(output);
	}
	else			//Sending data to proc 0
	{
		for(i=0;i<N;i++)
		{
			if(assignments[i] == my_rank)
			{
				for(j=0;j<dim;j++)	buffer[j] = sys[i]->list->tail->y_approx->ve[j];
				MPI_Send(buffer,dim,MPI_DOUBLE,0,i,MPI_COMM_WORLD);
			}
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);
}
*/


//This is the same as DataDump, but allows for a generic name
int DataDump2(Link** sys,unsigned int N,int* assignments,UnivVars* GlobalVars,char* name,ConnData* conninfo)
{
	unsigned int i,j;
	FILE* output;
	char filename[256];
	unsigned int dim = GlobalVars->dim;
	double buffer[dim];

	if(my_rank == 0)	//Creating the file
	{
		//sprintf(filename,"%s/%s.rec",GlobalVars->dump_location,name);
		if(name)	sprintf(filename,"%s%s.rec",GlobalVars->dump_loc_filename,name);
		else		sprintf(filename,"%s.rec",GlobalVars->dump_loc_filename);
		output = fopen(filename,"w");
		if(output == NULL)	printf("[%i]: Error opening file %s.\n",my_rank,filename);
		fprintf(output,"%hu\n%u\n0.0\n\n",GlobalVars->type,N);

		for(i=0;i<N;i++)
		{
			if(assignments[i] != 0)
				MPI_Recv(buffer,dim,MPI_DOUBLE,assignments[i],i,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			else
				for(j=0;j<dim;j++)	buffer[j] = sys[i]->list->tail->y_approx->ve[j];

			fprintf(output,"%u\n",sys[i]->ID);
			for(j=0;j<dim;j++)	fprintf(output,"%.6e ",buffer[j]);
			fprintf(output,"\n");
		}

		fclose(output);
	}
	else			//Sending data to proc 0
	{
		for(i=0;i<N;i++)
		{
			if(assignments[i] == my_rank)
			{
				for(j=0;j<dim;j++)	buffer[j] = sys[i]->list->tail->y_approx->ve[j];
				MPI_Send(buffer,dim,MPI_DOUBLE,0,i,MPI_COMM_WORLD);
			}
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);
	return 0;
}

//!!!! Upload to temp tables !!!!
int UploadDBDataDump(Link** sys,unsigned int N,int* assignments,UnivVars* GlobalVars,char* name,ConnData* conninfo)
{
	unsigned int i,j,init_length,nbytes;
	char query[GlobalVars->query_size],state_name[32],temptablename[GlobalVars->query_size];
	unsigned int dim = GlobalVars->dim,size = dim*(16+1)+16+1+16+1;
	char submission[size];	//Assumes 16 chars for each double
	static short int first_call = 0;
	PGresult* res;
	int result,error;
	for(i=0;i<size;i++)	submission[i] = 0;	//Just to keep valgrind from complaining

	//Since forecast_time does not change, just put them in submission once
	if(name)
	{
		strcpy(submission,name);
		init_length = strlen(name);
		submission[init_length++] = ',';
	}
	else
		init_length = 0;

	//Make a snapshot
	if(my_rank == 0)
	{
		ConnectPGDB(conninfo);
		error = 0;

		if(!first_call)	//!!!! Should this be done in an init routine? Does the file need an init routine? !!!!
		{
			//Trash table
			if(conninfo->num_queries > 1)
			{
				sprintf(query,conninfo->queries[1],GlobalVars->dump_table);
				res = PQexec(conninfo->conn,query);
				error |= CheckResError(res,"applying query 1 to snapshot table");
				PQclear(res);
			}

			//Create table
			if(conninfo->num_queries > 0)
			{
				sprintf(query,conninfo->queries[0],GlobalVars->dump_table);
				res = PQexec(conninfo->conn,query);
				error |= CheckResError(res,"applying query 0 to snapshot table");
				PQclear(res);
			}

			//Don't create the table again
			first_call = 1;
		}

		if(conninfo->num_queries > 2)	//For truncating the table
		{
			sprintf(query,conninfo->queries[2],GlobalVars->dump_table);
			res = PQexec(conninfo->conn,query);
			error |= CheckResError(res,"applying query 2 to snapshot table");
			PQclear(res);
		}

		//Create temp table
		sprintf(temptablename,"tmp_%s",GlobalVars->dump_table);

		sprintf(conninfo->query,"DROP TABLE IF EXISTS %s;",temptablename);
		res = PQexec(conninfo->conn,conninfo->query);
		CheckResError(res,"dropping temporary snapshot table");
		PQclear(res);

		sprintf(conninfo->query,conninfo->queries[0],temptablename);
		res = PQexec(conninfo->conn,conninfo->query);
		CheckResError(res,"creating temporary snapshot table");
		PQclear(res);

		//Tell database to prepare for copying
		if(!error)
		{
			sprintf(conninfo->query,"COPY %s FROM STDIN WITH DELIMITER ',' NULL AS 'NULL';",temptablename);
			//sprintf(conninfo->query,"COPY %s FROM STDIN WITH DELIMITER ',' NULL AS 'NULL';",GlobalVars->dump_table);
			res = PQexec(conninfo->conn,conninfo->query);
			error |= CheckResState(res,PGRES_COPY_IN);
			PQclear(res);
		}

		//Tell all procs if an error occured
		MPI_Bcast(&error,1,MPI_INT,0,MPI_COMM_WORLD);

		//Upload data
		if(!error)
		{
			for(i=0;i<N;i++)
			{
				nbytes = init_length;
				if(assignments[i] == my_rank)
				{
					nbytes += CatBinaryToString(&(submission[nbytes]),"%i",(char*)&(sys[i]->ID),ASYNCH_INT,",");
					for(j=0;j<dim;j++)
						nbytes += CatBinaryToString(&(submission[nbytes]),"%.6e",(char*)&(sys[i]->list->tail->y_approx->ve[j]),ASYNCH_DOUBLE,",");
					submission[nbytes-1] = '\n';
				}
				else
				{
					MPI_Recv(&(submission[init_length]),size-init_length,MPI_CHAR,assignments[i],sys[i]->ID,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
					nbytes += strlen(&(submission[init_length]));
				}

				result = PQputCopyData(conninfo->conn,submission,nbytes);
				if(result != 1)
				{
					printf("[%i]: Error: copy returned result %i.\n",my_rank,result);
					error = 1;
				}
			}
		}

		//Finish copy
		result = PQputCopyEnd(conninfo->conn,NULL);
		if(result != 1)
		{
			printf("[%i]: Returned %i while closing copy to peakflow database.\n",my_rank,result);
			error = 1;
		}

		//If no error, insert temporary table into permanent table
		if(!error)
		{
			sprintf(conninfo->query,"INSERT INTO %s (SELECT * FROM %s);",GlobalVars->dump_table,temptablename);
			res = PQexec(conninfo->conn,conninfo->query);
			if(CheckResError(res,"inserting temp snapshot table into permanent table"))
				error = 1;
			PQclear(res);
		}

		//Delete temporary table
		sprintf(conninfo->query,"DROP TABLE IF EXISTS %s;",temptablename);
		res = PQexec(conninfo->conn,conninfo->query);
		CheckResError(res,"dropping temporary snapshot table");
		PQclear(res);

		//Disconnect
		DisconnectPGDB(conninfo);
	}
	else			//Sending data to proc 0
	{
		//See if proc 0 encountered an error
		MPI_Bcast(&error,1,MPI_INT,0,MPI_COMM_WORLD);

		if(!error)
		{
			for(i=0;i<N;i++)
			{
				if(assignments[i] == my_rank)
				{
					nbytes = init_length;
					nbytes += CatBinaryToString(&(submission[nbytes]),"%i",(char*)&(sys[i]->ID),ASYNCH_INT,",");
					for(j=0;j<dim;j++)
						nbytes += CatBinaryToString(&(submission[nbytes]),"%.6e",(char*)&(sys[i]->list->tail->y_approx->ve[j]), ASYNCH_DOUBLE,",");
					submission[nbytes-1] = '\n';
					MPI_Send(&(submission[init_length]),size-init_length,MPI_CHAR,0,sys[i]->ID,MPI_COMM_WORLD);
				}
			}
		}
	}

	//See if proc 0 encountered an error
	MPI_Bcast(&error,1,MPI_INT,0,MPI_COMM_WORLD);

	return error;
}


FILE* PrepareTempFiles(Link** sys,unsigned int N,int* assignments,UnivVars* GlobalVars,unsigned int* save_list,unsigned int save_size,unsigned int my_save_size,char* additional,unsigned int** id_to_loc)
{
	unsigned int i,j,loc,start;
	Link* current;
	FILE* outputfile = NULL;
	//double* dummy_value;
	char filename[GlobalVars->string_size];
	VEC* dummy_y = v_get(GlobalVars->dim);	//Used for blanking lines in the temp files
	double dummy_t = 0.0;			//For blanking lines in the temp files
	//fpos_t holder1,holder2;
	long int current_pos = SEEK_SET,total;

	//Setup temporary output data file
	if(my_save_size > 0)
	{
		//dummy_value = (double*) malloc((1+GlobalVars->num_print)*sizeof(double));
		//for(i=0;i<1+GlobalVars->num_print;i++)	dummy_value[i] = 0.0;
		if(additional != NULL)
			sprintf(filename,"%s_%s",GlobalVars->temp_filename,additional);
			//sprintf(filename,"%s_%s_%.3i",GlobalVars->temp_filename,additional,my_rank);
		else
			sprintf(filename,"%s",GlobalVars->temp_filename);
			//sprintf(filename,"%s_%.3i",GlobalVars->temp_filename,my_rank);
		outputfile = fopen(filename,"w+b");
		if(!outputfile)	//If there's an error, try one more time...
		{		//This was added because of problems with the filesystem on Helium
			sleep(1);
			outputfile = fopen(filename,"w+b");
			if(outputfile)
				printf("[%i]: Notice: Needed two tries to create file %s.\n",my_rank,filename);
			else
			{
				printf("[%i]: Error: Could not create file %s.\n",my_rank,filename);
				MPI_Abort(MPI_COMM_WORLD,1);
			}
		}
		//if(GlobalVars->assim_flag == 1)	start = 0;
		//else				start = 1;
		start = 0;

		for(i=0;i<save_size;i++)
		{
			loc = find_link_by_idtoloc(save_list[i],id_to_loc,N);

			if(assignments[loc] == my_rank)
			{
				total = 0;
				current = sys[loc];
				fwrite(&(current->ID),sizeof(unsigned int),1,outputfile);

				//Calculate how many steps should be stored
				current->expected_file_vals = (unsigned int) rint(GlobalVars->maxtime / current->print_time) + 2;
				fwrite(&(current->expected_file_vals),sizeof(unsigned int),1,outputfile);
				//fgetpos(outputfile,&(current->pos));
				current_pos += 2 * sizeof(unsigned int);
				current->pos_offset = current_pos;

				//Fill out file
				for(j=start;j<current->expected_file_vals;j++)
					total += WriteStep(dummy_t,dummy_y,GlobalVars,current->params,current->iparams,current->state,outputfile,current->output_user,NULL);
					//WriteStep(dummy_t,dummy_y,GlobalVars,current->params,current->iparams,current->state,outputfile,current->output_user,&current_pos);

				//Update current_pos
				current_pos += total;
			}
		}

		//Add a few padding bytes to the end of the file.
		//This is to fix an issue with having the temp files open while by proc p while proc 0 reads them.
		//On Helium (at least), the last few bytes are not readable by 0 until p closes the file.
		//for(i=0;i<4;i++)	fwrite(&dummy_t,sizeof(double),1,outputfile);
	}

	v_free(dummy_y);
	MPI_Barrier(MPI_COMM_WORLD);
	return outputfile;
}

//Deletes this process's temporary file.
//Returns 0 if file delete, 1 if file does not exist, 2 if an error occurred.
int RemoveTemporaryFiles(UnivVars* GlobalVars,unsigned int my_save_size,char* additional_temp)
{
	int ret_val;
	char filename[GlobalVars->string_size];

	if(my_save_size > 0)
	{
		//Open the temp file
		if(!additional_temp)
			sprintf(filename,"%s",GlobalVars->temp_filename);
			//sprintf(filename,"%s_%.3i",GlobalVars->temp_filename,my_rank);
		else
			sprintf(filename,"%s_%s",GlobalVars->temp_filename,additional_temp);
			//sprintf(filename,"%s_%s_%.3i",GlobalVars->temp_filename,additional_temp,my_rank);
		ret_val = remove(filename);
		if(!ret_val)	return 0;
		else		return 2;
	}
	else
		return 1;
	
}

//Resets all temp files to the beginning of each link. This will allow all data in the temp files to be overwritten.
//Returns 2 if there is an error, 1 if a warning, 0 if good.
int ResetTempFiles(double set_time,Link** sys,unsigned int N,FILE* tempfile,UnivVars* GlobalVars,unsigned int my_save_size,unsigned int** id_to_loc)
{
	unsigned int i,j,id;
	long int current_pos = SEEK_SET;
	Link* current;

	//If nothing to write, don't reset!
	if(my_save_size == 0)	return 0;

	//Make sure file is open
	if(!tempfile)
	{
		printf("[%i]: Error: Cannot reset temporary file. A temporary file is not open.\n",my_rank);
		return 2;
	}

	//Restart the file
	rewind(tempfile);

	//Compute the number of bytes in for the component_idx
	unsigned int line_size = CalcTotalOutputSize(GlobalVars);

	//Find set_time at each link
	for(i=0;i<my_save_size;i++)
	{
		//Get linkid from file
		fread(&id,sizeof(unsigned int),1,tempfile);
		fseek(tempfile,sizeof(unsigned int),SEEK_CUR);	//Skip over number of stored values (I don't think it's needed here)
		current = sys[find_link_by_idtoloc(id,id_to_loc,N)];
		//current_pos = fgetpos(tempfile,&(current->pos));
		current_pos += 2 * sizeof(unsigned int);
		current->pos_offset = current_pos;

		//Set file position
		//fgetpos(tempfile,&(current->pos));
		current->disk_iterations = 0;
		current->next_save = set_time;		//!!!! This forces the print times to match up with the assimilation times !!!!

		//Get to next link in file
		fseek(tempfile,(current->expected_file_vals)*line_size,SEEK_CUR);
		current_pos += (current->expected_file_vals)*line_size;
	}

	return 0;
}


//Sets all temp files to set_value. The set_value is in the component_idx spot. All data that has been stored later can be overwritten.
//set_time is the time that is used for the next_save times.
//Returns 2 if there is an error, 1 if a warning, 0 if good.
int SetTempFiles(double set_time,void* set_value,short int data_type,unsigned int component_idx,Link** sys,unsigned int N,FILE* tempfile,UnivVars* GlobalVars,unsigned int my_save_size,unsigned int** id_to_loc,data_types* dt_info)
{
	unsigned int i,j,k,id;
	int warning = 0;
	long int current_pos = SEEK_SET;
	//double current_time;
	void* current_value;
	Link* current;

//printf("!!!! Note: SetTempFiles has not been tested. Be careful... !!!!\n");

	if(!tempfile)
	{
		printf("[%i]: Error: Cannot reset temporary file. A temporary file is not open.\n",my_rank);
		return 2;
	}

	if(component_idx >= GlobalVars->num_print)
	{
		printf("[%i]: Error: Cannot reset file to index %u. Only %u values are being outputed.\n",my_rank,component_idx,GlobalVars->num_print);
		return 2;
	}	

	//Restart the file
	rewind(tempfile);

	//Compute the number of bytes in for the component_idx
	unsigned int line_size = CalcTotalOutputSize(GlobalVars);
	unsigned int start_offset = 0;
	for(i=0;i<component_idx;i++)
		start_offset += GlobalVars->output_sizes[i];
	current_value = malloc(GlobalVars->output_sizes[component_idx]);
//unsigned int tempy,obtained;
//int got;
//double tempy_d;
	//Find set_time at each link
	for(i=0;i<my_save_size;i++)
	{
		//Get linkid from file
		fread(&id,sizeof(unsigned int),1,tempfile);
/*
printf("read id %u (%u)\n",id,obtained);
if(obtained == 0)
{
printf("Error function: %i %i\n",feof(tempfile),ferror(tempfile));
}
*/
		fseek(tempfile,sizeof(unsigned int),SEEK_CUR);	//Skip over number of stored values (I don't think it's needed here)
/*
obtained = fread(&tempy,sizeof(unsigned int),1,tempfile);
printf("read steps %u (%u)\n",tempy,obtained);
*/
		current = sys[find_link_by_idtoloc(id,id_to_loc,N)];
		current_pos += 2 * sizeof(unsigned int);
//printf("Got id = %u, stored %u, disk_iters = %u, start offset = %u\n",id,tempy,current->disk_iterations,start_offset);
/*
{
obtained = fread(&tempy_d,sizeof(double),1,tempfile);
printf("read %e (%u)\n",tempy_d,obtained);
obtained = fread(&tempy_d,sizeof(double),1,tempfile);
printf("read %e (%u)\n",tempy_d,obtained);
obtained = fread(&tempy_d,sizeof(double),1,tempfile);
printf("read %e (%u)\n",tempy_d,obtained);
obtained = fread(&tempy_d,sizeof(double),1,tempfile);
printf("read %e (%u)\n",tempy_d,obtained);
obtained = fread(&tempy_d,sizeof(double),1,tempfile);
printf("read %e (%u)\n",tempy_d,obtained);
obtained = fread(&tempy_d,sizeof(double),1,tempfile);
printf("read %e (%u)\n",tempy_d,obtained);
fseek(tempfile,-6*sizeof(double),SEEK_CUR);
}
*/
		//Move to component index
		fseek(tempfile,start_offset,SEEK_CUR);

		for(j=0;j<current->disk_iterations;j++)
		{
			fread(current_value,GlobalVars->output_sizes[component_idx],1,tempfile);
/*
if(id == 0)
{
printf("Comparing %e to %e. (%u) Got %i.\n",*(double*)set_value,*(double*)current_value,obtained,dt_info->is_equal[data_type](set_value,current_value));
printf("i = %i j = %i iters = %u\n",i,j,current->disk_iterations);
//getchar();
}
*/
			if(dt_info->is_equal[data_type](set_value,current_value))
			{
/*
if(id == 0)
{
printf("breaking...\n");
printf("i = %i j = %i iters = %u\n",i,j,current->disk_iterations);
//getchar();
}
*/
				//fseek(tempfile,-start_offset,SEEK_CUR);	//Backup to start of step
				fseek(tempfile,(long int) -(int)(start_offset+GlobalVars->output_sizes[component_idx]),SEEK_CUR);	//Backup to start of step
/*
if(id == 0)
{
long int val = (long int) -(int)(start_offset+GlobalVars->output_sizes[component_idx]);
printf("Went back %li, got = %i, ferror = %i, feof = %i\n",val,got,ferror(tempfile),feof(tempfile));
obtained = fread(&tempy_d,sizeof(double),1,tempfile);
printf("checky: read %e (%u %i)\n",tempy_d,obtained,feof(tempfile));
obtained = fread(&tempy_d,sizeof(double),1,tempfile);
printf("checky: read %e (%u %i)\n",tempy_d,obtained,feof(tempfile));
fseek(tempfile,-2*sizeof(double),SEEK_CUR);
}
*/
				break;
			}

			fseek(tempfile,line_size-GlobalVars->output_sizes[component_idx],SEEK_CUR);	//Skip to next step
			current_pos += line_size;
		}

		if(j == current->disk_iterations)
		{
			printf("[%i]: Warning: Cannot set file position for link %u.\n",my_rank,id);
			warning = 1;
			//continue;
/*
if(id == 0)
{
printf("Failed on %e\n",*(double*)set_value);
printf("i = %i j = %i iters = %u\n",i,j,current->disk_iterations);
//getchar();
}
*/
		}

		//Set file position
		//fgetpos(tempfile,&(current->pos));
		current->pos_offset = current_pos;
		current->disk_iterations = j;
		current->next_save = set_time;		//!!!! This forces the print times to match up with the assimilation times !!!!

/*
if(id == 0)
{
printf("i = %i j = %i line_size = %u expected = %u shift = %i pos_offset = %i\n",i,j,line_size,current->expected_file_vals,(current->expected_file_vals-j)*line_size,current_pos);
//getchar();
}
*/
		//Get to next link in file
		fseek(tempfile,(current->expected_file_vals-j)*line_size,SEEK_CUR);
		current_pos += (current->expected_file_vals-j)*line_size;
//printf("Seek gave %i, current_pos = %i\n",got,current_pos);
	}

	free(current_value);

	return warning;
}

//!!!! This assumes the first value in the file is time. !!!!
/*
int SetTempFiles(double set_time,Link** sys,unsigned int N,FILE* tempfile,UnivVars* GlobalVars,unsigned int my_save_size,unsigned int** id_to_loc)
{
	unsigned int i,j,id;
	int current_pos,warning = 0;
	double current_time;
	Link* current;

	if(!tempfile)
	{
		printf("[%i]: Error: Cannot reset temporary file. A temporary file is not open.\n",my_rank);
		return 2;
	}

	//Restart the file
	rewind(tempfile);

	//Find set_time at each link
	for(i=0;i<my_save_size;i++)
	{
		//Get linkid from file
		fread(&id,sizeof(unsigned int),1,tempfile);
		fseek(tempfile,sizeof(unsigned int),SEEK_CUR);	//Skip over number of stored values (I don't think it's needed here)
		current = sys[find_link_by_idtoloc(id,id_to_loc,N)];
		current_pos = fgetpos(tempfile,&(current->pos));
		for(j=0;j<current->disk_iterations;j++)
		{
			fread(&current_time,sizeof(double),1,tempfile);
			if( fabs(set_time - current_time) < 1e-10 )		//!!!! Is this a good bound? Maybe use a relative error. !!!!
			{
				fseek(tempfile,-sizeof(double),SEEK_CUR);	//Backup one double
				break;
			}

			//fseek(tempfile,GlobalVars->num_print*sizeof(double),SEEK_CUR);	//Skip ahead to next timestep
			fseek(tempfile,(GlobalVars->num_print-1)*sizeof(double),SEEK_CUR);	//Skip ahead to next timestep
		}

		if(j == current->disk_iterations)
		{
			printf("[%i]: Warning: Cannot set file position for link %u to time %f.\n",my_rank,id,set_time);
			warning = 1;
			continue;
		}

		//Set file position
		fgetpos(tempfile,&(current->pos));
		current->disk_iterations = j;
		current->next_save = set_time;		//!!!! This forces the print times to match up with the assimilation times !!!!

		//Get to next link in file
		//fseek(tempfile,(current->expected_file_vals-j)*(1+GlobalVars->num_print)*sizeof(double),SEEK_CUR);
		fseek(tempfile,(current->expected_file_vals-j)*GlobalVars->num_print*sizeof(double),SEEK_CUR);
	}

	return warning;
}
*/

//Read in a .rec file from disk and loads it into the intial condition for sys (tail).
//This does NOT set the current time at each link.
void LoadRecoveryFile(char* filename,Link** sys,unsigned int N,unsigned int my_N,unsigned int* assignments,UnivVars* GlobalVars)
{
	FILE* input;
	unsigned int i,j,read_type,read_N,id,counter=0;
	VEC* buffer = v_get(GlobalVars->problem_dim);
	MPI_Status status;

	MPI_Barrier(MPI_COMM_WORLD);

	if(my_rank == 0)
	{
		input = fopen(filename,"r");
		if(!input)
		{
			printf("[%i]: Error opening recovery file %s.\n",my_rank,filename);
			MPI_Abort(MPI_COMM_WORLD,1);
		}

		fscanf(input,"%u %u %*f",&read_type,&read_N);
		if(N != read_N)		//!!!! Ignoring read_type for now (190 vs 19) !!!!
		{
			printf("[%i]: Error reading recovery file: bad model type (%i) or wrong number of links (%i).\n",my_rank,read_type,read_N);
			MPI_Abort(MPI_COMM_WORLD,1);
		}

		for(i=0;i<N;i++)
		{
			fscanf(input,"%u",&id);
			if(sys[i]->ID != id)
			{
				printf("[%i]: Error reading recovery file: bad link id (%i); expected %i.\n",my_rank,id,sys[i]->ID);
				MPI_Abort(MPI_COMM_WORLD,1);
			}

			for(j=0;j<GlobalVars->problem_dim;j++)	fscanf(input,"%lf",&(buffer->ve[j]));

			if(assignments[i] == my_rank)	v_copy(buffer,sys[i]->list->tail->y_approx);
			else	MPI_Send(buffer->ve,buffer->dim,MPI_DOUBLE,assignments[i],i,MPI_COMM_WORLD);
		}
	}
	else
	{
		while(counter < my_N)
		{
			MPI_Recv(buffer->ve,buffer->dim,MPI_DOUBLE,0,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
			v_copy(buffer,sys[status.MPI_TAG]->list->tail->y_approx);
			counter++;
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);
}

//Rewrites the previous step written at link_i with the current step. The current time and state is used.
//!!!! This assumes we are writting into files. !!!!
//Returns 0 if all is well
//Returns 1 if there is no previous iteration to overwrite
//Returns 2 a step as been previously written, but it is not the expected number of bytes
int overwrite_last_step(Link* link_i,UnivVars *GlobalVars,FILE* outputfile)
{
	unsigned int i,step_byte_size = 0;

	//Check that something has actually been written for this link
	if(link_i->disk_iterations == 0)	return 1;

	//Backup a step in the file
	for(i=0;i<GlobalVars->num_print;i++)
		step_byte_size += GlobalVars->output_sizes[i];
	if(link_i->pos_offset < step_byte_size)	return 2;
	link_i->pos_offset -= step_byte_size;

	//Write the current step
	WriteStep(link_i->last_t,link_i->list->tail->y_approx,GlobalVars,link_i->params,link_i->iparams,link_i->state,outputfile,link_i->output_user,&(link_i->pos_offset));
	return 0;
}

/*
int ConvertBinaryToString(double* data_storage,char* submission,unsigned int blocks,unsigned int dimp1,unsigned int id)
{
	unsigned int i,j,written,count = 0;

	for(i=0;i<blocks;i++)
	{
		written = sprintf(&(submission[count]),"%u,",id);
		count += written;
		for(j=0;j<dimp1;j++)
		{
			written = sprintf(&(submission[count]),"%.8e,",data_storage[i*dimp1 + j]);
			count += written;
		}
		submission[count-1] = '\n';
	}

	return count;
}
*/

