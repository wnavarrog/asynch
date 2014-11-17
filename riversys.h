#ifndef RIVERSYS_H
#define RIVERSYS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <libpq-fe.h>
#include "structs.h"
#include "rkmethods.h"
#include "partition.h"
#include "comm.h"
#include "sort.h"
#include "definetype.h"
#include "misc.h"
#include "io.h"
#include "forcings.h"
#include "modeloutputs.h"

extern int np;
extern int my_rank;

Link** Create_River_System_parallel(char rk_filename[],unsigned int* N,unsigned int** my_sys,unsigned int* my_N,unsigned int* my_max_nodes,TransData** my_data,int** assignments,short int** getting,RKMethod*** AllMethods,unsigned int* nummethods,UnivVars* GlobalVars,ErrorData* GlobalErrors,unsigned int** save_list,unsigned int* my_save_size,unsigned int* save_size,unsigned int* peaksave_size,unsigned int*** id_to_loc,Forcing** forcings,ConnData** db_connections,TempStorage** workspace,model* custom_model);

UnivVars* Read_Global_Data(char globalfilename[],ErrorData** GlobalErrors,Forcing** forcings,ConnData** db_connections,char* rkdfilename,model* custom_model);

unsigned int* Create_SAV_Data(char filename[],Link** sys,unsigned int N,unsigned int* size,ConnData *conninfo,unsigned short int flag);

void ReadLineGlobal(FILE* globalfile,char* linebuffer,unsigned int size,unsigned int string_size);

int ReadLineError(int valsread,int valswant,char message[]);

int RemoveSuffix(char* filename,char suffix[]);

int AttachParameters(char* filename,unsigned int max_size,VEC* v,unsigned int string_size);

#endif

