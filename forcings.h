#ifndef FORCINGS_H
#define FORCINGS_H

#include <stdio.h>
#include <stdlib.h>
#include "structs.h"
#include "rainfall.h"
#include "io.h"

extern int np;
extern int my_rank;

Forcing* InitializeForcings();
void FreeForcing(Forcing** forcings);

unsigned int PassesOther(Forcing* forcing);
unsigned int PassesBinaryFiles(Forcing* forcing);
unsigned int PassesDatabase(Forcing* forcing);
unsigned int PassesRecurring(Forcing* forcing);

double NextForcingOther(Link** sys,unsigned int N,unsigned int* my_sys,unsigned int my_N,int* assignments,UnivVars* GlobalVars,Forcing* forcing,ConnData** db_connections,unsigned int** id_to_loc,unsigned int forcing_idx);
double NextForcingBinaryFiles(Link** sys,unsigned int N,unsigned int* my_sys,unsigned int my_N,int* assignments,UnivVars* GlobalVars,Forcing* forcing,ConnData** db_connections,unsigned int** id_to_loc,unsigned int forcing_idx);
double NextForcingGZBinaryFiles(Link** sys,unsigned int N,unsigned int* my_sys,unsigned int my_N,int* assignments,UnivVars* GlobalVars,Forcing* forcing,ConnData** db_connections,unsigned int** id_to_loc,unsigned int forcing_idx);
double NextForcingDatabase(Link** sys,unsigned int N,unsigned int* my_sys,unsigned int my_N,int* assignments,UnivVars* GlobalVars,Forcing* forcing,ConnData** db_connections,unsigned int** id_to_loc,unsigned int forcing_idx);
double NextForcingRecurring(Link** sys,unsigned int N,unsigned int* my_sys,unsigned int my_N,int* assignments,UnivVars* GlobalVars,Forcing* forcing,ConnData** db_connections,unsigned int** id_to_loc,unsigned int forcing_idx);

#endif

