#ifndef SOLVERS_H
#define SOLVERS_H

#include "structs.h"
#include "rainfall.h"
#include "comm.h"
#include "rkmethods.h"
#include "system.h"
#include "forcings.h"
#include <stdio.h>

extern int np;
extern int my_rank;

void AsynchSolver(Link** sys,unsigned int N,unsigned int* my_sys,unsigned int my_N,unsigned int my_max_nodes,UnivVars* GlobalVars,int* assignments,unsigned int** id_to_loc,TempStorage* workspace,Forcing** forcings,ConnData** db_connections,TransData* my_data,short int print_flag,FILE* outputfile);
void AsynchSolverPersis(Link** sys,unsigned int* my_sys,unsigned int my_N,unsigned int my_max_nodes,UnivVars* GlobalVars,int* assignments,TempStorage* workspace,Forcing** forcings,ConnData** db_connections,TransData* my_data,short int print_flag,FILE* outputfile);

#endif

