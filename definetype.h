#ifndef DEFINETYPE_H
#define DEFINETYPE_H

#include <stdio.h>
#include <stdlib.h>
#include "structs.h"
#include "system.h"
#include "rkmethods.h"
#include "problems.h"
#include <unistd.h>
#include "math.h"

void SetParamSizes(UnivVars* GlobalVars);
void ConvertParams(VEC* params,unsigned int type);
void InitRoutines(Link* link,unsigned int type,unsigned int exp_imp,unsigned short int dam);
void Precalculations(Link* link_i,VEC* global_params,VEC* params,IVEC* iparams,unsigned int disk_params,unsigned int params_size,unsigned short int dam,unsigned int type);
unsigned int ReadInitData(VEC* global_params,VEC* params,IVEC* iparams,QVSData* qvs,unsigned short int dam,VEC* y_0,unsigned int type);
//void SetRKInitDim(unsigned int rktype,unsigned int inittype,unsigned int* rkdim,unsigned int* initdim);
void AssimError(unsigned int N,UnivVars* GlobalVars,ErrorData* GlobalErrors);
//void CheckConsistency(VEC* y,unsigned int dim,unsigned int type,VEC* params,VEC* global_params);

#endif
