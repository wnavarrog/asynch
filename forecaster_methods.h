#ifndef FORECASTER_METHODS_H
#define FORECASTER_METHODS_H

#include "structs.h"
#include "comm.h"
#include <time.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

typedef struct ForecastData
{
	char* model_name;
	char* halt_filename;
	unsigned int forecasting_forcing;
	unsigned int num_rainsteps;
	short int ifis_display;
} ForecastData;

int DeleteFutureValues(ConnData* conninfo,unsigned int num_tables,UnivVars* GlobalVars,char* table_name,char* model_name,unsigned int clear_after,unsigned int equal);
void PerformTableMaintainance(ConnData* conninfo_hydros,UnivVars* GlobalVars,ForecastData* Forecaster,short int* vac,short unsigned int hr1,unsigned int num_tables,char* tablename);
void CheckPartitionedTable(ConnData* conninfo,UnivVars* GlobalVars,ForecastData* Forecaster,unsigned int num_tables,char* tablename,char* colname);
void CreateHaltFile(char* filename);
short int CheckFinished(char* filename);
ForecastData* Init_ForecastData(char* fcst_filename,unsigned int string_size);
void Free_ForecastData(ForecastData** Forecaster);

#endif

