#include "definetype.h"

//Sets the various sizes and flags for the model. This method should set the following fields:
//dim:			The number of unknowns in the differential equations (or the number of ODEs at each link).
//diff_start:		The starting index of the differential unknowns.
//str_flag:		1 if reading rainfall data from a .str file, 0 else.
//binrain_flag:		1 if reading rainfall data from binary files, 0 else.
//uses_dam:		1 if dams are compatible with the model given by type.
//params_size:		The total number of parameters (including precalculations) at links with no dam.
//dam_params_size:	The total number of parameters (including precalculations) at links with a dam.
//area_idx:		The entry in params where the upstream area is stored.
//disk_params:		The number of enries in param that are read from DEM data.
//Currently, this program assumes the same number of differential equations at each link.
//UnivVars* GlobalVars:	Contains the global variables for the system.
void SetParamSizes(UnivVars* GlobalVars)
{
	unsigned short int type = GlobalVars->type;
	unsigned int num_global_params;

	//Set dim and start of differential variables
	switch(type)
	{
//--------------------------------------------------------------------------------------------
		case 0:	GlobalVars->dim = GlobalVars->problem_dim = 1;
			GlobalVars->template_flag = 0;
			GlobalVars->assim_flag = 0;
			GlobalVars->diff_start = 0;
			GlobalVars->no_ini_start = GlobalVars->dim;
			num_global_params = 6;
			GlobalVars->uses_dam = 0;
			GlobalVars->params_size = 20;
			GlobalVars->iparams_size = 0;
			GlobalVars->dam_params_size = 0;
			GlobalVars->area_idx = 2;
			GlobalVars->areah_idx = 1;
			GlobalVars->disk_params = 12;
			GlobalVars->num_dense = 1;
			GlobalVars->convertarea_flag = 1;
			GlobalVars->num_forcings = 0;
			break;
//--------------------------------------------------------------------------------------------
		case 1: GlobalVars->dim = GlobalVars->problem_dim = 2;
			GlobalVars->template_flag = 0;
			GlobalVars->assim_flag = 0;
			GlobalVars->diff_start = 0;
			GlobalVars->no_ini_start = GlobalVars->dim;
			num_global_params = 6;
			GlobalVars->uses_dam = 0;
			GlobalVars->params_size = 20;
			GlobalVars->iparams_size = 0;
			GlobalVars->dam_params_size = 0;
			GlobalVars->area_idx = 2;
			GlobalVars->areah_idx = 1;
			GlobalVars->disk_params = 12;
			GlobalVars->num_dense = 1;
			GlobalVars->convertarea_flag = 1;
			GlobalVars->num_forcings = 1;
			break;
//--------------------------------------------------------------------------------------------
		case 2: GlobalVars->dim = GlobalVars->problem_dim = 2;
			GlobalVars->template_flag = 0;
			GlobalVars->assim_flag = 0;
			GlobalVars->diff_start = 0;
			GlobalVars->no_ini_start = GlobalVars->dim;
			num_global_params = 6;
			GlobalVars->uses_dam = 0;
			GlobalVars->params_size = 20;
			GlobalVars->iparams_size = 0;
			GlobalVars->dam_params_size = 0;
			GlobalVars->area_idx = 2;
			GlobalVars->areah_idx = 1;
			GlobalVars->disk_params = 12;
			GlobalVars->num_dense = 1;
			GlobalVars->convertarea_flag = 1;
			GlobalVars->num_forcings = 1;
			break;
//--------------------------------------------------------------------------------------------
		case 3: GlobalVars->dim = GlobalVars->problem_dim = 2;
			GlobalVars->template_flag = 0;
			GlobalVars->assim_flag = 0;
			GlobalVars->diff_start = 0;
			GlobalVars->no_ini_start = GlobalVars->dim;
			num_global_params = 6;
			GlobalVars->uses_dam = 0;
			GlobalVars->params_size = 20;
			GlobalVars->iparams_size = 0;
			GlobalVars->dam_params_size = 0;
			GlobalVars->area_idx = 2;
			GlobalVars->areah_idx = 1;
			GlobalVars->disk_params = 12;
			GlobalVars->num_dense = 1;
			GlobalVars->convertarea_flag = 1;
			GlobalVars->num_forcings = 1;
			break;
//--------------------------------------------------------------------------------------------
		case 4: GlobalVars->dim = GlobalVars->problem_dim = 4;
			GlobalVars->template_flag = 0;
			GlobalVars->assim_flag = 0;
			GlobalVars->diff_start = 0;
			GlobalVars->no_ini_start = GlobalVars->dim;
			num_global_params = 6;
			GlobalVars->uses_dam = 0;
			GlobalVars->params_size = 20;
			GlobalVars->iparams_size = 0;
			GlobalVars->dam_params_size = 0;
			GlobalVars->area_idx = 2;
			GlobalVars->areah_idx = 1;
			GlobalVars->disk_params = 12;
			GlobalVars->num_dense = 1;
			GlobalVars->convertarea_flag = 1;
			GlobalVars->num_forcings = 1;
			break;
//--------------------------------------------------------------------------------------------
		case 5: GlobalVars->dim = GlobalVars->problem_dim = 4;
			GlobalVars->template_flag = 0;
			GlobalVars->assim_flag = 0;
			GlobalVars->diff_start = 0;
			GlobalVars->no_ini_start = GlobalVars->dim;
			num_global_params = 6;
			GlobalVars->uses_dam = 0;
			GlobalVars->params_size = 20;
			GlobalVars->iparams_size = 0;
			GlobalVars->dam_params_size = 0;
			GlobalVars->area_idx = 2;
			GlobalVars->areah_idx = 1;
			GlobalVars->disk_params = 12;
			GlobalVars->num_dense = 1;
			GlobalVars->convertarea_flag = 1;
			GlobalVars->num_forcings = 1;
			break;
//--------------------------------------------------------------------------------------------
		case 6: GlobalVars->dim = GlobalVars->problem_dim = 4;
			GlobalVars->template_flag = 0;
			GlobalVars->assim_flag = 0;
			GlobalVars->diff_start = 0;
			GlobalVars->no_ini_start = GlobalVars->dim;
			num_global_params = 5;
			GlobalVars->uses_dam = 0;
			GlobalVars->params_size = 20;
			GlobalVars->iparams_size = 0;
			GlobalVars->dam_params_size = 0;
			GlobalVars->area_idx = 2;
			GlobalVars->areah_idx = 1;
			GlobalVars->disk_params = 14;
			GlobalVars->num_dense = 1;
			GlobalVars->convertarea_flag = 1;
			GlobalVars->num_forcings = 1;
			break;
//--------------------------------------------------------------------------------------------
		case 15: GlobalVars->dim = GlobalVars->problem_dim = 2;
			GlobalVars->template_flag = 0;
			GlobalVars->assim_flag = 0;
			GlobalVars->diff_start = 0;
			GlobalVars->no_ini_start = GlobalVars->dim;
			num_global_params = 6;
			GlobalVars->uses_dam = 0;
			GlobalVars->params_size = 20;
			GlobalVars->iparams_size = 0;
			GlobalVars->dam_params_size = 0;
			GlobalVars->area_idx = 2;
			GlobalVars->areah_idx = 1;
			GlobalVars->disk_params = 12;
			GlobalVars->num_dense = 1;
			GlobalVars->convertarea_flag = 0;
			GlobalVars->num_forcings = 1;
			break;
//--------------------------------------------------------------------------------------------
		case 19:	GlobalVars->dim = GlobalVars->problem_dim = 3;
				GlobalVars->template_flag = 0;
				GlobalVars->assim_flag = 0;
				GlobalVars->diff_start = 0;
				GlobalVars->no_ini_start = GlobalVars->dim;
				num_global_params = 7;
				GlobalVars->uses_dam = 0;
				GlobalVars->params_size = 8;
				GlobalVars->iparams_size = 0;
				GlobalVars->dam_params_size = 0;
				GlobalVars->area_idx = 0;
				GlobalVars->areah_idx = 2;
				GlobalVars->disk_params = 3;
				GlobalVars->num_dense = 1;
				GlobalVars->convertarea_flag = 0;
				GlobalVars->num_forcings = 1;
				break;
//--------------------------------------------------------------------------------------------
		case 20:	GlobalVars->dim = GlobalVars->problem_dim = 3;
				GlobalVars->template_flag = 0;
				GlobalVars->assim_flag = 0;
				GlobalVars->diff_start = 0;
				GlobalVars->no_ini_start = GlobalVars->dim;
				num_global_params = 9;
				GlobalVars->uses_dam = 0;
				GlobalVars->params_size = 6;
				GlobalVars->iparams_size = 0;
				GlobalVars->dam_params_size = 0;
				GlobalVars->area_idx = 0;
				GlobalVars->areah_idx = 2;
				GlobalVars->disk_params = 3;
				GlobalVars->num_dense = 1;
				GlobalVars->convertarea_flag = 0;
				GlobalVars->num_forcings = 1;
				break;
//--------------------------------------------------------------------------------------------
		case 21:	GlobalVars->dim = GlobalVars->problem_dim = 4;
				GlobalVars->template_flag = 0;
				GlobalVars->assim_flag = 0;
				GlobalVars->diff_start = 1;
				GlobalVars->no_ini_start = 2;
				num_global_params = 7;
				GlobalVars->uses_dam = 1;
				GlobalVars->params_size = 7;	//Need 1 extra for orifice_area
				GlobalVars->iparams_size = 0;
				GlobalVars->dam_params_size = 15;
				GlobalVars->area_idx = 0;
				GlobalVars->areah_idx = 2;
				GlobalVars->disk_params = 3;
				GlobalVars->num_dense = 1;
				GlobalVars->convertarea_flag = 0;
				GlobalVars->num_forcings = 1;
				break;
//--------------------------------------------------------------------------------------------
		case 22:	GlobalVars->dim = GlobalVars->problem_dim = 4;
				GlobalVars->template_flag = 0;
				GlobalVars->assim_flag = 0;
				GlobalVars->diff_start = 1;
				GlobalVars->no_ini_start = 2;
				num_global_params = 4;
				GlobalVars->uses_dam = 1;
				GlobalVars->params_size = 10;		//Need 1 extra for orifice_area
				GlobalVars->iparams_size = 0;
				GlobalVars->dam_params_size = 18;
				GlobalVars->area_idx = 0;
				GlobalVars->areah_idx = 2;
				GlobalVars->disk_params = 6;
				GlobalVars->num_dense = 1;
				GlobalVars->convertarea_flag = 0;
				GlobalVars->num_forcings = 1;
				break;
//--------------------------------------------------------------------------------------------
		case 23:	GlobalVars->dim = GlobalVars->problem_dim = 4;
				GlobalVars->template_flag = 0;
				GlobalVars->assim_flag = 0;
				GlobalVars->diff_start = 1;
				GlobalVars->no_ini_start = 2;
				num_global_params = 4;
				GlobalVars->uses_dam = 1;
				GlobalVars->params_size = 10;		//Need 1 extra for orifice_area
				GlobalVars->iparams_size = 0;
				GlobalVars->dam_params_size = 18;
				GlobalVars->area_idx = 0;
				GlobalVars->areah_idx = 2;
				GlobalVars->disk_params = 6;
				GlobalVars->num_dense = 1;
				GlobalVars->convertarea_flag = 0;
				GlobalVars->num_forcings = 1;
				break;
//--------------------------------------------------------------------------------------------
		case 30:	GlobalVars->dim = GlobalVars->problem_dim = 4;
				GlobalVars->template_flag = 0;
				GlobalVars->assim_flag = 0;
				GlobalVars->diff_start = 0;
				GlobalVars->no_ini_start = GlobalVars->dim;
				num_global_params = 7;
				GlobalVars->uses_dam = 0;
				GlobalVars->params_size = 21;
				GlobalVars->iparams_size = 0;
				GlobalVars->dam_params_size = 0;
				GlobalVars->area_idx = 2;
				GlobalVars->areah_idx = 1;
				GlobalVars->disk_params = 12;
				GlobalVars->num_dense = 1;
				GlobalVars->convertarea_flag = 0;
				GlobalVars->num_forcings = 2;
				break;
//--------------------------------------------------------------------------------------------
		case 40:	GlobalVars->dim = GlobalVars->problem_dim = 4;
				GlobalVars->template_flag = 0;
				GlobalVars->assim_flag = 0;
				GlobalVars->diff_start = 1;
				GlobalVars->no_ini_start = 2;
				num_global_params = 4;
				GlobalVars->uses_dam = 1;
				GlobalVars->params_size = 9;
				GlobalVars->iparams_size = 0;
				GlobalVars->dam_params_size = 6;
				GlobalVars->area_idx = 0;
				GlobalVars->areah_idx = 2;
				GlobalVars->disk_params = 6;
				GlobalVars->num_dense = 1;
				GlobalVars->convertarea_flag = 0;
				GlobalVars->num_forcings = 1;
				break;
//--------------------------------------------------------------------------------------------
		case 60:	GlobalVars->dim = GlobalVars->problem_dim = 3;
				GlobalVars->template_flag = 0;
				GlobalVars->assim_flag = 0;
				GlobalVars->diff_start = 0;
				GlobalVars->no_ini_start = GlobalVars->dim;
				num_global_params = 6;
				GlobalVars->uses_dam = 0;
				GlobalVars->params_size = 8;
				GlobalVars->iparams_size = 0;
				GlobalVars->dam_params_size = 0;
				GlobalVars->area_idx = 0;
				GlobalVars->areah_idx = 2;
				GlobalVars->disk_params = 4;
				GlobalVars->num_dense = 1;
				GlobalVars->convertarea_flag = 0;
				GlobalVars->num_forcings = 1;
				break;
//--------------------------------------------------------------------------------------------
		case 101:	GlobalVars->dim = GlobalVars->problem_dim = 3;
				GlobalVars->template_flag = 0;
				GlobalVars->assim_flag = 0;
				GlobalVars->diff_start = 0;
				GlobalVars->no_ini_start = GlobalVars->dim;
				num_global_params = 0;
				GlobalVars->num_dense = 3;
				GlobalVars->convertarea_flag = 1;
				GlobalVars->num_forcings = 1;
				break;
//--------------------------------------------------------------------------------------------
		case 105:	GlobalVars->dim = GlobalVars->problem_dim = 2;
				GlobalVars->template_flag = 0;
				GlobalVars->assim_flag = 0;
				GlobalVars->diff_start = 0;
				GlobalVars->no_ini_start = GlobalVars->dim;
				num_global_params = 0;
				GlobalVars->uses_dam = 0;
				GlobalVars->params_size = 16;
				GlobalVars->iparams_size = 0;
				GlobalVars->dam_params_size = 0;
				GlobalVars->area_idx = 2;
				GlobalVars->areah_idx = 1;
				GlobalVars->disk_params = 12;
				GlobalVars->num_dense = 1;
				GlobalVars->convertarea_flag = 1;
				GlobalVars->num_forcings = 1;
				break;
//--------------------------------------------------------------------------------------------
		case 190:	GlobalVars->dim = GlobalVars->problem_dim = 3;
				GlobalVars->template_flag = 0;
				GlobalVars->assim_flag = 0;
				GlobalVars->diff_start = 0;
				GlobalVars->no_ini_start = GlobalVars->dim;
				num_global_params = 6;
				GlobalVars->uses_dam = 0;
				GlobalVars->params_size = 8;
				GlobalVars->iparams_size = 0;
				GlobalVars->dam_params_size = 0;
				GlobalVars->area_idx = 0;
				GlobalVars->areah_idx = 2;
				GlobalVars->disk_params = 3;
				GlobalVars->num_dense = 1;
				GlobalVars->convertarea_flag = 0;
				GlobalVars->num_forcings = 2;
				break;
//--------------------------------------------------------------------------------------------
		case 200:	GlobalVars->dim = GlobalVars->problem_dim = 2;
				GlobalVars->template_flag = 0;
				GlobalVars->assim_flag = 0;
				GlobalVars->diff_start = 0;
				GlobalVars->no_ini_start = 1;
				num_global_params = 10;
				GlobalVars->uses_dam = 0;
				GlobalVars->params_size = 20;
				GlobalVars->iparams_size = 0;
				GlobalVars->dam_params_size = 0;
				GlobalVars->area_idx = 2;
				GlobalVars->areah_idx = 1;
				GlobalVars->disk_params = 12;
				GlobalVars->num_dense = 2;
				GlobalVars->convertarea_flag = 1;
				GlobalVars->num_forcings = 1;
				break;
//--------------------------------------------------------------------------------------------
		case 219:	GlobalVars->dim = GlobalVars->problem_dim = 1;
				GlobalVars->template_flag = 0;
				GlobalVars->assim_flag = 0;
				GlobalVars->diff_start = 0;
				GlobalVars->no_ini_start = GlobalVars->dim;
				num_global_params = 3;
				GlobalVars->uses_dam = 0;
				GlobalVars->params_size = 5;
				GlobalVars->iparams_size = 0;
				GlobalVars->dam_params_size = 0;
				GlobalVars->area_idx = 0;
				GlobalVars->areah_idx = 2;
				GlobalVars->disk_params = 3;
				GlobalVars->num_dense = 1;
				GlobalVars->convertarea_flag = 0;
				GlobalVars->num_forcings = 1;
				break;
//--------------------------------------------------------------------------------------------
		case 250:	GlobalVars->dim = GlobalVars->problem_dim = 3;
				GlobalVars->template_flag = 0;
				GlobalVars->assim_flag = 0;
				GlobalVars->diff_start = 0;
				GlobalVars->no_ini_start = GlobalVars->dim;
				num_global_params = 9;
				GlobalVars->uses_dam = 0;
				GlobalVars->params_size = 9;
				GlobalVars->iparams_size = 0;
				GlobalVars->dam_params_size = 0;
				GlobalVars->area_idx = 0;
				GlobalVars->areah_idx = 2;
				GlobalVars->disk_params = 4;
				GlobalVars->num_dense = 1;
				GlobalVars->convertarea_flag = 0;
				GlobalVars->num_forcings = 1;
				break;
//--------------------------------------------------------------------------------------------
		case 252:	GlobalVars->dim = GlobalVars->problem_dim = 4;
				GlobalVars->template_flag = 0;
				GlobalVars->assim_flag = 0;
				GlobalVars->diff_start = 0;
				GlobalVars->no_ini_start = GlobalVars->dim;
				num_global_params = 11;
				GlobalVars->uses_dam = 0;
				GlobalVars->params_size = 8;
				GlobalVars->iparams_size = 0;
				GlobalVars->dam_params_size = 0;
				GlobalVars->area_idx = 0;
				GlobalVars->areah_idx = 2;
				GlobalVars->disk_params = 3;
				GlobalVars->num_dense = 1;
				GlobalVars->convertarea_flag = 0;
				GlobalVars->num_forcings = 2;
				break;
//--------------------------------------------------------------------------------------------
		case 253:	GlobalVars->dim = GlobalVars->problem_dim = 4;
				GlobalVars->template_flag = 0;
				GlobalVars->assim_flag = 0;
				GlobalVars->diff_start = 0;
				GlobalVars->no_ini_start = GlobalVars->dim;
				num_global_params = 11;
				GlobalVars->uses_dam = 0;
				GlobalVars->params_size = 8;
				GlobalVars->iparams_size = 0;
				GlobalVars->dam_params_size = 0;
				GlobalVars->area_idx = 0;
				GlobalVars->areah_idx = 2;
				GlobalVars->disk_params = 3;
				GlobalVars->num_dense = 1;
				GlobalVars->convertarea_flag = 0;
				GlobalVars->num_forcings = 3;
				break;
//--------------------------------------------------------------------------------------------
		case 260:	GlobalVars->dim = GlobalVars->problem_dim = 4;
				GlobalVars->template_flag = 0;
				GlobalVars->assim_flag = 0;
				GlobalVars->diff_start = 0;
				GlobalVars->no_ini_start = GlobalVars->dim;
				num_global_params = 11;
				GlobalVars->uses_dam = 0;
				GlobalVars->params_size = 6;
				GlobalVars->iparams_size = 0;
				GlobalVars->dam_params_size = 0;
				GlobalVars->area_idx = 0;
				GlobalVars->areah_idx = 2;
				GlobalVars->disk_params = 3;
				GlobalVars->num_dense = 1;
				GlobalVars->convertarea_flag = 0;
				GlobalVars->num_forcings = 2;
				break;
//--------------------------------------------------------------------------------------------
		case 300:	GlobalVars->dim = 2;		//Make this a little bigger for the variational errors
				GlobalVars->problem_dim = 1;
				GlobalVars->template_flag = 0;
				GlobalVars->assim_flag = 1;
				GlobalVars->diff_start = 0;
				GlobalVars->no_ini_start = 1;
				num_global_params = 6;
				GlobalVars->uses_dam = 0;
				GlobalVars->params_size = 20;
				GlobalVars->iparams_size = 1;
				GlobalVars->dam_params_size = 0;
				GlobalVars->area_idx = 2;
				GlobalVars->areah_idx = 1;
				GlobalVars->disk_params = 12;
				GlobalVars->num_dense = 1;
				GlobalVars->convertarea_flag = 0;
				GlobalVars->num_forcings = 1;
				break;
//--------------------------------------------------------------------------------------------
		case 301:	GlobalVars->dim = 3;		//Make this a little bigger for the variational errors
				GlobalVars->problem_dim = 2;
				GlobalVars->template_flag = 0;
				GlobalVars->assim_flag = 1;
				GlobalVars->diff_start = 0;
				GlobalVars->no_ini_start = 2;
				num_global_params = 6;
				GlobalVars->uses_dam = 0;
				GlobalVars->params_size = 20;
				GlobalVars->iparams_size = 1;
				GlobalVars->dam_params_size = 0;
				GlobalVars->area_idx = 2;
				GlobalVars->areah_idx = 1;
				GlobalVars->disk_params = 12;
				GlobalVars->num_dense = 1;
				GlobalVars->convertarea_flag = 0;
				GlobalVars->num_forcings = 1;
				break;
//--------------------------------------------------------------------------------------------
		case 315:	GlobalVars->dim = 3;
				GlobalVars->problem_dim = 2;
				GlobalVars->template_flag = 0;
				GlobalVars->assim_flag = 1;
				GlobalVars->diff_start = 0;
				GlobalVars->no_ini_start = GlobalVars->dim;
				num_global_params = 6;
				GlobalVars->uses_dam = 0;
				GlobalVars->params_size = 20;
				GlobalVars->iparams_size = 1;
				GlobalVars->dam_params_size = 0;
				GlobalVars->area_idx = 2;
				GlobalVars->areah_idx = 1;
				GlobalVars->disk_params = 12;
				GlobalVars->num_dense = 1;
				GlobalVars->convertarea_flag = 0;
				GlobalVars->num_forcings = 1;
				break;
//--------------------------------------------------------------------------------------------
		case 2000: 	GlobalVars->dim = GlobalVars->problem_dim = 2;
				GlobalVars->template_flag = 1;
				GlobalVars->assim_flag = 0;
				GlobalVars->diff_start = 0;
				GlobalVars->no_ini_start = GlobalVars->dim;
				num_global_params = 6;
				GlobalVars->uses_dam = 0;
				GlobalVars->params_size = 20;
				GlobalVars->iparams_size = 0;
				GlobalVars->dam_params_size = 0;
				GlobalVars->area_idx = 2;
				GlobalVars->areah_idx = 1;
				GlobalVars->disk_params = 12;
				GlobalVars->num_dense = 1;
				GlobalVars->convertarea_flag = 1;
				GlobalVars->num_forcings = 1;
				break;
//--------------------------------------------------------------------------------------------
		default:	printf("Error: Invalid type (%u) in SetParamSizes.\n",type);
				abort();
//--------------------------------------------------------------------------------------------
	}

	//Create the dense_indices list
	//Note: If assim_flag is set, the variational eqs stuff is set outside here
	GlobalVars->dense_indices = (unsigned int*) malloc(GlobalVars->num_dense * sizeof(unsigned int));
	if(type == 21 || type == 22 || type == 23 || type == 40)
	{
		GlobalVars->dense_indices[0] = 1;	//Storage
	}
	else	//For all other types, just needed the first entry
	{
		GlobalVars->dense_indices[0] = 0;	//Discharge
	}

	//Make sure the appropriate number of global parameters are given
	if(GlobalVars->global_params->dim < num_global_params)
	{
		printf("\nError: Obtained %u parameters from .gbl file. Expected %u for model type %hu.\n",GlobalVars->global_params->dim,num_global_params,type);
		abort();
	}
	if(GlobalVars->global_params->dim > num_global_params)
		printf("\nWarning: Obtained %u parameters from .gbl file. Expected %u for model type %hu.\n",GlobalVars->global_params->dim,num_global_params,type);
}


//Performs some unit conversions on the data in params. This takes place immediately after reading in the DEM data,
//so these changes are available in all routines of user.c, if params is available. Note that dam data
//and precalculations are not available here.
//VEC* params:		Vector of parameters to convert.
//unsigned int type:	The index of the model.
void ConvertParams(VEC* params,unsigned int type)
{
	if(type == 19)
	{
		params->ve[1] *= 1000;	//L: km -> m
		params->ve[2] *= 1e6;	//A_h: km^2 -> m^2
	}
	else if(type == 190)
	{
		params->ve[1] *= 1000;	//L: km -> m
		params->ve[2] *= 1e6;	//A_h: km^2 -> m^2
	}
	else if(type == 20)
	{
		//params->ve[0] *= 1e6;	//km^2 -> m^2
		params->ve[1] *= 1000;	//km -> m
		params->ve[2] *= 1e6;	//km^2 -> m^2
	}
	else if(type == 60)
	{
		params->ve[1] *= 1000;	//L: km -> m
		params->ve[2] *= 1e6;	//A_h: km^2 -> m^2
		params->ve[3] *= 1e-3;	//h_b: mm->m
	}
	else if(type == 21)
	{
		//params->ve[0] *= 1e6;	//km^2 -> m^2
		params->ve[1] *= 1000;	//km -> m
		params->ve[2] *= 1e6;	//km^2 -> m^2
	}
	else if(type == 22 || type == 23 || type == 40)
	{
		//params->ve[0] *= 1e6;	//km^2 -> m^2
		params->ve[1] *= 1000;	//km -> m
		params->ve[2] *= 1e6;	//km^2 -> m^2
	}
	else if(type <= 5)
	{
		params->ve[0] *= 1000;	//km -> m
		params->ve[3] *= .001;	//mm -> m
		params->ve[4] *= .001;	//mm -> m
	}
	else if(type == 6)
	{
		params->ve[0] *= 1000;	//km -> m
		params->ve[3] *= .001;	//mm -> m
	}
	else if(type == 15 || type == 315)
	{
		params->ve[0] *= 1000;	//L: km -> m
		params->ve[3] *= .001;	//h_b: mm -> m
		params->ve[4] *= .001;	//h_H: mm -> m
	}
	else if(type == 30)
	{
		params->ve[0] *= 1000;		//L_h:  km -> m
		params->ve[4] *= .001;		//H_h:  mm -> m
		params->ve[5] *= 1000.0; 	//MaxInfRate:  m/hr -> mm/hr
	}
	else if(type == 105)
	{
		params->ve[0] *= 1000;	//km -> m
		params->ve[3] *= .001;	//mm -> m
		params->ve[4] *= .001;	//mm -> m
	}
	else if(type == 200)	//!!!! Did I screw these up on accident?!!!!
	{
		params->ve[0] *= 1000;	//L_h:  km -> m
/*
		//params->ve[3] *= .001;	//mm -> m
		params->ve[4] *= .001;	//H_h:  mm -> m
		params->ve[5] *= 1000.0; //MaxInfRate:  m/hr -> mm/hr
*/
	}
	if(type == 219)
	{
		params->ve[1] *= 1000;	//L: km -> m
		params->ve[2] *= 1e6;	//A_h: km^2 -> m^2
	}
	else if(type == 250)
	{
		params->ve[1] *= 1000;		//L_h: km -> m
		params->ve[2] *= 1e6;		//A_h: km^2 -> m^2
		params->ve[4] *= .001;		//H_h: mm -> m
	}
	else if(type == 252 || type == 260)
	{
		params->ve[1] *= 1000;		//L_h: km -> m
		params->ve[2] *= 1e6;		//A_h: km^2 -> m^2
		//params->ve[4] *= .001;		//H_h: mm -> m
	}
	else if(type == 253)
	{
		params->ve[1] *= 1000;		//L_h: km -> m
		params->ve[2] *= 1e6;		//A_h: km^2 -> m^2
	}
	else if(type == 300 || type == 301)
	{
		params->ve[0] *= 1000;	//km -> m
		params->ve[3] *= .001;	//mm -> m
		params->ve[4] *= .001;	//mm -> m
	}
	else if(type == 2000)
	{
		params->ve[0] *= 1000;	//km -> m
		params->ve[3] *= .001;	//mm -> m
		params->ve[4] *= .001;	//mm -> m
	}
}


//Sets the system of ODEs and the Runge-Kutta solver for link. This method MUST set both link->f
//	and link->RKSolver. The Jacobian of f (link->Jacobian) may be set here, if using an
//	implicit solver.
//Link* link: 		The link at which the ODEs and Runge-Kutta solver are selected.
//unsigned int type: 	The index of the model to be set.
//unsigned int exp_imp: 0 if using an explicit solver, 1 if implicit.
//unsigned int dam: 	0 if no dam is present at link, 1 if a dam is present.
void InitRoutines(Link* link,unsigned int type,unsigned int exp_imp,unsigned short int dam)
{
	//Select appropriate RK Solver for the numerical method (link->RKSolver)
//	if(type == 19)
//		link->RKSolver = &ExplicitRKSolverDiscont;
	if( (type == 21 || type == 22 || type == 23 || type == 40) && dam == 1)
		link->RKSolver = &ExplicitRKIndex1SolverDam;
	else if( (type == 21 || type == 22 || type ==23 || type == 40) && dam == 0)
		link->RKSolver = &ExplicitRKIndex1Solver;
	else if(exp_imp == 0)
		link->RKSolver = &ExplicitRKSolver;
//	else if(link->method->exp_imp == 1)
//		link->RKSolver = &RadauRKSolver;
	else
		printf("Warning: No solver selected for link ID %u.\n",link->ID);

	//Select the RHS function of the ODE (link->f, link->Jacobian)
	if(type == 0)
	{
		link->f = &simple_river;
		link->Jacobian = &Jsimple;
		link->alg = NULL;
		link->state_check = NULL;
		link->CheckConsistency = &CheckConsistency_Nonzero_1States;
	}
	else if(type == 1 || type == 3)
	{
		link->f = &river_rainfall;
		link->alg = NULL;
		link->state_check = NULL;
		link->CheckConsistency = &CheckConsistency_Nonzero_2States;
	}
	else if(type == 2)
	{
		link->f = &simple_hillslope;
		link->alg = NULL;
		link->state_check = NULL;
		link->CheckConsistency = &CheckConsistency_Nonzero_2States;
	}
	else if(type == 4)
	{
		link->f = &simple_soil;
		link->Jacobian = &Jsimple_soil;
		link->alg = NULL;
		link->state_check = NULL;
		link->CheckConsistency = &CheckConsistency_Nonzero_4States;
	}
	else if(type == 5)
	{
		link->f = &soil_rainfall;
		link->alg = NULL;
		link->state_check = NULL;
		link->CheckConsistency = &CheckConsistency_Model5;
	}
	else if(type == 6)
	{
		link->f = &qsav_rainfall;
		link->alg = NULL;
		link->state_check = NULL;
		link->CheckConsistency = &CheckConsistency_Nonzero_4States;
	}
	else if(type == 15)
	{
		link->f = &river_rainfall_adjusted;
		link->alg = NULL;
		link->state_check = NULL;
		link->CheckConsistency = &CheckConsistency_Nonzero_2States;
	}
	else if(type == 19)
	{
		link->f = &LinearHillslope_Evap;
		link->alg = NULL;
		link->state_check = NULL;
		link->CheckConsistency = &CheckConsistency_Nonzero_3States;
	}
	else if(type == 20)
	{
		link->f = &Hillslope_Toy;
		link->alg = NULL;
		link->state_check = NULL;
		link->CheckConsistency = &CheckConsistency_Nonzero_3States;
	}
	else if(type == 21)
	{
		if(dam)	link->f = &dam_rain_hillslope;
		else	link->f = &nodam_rain_hillslope;
		link->alg = &dam_q;
		link->state_check = &dam_check;
		link->CheckConsistency = &CheckConsistency_Nonzero_4States;
	}
	else if(type == 22)
	{
		if(dam)	link->f = &dam_rain_hillslope2;
		else	link->f = &nodam_rain_hillslope2;
		link->alg = &dam_q2;
		link->state_check = &dam_check2;
		link->CheckConsistency = &CheckConsistency_Nonzero_4States;
	}
	else if(type == 23)
	{
		if(dam)	link->f = &dam_rain_hillslope3;
		else	link->f = &nodam_rain_hillslope3;
		link->alg = &dam_q3;
		link->state_check = &dam_check3;
		link->CheckConsistency = &CheckConsistency_Nonzero_4States;
	}
	else if(type == 30)
	{
		link->f = &lcuencas_soilrain;
		link->alg = NULL;
		link->state_check = NULL;
		link->CheckConsistency = &CheckConsistency_Model30;
	}
	else if(type == 40)
	{
		if(dam)	link->f = &dam_rain_hillslope_qsv;
		else	link->f = &nodam_rain_hillslope_qsv;
		link->alg = &dam_q_qvs;
		link->state_check = &dam_check_qvs;
		link->CheckConsistency = &CheckConsistency_Nonzero_4States;
	}
	else if(type == 60)
	{
		link->f = &LinearHillslope_Evap_RC;
		link->alg = NULL;
		link->state_check = NULL;
		link->CheckConsistency = &CheckConsistency_Nonzero_3States;
	}
	else if(type == 101)
	{
		link->f = &Robertson;
		link->Jacobian = &JRobertson;
		link->alg = NULL;
		link->state_check = NULL;
		link->CheckConsistency = &CheckConsistency_Nonzero_3States;
	}
	else if(type == 105)
	{
		link->f = &river_rainfall_summary;
		link->alg = NULL;
		link->state_check = NULL;
		link->CheckConsistency = &CheckConsistency_Nonzero_2States;
	}
	else if(type == 190)
	{
		link->f = &LinearHillslope_MonthlyEvap;
		link->alg = NULL;
		link->state_check = NULL;
		link->CheckConsistency = &CheckConsistency_Nonzero_3States;
	}
	else if(type == 200)	//This is for use with SIMPLE only
	{
		link->f = NULL;
		link->alg = NULL;
		link->state_check = NULL;
		link->CheckConsistency = &CheckConsistency_Nonzero_2States;
	}
	else if(type == 219)
	{
		link->f = &Tiling;
		link->alg = NULL;
		link->state_check = NULL;
		link->CheckConsistency = &CheckConsistency_Nonzero_1States;
	}
	else if(type == 250)
	{
		link->f = &NonLinearHillslope;
		link->alg = NULL;
		link->state_check = NULL;
		link->CheckConsistency = &CheckConsistency_Nonzero_3States;
	}
	else if(type == 252)
	{
		link->f = &TopLayerHillslope;
		link->alg = NULL;
		link->state_check = NULL;
		link->CheckConsistency = &CheckConsistency_Nonzero_4States;
	}
	else if(type == 253)
	{
		//if(link->ID == 318208)
		if(link->res)
		{
			link->f = &TopLayerHillslope_Reservoirs;
			link->RKSolver = &ForcedSolutionSolver;
		}
		else			link->f = &TopLayerHillslope;
		link->alg = NULL;
		link->state_check = NULL;
		link->CheckConsistency = &CheckConsistency_Nonzero_4States;
	}
	else if(type == 260)
	{
		link->f = &TopLayerNonlinearExp;
		link->alg = NULL;
		link->state_check = NULL;
		link->CheckConsistency = &CheckConsistency_Nonzero_4States;
	}
	else if(type == 300)
	{
		link->f = &assim_simple_river;
		link->alg = NULL;
		link->state_check = NULL;
		link->CheckConsistency = &CheckConsistency_Nonzero_1States;
	}
	else if(type == 301)
	{
		link->f = &assim_river_rainfall;
		link->alg = NULL;
		link->state_check = NULL;
		link->CheckConsistency = &CheckConsistency_Nonzero_2States;
	}
	else if(type == 315)
	{
		link->f = &assim_river_rainfall_adjusted;
		link->alg = NULL;
		link->state_check = NULL;
		link->CheckConsistency = &CheckConsistency_Nonzero_2States;
	}
/*
	else if(type == 2000)
	{
		link->f = &parser_test;
		link->alg = NULL;
		link->state_check = NULL;
	}
*/
	else
		printf("Warning: No ODE selected for link ID %u.\n",link->ID);
}


//Perform precalculations needed for the differential equation.  These should be stored in params after the DEM
//	data and after the dam data (i.e. params->ve[disk_params] is the first precalcuation, params->ve[params_size]
//	is the first dam datum). This method is run at each link. This method can do nothing, if no precalculations are
//	expected for the model. This assumes the same number of precalculations regardless if there is a dam or not.
//VEC* global_params:		Vector of the global parameters of the system. These are already set and are available for use.
//VEC* params:			Vector of the parameters at a link. The DEM data and dam data are already set (but may be
//					altered here). Only the entries for precalculations need to be set.
//IVEC* iparams:		Integer vector of parameters at a link. This should be set entirely here.
//unsigned int disk_params:	The first entry of params that should be set here.
//unsigned int params_size:	First entry of the dam data. Don't change this entry or later unless you want to modify the dam!
//unsigned int type:		The index of the model.
void Precalculations(Link* link_i,VEC* global_params,VEC* params,IVEC* iparams,unsigned int disk_params,unsigned int params_size,unsigned short int dam,unsigned int type)
{
	if(type == 19)
	{
		//Order of parameters: A_i,L_i,A_h,k2,k3,invtau,c_1,c_2
		//The numbering is:	0   1   2   3  4    5    6   7
		//Order of global_params: v_r,lambda_1,lambda_2,RC,v_h,v_g,e_pot
		//The numbering is:        0      1        2     3  4   5   6
		double* vals = params->ve;
		double A_i = params->ve[0];
		double L_i = params->ve[1];
		double A_h = params->ve[2];
		double v_r = global_params->ve[0];
		double lambda_1 = global_params->ve[1];
		double lambda_2 = global_params->ve[2];
		double RC = global_params->ve[3];
		double v_h = global_params->ve[4];
		double v_g = global_params->ve[5];

		vals[3] = v_h * L_i / A_h * 60.0;	//[1/min]  k2
		vals[4] = v_g * L_i / A_h * 60.0;	//[1/min]  k3
		vals[5] = 60.0*v_r*pow(A_i,lambda_2) / ((1.0-lambda_1)*L_i);	//[1/min]  invtau
		vals[6] = RC*(0.001/60.0);		//(mm/hr->m/min)  c_1
		vals[7] = (1.0-RC)*(0.001/60.0);	//(mm/hr->m/min)  c_2
	}
	else if(type == 190)
	{
		//Order of parameters: A_i,L_i,A_h,k2,k3,invtau,c_1,c_2
		//The numbering is:	0   1   2   3  4    5    6   7
		//Order of global_params: v_r,lambda_1,lambda_2,RC,v_h,v_g
		//The numbering is:        0      1        2     3  4   5 
		double* vals = params->ve;
		double A_i = params->ve[0];
		double L_i = params->ve[1];
		double A_h = params->ve[2];
		double v_r = global_params->ve[0];
		double lambda_1 = global_params->ve[1];
		double lambda_2 = global_params->ve[2];
		double RC = global_params->ve[3];
		double v_h = global_params->ve[4];
		double v_g = global_params->ve[5];

		vals[3] = v_h * L_i / A_h * 60.0;	//[1/min]  k2
		vals[4] = v_g * L_i / A_h * 60.0;	//[1/min]  k3
		vals[5] = 60.0*v_r*pow(A_i,lambda_2) / ((1.0-lambda_1)*L_i);	//[1/min]  invtau
		vals[6] = RC*(0.001/60.0);		//(mm/hr->m/min)  c_1
		vals[7] = (1.0-RC)*(0.001/60.0);	//(mm/hr->m/min)  c_2
	}
	else if(type == 20)
	{
		//Order of parameters: A_i,L_i,A_h,invtau,c_1,c_2
		//The numbering is:	0   1   2    3     4   5
		//Order of global_params: v_r,lambda_1,lambda_2,beta,k_p,k_a,theta_p,theta_a,scale_p,scale_a
		//The numbering is:        0      1        2     3    4   5     6       7      8        9
		double* vals = params->ve;
		double A_i = params->ve[0];
		double L_i = params->ve[1];
		double A_h = params->ve[2];
		double v_r = global_params->ve[0];
		double lambda_1 = global_params->ve[1];
		double lambda_2 = global_params->ve[2];


		vals[3] = 60.0*v_r*pow(A_i,lambda_2) / ((1.0-lambda_1)*L_i); //[1/min]  invtau

/*
		//For gamma
		double k_p = global_params->ve[4];
		double k_a = global_params->ve[5];
		double theta_p = global_params->ve[6];
		double theta_a = global_params->ve[7];
		vals[4] = 1.0/(tgamma(k_p)*pow(theta_p,k_p)); //c_1
		vals[5] = 1.0/(tgamma(k_a)*pow(theta_a,k_a)); //c_2
*/
/*
		//For log normal
		double mu = global_params->ve[3];
		double sigma2 = global_params->ve[4];
		double scale = global_params->ve[5];
		vals[4] = scale/(pow(2.0*3.141592653589*sigma2,0.5));
		vals[5] = 0.0;
*/

		//Set some random values
		//vals[4] = (double)(rand()%996) / 10.0 + 0.5;	//0.5 to 100.0
		//vals[5] = (double)(rand()%991) + 10.0;	//10 to 1000

/*
		//Order of parameters: A_i,L_i,A_h,k2,k3,invtau,c_1,c_2
		//The numbering is:	0   1   2   3  4    5    6   7
		//Order of global_params: v_r,lambda_1,lambda_2,RC,v_h,v_g
		//The numbering is:        0      1        2     3  4   5 
		double* vals = params->ve;
		double A_i = params->ve[0];
		double L_i = params->ve[1];
		double A_h = params->ve[2];
		double v_r = global_params->ve[0];
		double lambda_1 = global_params->ve[1];
		double lambda_2 = global_params->ve[2];
		double RC = global_params->ve[3];
		double v_h = global_params->ve[4];
		double v_g = global_params->ve[5];
		double F_et = 0.05;

		vals[3] = v_h * L_i / A_h * 60.0;	//[1/min]  k2
		vals[4] = v_g * L_i / A_h * 60.0;	//[1/min]  k3
		vals[5] = 60.0*v_r*pow(A_i,lambda_2) / ((1.0-lambda_1)*L_i);	//[1/min]  invtau
		vals[6] = RC*(0.001/60.0);		//(mm/hr->m/min)  c_1
		vals[7] = F_et*(1.0-RC)*(0.001/60.0);	//(mm/hr->m/min)  c_2
*/
	}
	else if(type == 60)
	{
		//Order of parameters: A_i,L_i,A_h,h_b,k2,k3,invtau,c_1
		//The numbering is:	0   1   2   3   4  5    6    7 
		//Order of global_params: v_r,lambda_1,lambda_2,v_h,v_g,e_pot
		//The numbering is:        0      1        2     3   4   5  
		double* vals = params->ve;
		double A_i = params->ve[0];
		double L_i = params->ve[1];
		double A_h = params->ve[2];
		double v_r = global_params->ve[0];
		double lambda_1 = global_params->ve[1];
		double lambda_2 = global_params->ve[2];
		double v_h = global_params->ve[3];
		double v_g = global_params->ve[4];

		vals[4] = v_h * L_i / A_h * 60.0;	//[1/min]  k2
		vals[5] = v_g * L_i / A_h * 60.0;	//[1/min]  k3
		vals[6] = 60.0*v_r*pow(A_i,lambda_2) / ((1.0-lambda_1)*L_i);	//[1/min]  invtau
		vals[7] = (0.001/60.0);		//(mm/hr->m/min)  c_1
	}
	else if(type == 21)
	{
		//Order of parameters: A_i,L_i,A_h,k2,k3,invtau,orifice_area,H_spill,H_max,S_max,alpha,orifice_diam,c_1,c_2,L_spill
		//The numbering is:	0   1   2  3  4    5	       6      7       8     9	  10	    11       12  13  14
		//Order of global_params: v_r,lambda_1,lambda_2,RC,S_0,v_h,v_g
		//The numbering is:        0      1        2     3  4   5   6
		//Need to set entries 3, 4, 5, and 6 of params.
		double* vals = params->ve;
		double A_i = params->ve[0];
		double L_i = params->ve[1];
		double A_h = params->ve[2];
		double v_r = global_params->ve[0];
		double lambda_1 = global_params->ve[1];
		double lambda_2 = global_params->ve[2];
		double v_h = global_params->ve[5];
		double v_g = global_params->ve[6];

		vals[3] = v_h * L_i / A_h * 60.0;	//[1/min]  k2
		//vals[4] = vals[3] / 20.0;		//[1/min]  k3
		vals[4] = v_g * L_i / A_h * 60.0;	//[1/min]  k3
		vals[5] = pow(v_r*pow(A_i,lambda_2)/L_i,1.0/(1.0-lambda_1)) * 60.0; //[1/min]  invtau

		if(dam)		vals[6] = 3.1415926535897932 * vals[11]*vals[11] / 4.0;
	}
	else if(type == 22 || type == 23)
	{
		//Order of parameters: A_i,L_i,A_h,RC,v_h,v_r,k2,k3,invtau,orifice_area,H_spill,H_max,S_max,alpha,orifice_diam,c_1,c_2,L_spill
		//The numbering is:	0   1   2  3   4   5   6  7   8          9	  10	  11   12     13      14        15  16   17
		//Order of global_params: lambda_1,lambda_2,S_0,v_g
		//The numbering is:         0        1       2   3
		//Need to set entries 6, 7, 8, and 9 of params.
		double* vals = params->ve;
		double A_i = params->ve[0];
		double L_i = params->ve[1];
		double A_h = params->ve[2];
		double v_h = params->ve[4];
		double v_r = params->ve[5];
		double lambda_1 = global_params->ve[0];
		double lambda_2 = global_params->ve[1];
		double v_g = global_params->ve[3];

		vals[6] = v_h * L_i / A_h * 60.0;	//[1/min]  k2
		//vals[7] = vals[6] / 20.0;		//[1/min]  k3
		vals[7] = v_g * L_i / A_h * 60.0;	//[1/min]  k3
		vals[8] = pow(v_r*pow(A_i,lambda_2)/L_i,1.0/(1.0-lambda_1)) * 60.0; //[1/min]  invtau

		if(dam)		vals[9] = 3.1415926535897932 * vals[14]*vals[14] / 4.0; //orifice_area
	}
	else if(type == 40)
	{
		//Order of parameters: A_i,L_i,A_h,RC,v_h,v_r,k2,k3,invtau
		//The numbering is:	0   1   2  3   4   5   6  7   8
		//Order of global_params: lambda_1,lambda_2,S_0,v_g
		//The numbering is:         0        1       2   3
		//Need to set entries 6, 7, and 8 of params.
		double* vals = params->ve;
		double A_i = params->ve[0];
		double L_i = params->ve[1];
		double A_h = params->ve[2];
		double v_h = params->ve[4];
		double v_r = params->ve[5];
		double lambda_1 = global_params->ve[0];
		double lambda_2 = global_params->ve[1];
		double v_g = global_params->ve[3];

		vals[6] = v_h * L_i / A_h * 60.0;	//[1/min]  k2
		//vals[7] = vals[6] / 20.0;		//[1/min]  k3
		vals[7] = v_g * L_i / A_h * 60.0;	//[1/min]  k3
		vals[8] = pow(v_r*pow(A_i,lambda_2)/L_i,1.0/(1.0-lambda_1)) * 60.0; //[1/min]  invtau
	}
	else if(type <= 5 || type == 200)
	{
		//Order of parameters: L_i,A_h,A_i,h_b,h_H,max_inf_rate,K_sat,S_h,eta,b_H,c_H,d_H,invtau,epsilon,c_1,c_2,c_3,c_4,c_5,c_6
		//The numbering is:     0   1   2   3   4       5         6    7   8   9   10  11  12    13      14  15  16  17  18  19
		//Order of global_params: v_r,lambda_1,lambda_2,Q_r,A_r,RC, (for 200) u_0,U,DL,scale
		//The numbering is:        0      1        2     3   4   5   	       6  7  8   9
		//Need to set entries 12-19 of params.
		double* vals = params->ve;
		double K_T = 1.0;
		double s_r = 1.0;
		double rootS_h = pow(vals[7],.5);
		double v_r = global_params->ve[0];
		double lambda_1 = global_params->ve[1];
		double lambda_2 = global_params->ve[2];
		double Q_r = global_params->ve[3];
		double A_r = global_params->ve[4];
		double RC = global_params->ve[5];

		vals[12] = 60.0*v_r*pow(vals[2]/A_r,lambda_2)/((1.0-lambda_1)*vals[0]);
		vals[13] = vals[3] / s_r;
		vals[14] = 2.0 / .6 * vals[0] * s_r / Q_r * 1.0 / vals[8] * rootS_h * pow(vals[3],2.0/3.0);	//c_1
		vals[15] = vals[6] * vals[0] * vals[3] / (3600.0 * Q_r);
		vals[16] = 1e-3/(60*s_r) * RC;	//c_3
		vals[17] = 60 * 1e-6 * 2.0 / .6 * vals[0] / vals[1] * 1.0 / vals[8] * rootS_h;	//c_4
		vals[18] = K_T/60.0;
		vals[19] = vals[6]/(60.0*s_r);
	}
	else if(type == 6)
	{
		//Order of parameters: L_i,A_h,A_i,h_b,K_sat,K_sp,d_0,d_1,d_2,invbeta,alpha_N,alpha_soil,a_res,v_res,invtau,gamma,c_1,c_2,c_3,c_4
		//The numbering is:     0   1   2   3   4     5    6   7   8     9      10       11       12    13     14    15   16  17  18  19
		//Order of global_params: v_r,lambda_1,lambda_2,Q_r,A_r
		//The numbering is:        0      1        2     3   4
		//Need to set entries 14 - 19 and 9 of params.
		double v_r = global_params->ve[0];
		double lambda_1 = global_params->ve[1];
		double lambda_2 = global_params->ve[2];
		double Q_r = global_params->ve[3];
		double A_r = global_params->ve[4];

		double* vals = params->ve;
		double L_i = vals[0];
		double A_h = vals[1];
		double A_i = vals[2];
		double h_b = vals[3];
		double K_sat = vals[4];
		double K_sp = vals[5];
		double beta = vals[9];
		double alpha_N = vals[10];
		double alpha_soil = vals[11];

		vals[9] = 1.0 / beta;		//invbeta

		vals[14] = 60.0 * v_r * pow(A_i/A_r,lambda_2) / ( (1.0 - lambda_1) * L_i );	//invtau
		vals[15] = (1e6/60.0) * A_h / Q_r;	//gamma
		vals[16] = K_sp / (60.0 * h_b);	//c_1
		vals[17] = (1e-6/60.0) * alpha_soil * K_sat * L_i / (alpha_N * A_h); //c_2
		vals[18] = 1e-3/60.0; //c_3
		vals[19] = alpha_N * h_b; //c_4
	}
	else if(type == 15)
	{
		//Order of parameters: L_i,A_h,A_i,h_b,h_H,max_inf_rate,K_sat,S_h,eta,b_H,c_H,d_H,invtau,epsilon,c_1,c_2,c_3,c_4,c_5,c_6
		//The numbering is:     0   1   2   3   4       5         6    7   8   9   10  11  12    13      14  15  16  17  18  19
		//Order of global_params: v_r,lambda_1,lambda_2,v_h,A_r,RC
		//The numbering is:        0      1        2     3   4   5
		//Need to set entries 12-19 of params.
		double* vals = params->ve;
		double K_T = 1.0;
		double s_r = 1.0;
		double rootS_h = pow(vals[7],.5);
		double L = params->ve[0];
		double A_h = params->ve[1] * 1e6;	//Put into m^2
		double eta = params->ve[8];
		double v_r = global_params->ve[0];
		double lambda_1 = global_params->ve[1];
		double lambda_2 = global_params->ve[2];
		double v_h = global_params->ve[3];
		double A_r = global_params->ve[4];
		double RC = global_params->ve[5];


		vals[12] = 60.0*v_r*pow(vals[2]/A_r,lambda_2)/((1.0-lambda_1)*vals[0]);	//invtau
		vals[13] = vals[3] / s_r; //epsilon
		vals[14] = v_h*L;	//c_1 [m^2/s]
		vals[15] = vals[6] * vals[0] * vals[3] / 3600.0; //c_2
		vals[16] = (1e-3/60.0) * RC;	//c_3
		vals[17] = 60.0*v_h*L/A_h;	//c_4 [1/min], A_h converted above
		vals[18] = K_T/60.0;
		vals[19] = vals[6]/(60.0*s_r);

	}
	else if(type == 30)
	{
		//y_i has q (y_i[0]), s_p (y_i[1]), h_w (y_i[2]), theta (y_i[3])
		//Order of parameters: L_h,A_h,A_up,H_b,H_h,max_inf_rate,K_SAT,S_H,n_vh, b, c, d,K_Q,V_T,c_1,c_2,c_3,c_4,c_5,c_6,c_7
		//The numbering is:     0   1   2    3   4       5         6    7   8    9 10 11  12  13  14  15  16  17  18  19  20
		//Order of global_params: v_0,lambda_1,lambda_2,Q_r,A_r,K_T,C_r
		//The numbering is:        0      1        2     3   4   5   6
		//Need to set entries 12-19 of params.

		//Global params
		double v_0 = global_params->ve[0];
		double lambda_1 = global_params->ve[1];
		double lambda_2 = global_params->ve[2];
		double A_r = global_params->ve[4];
		double C_r = global_params->ve[6];

		//Local params
		double L_H = params->ve[0];
		double A_H = params->ve[1];
		double A_up = params->ve[2];
		double H_b = params->ve[3];
		double H_h = params->ve[4];
		double K_SAT = params->ve[6];
		double S_H = params->ve[7];
		double n_vh = params->ve[8];
		double H_relmax = H_h;


		params->ve[12] = 60.0*C_r*v_0*pow(A_up/A_r,lambda_2)/((1.0-lambda_1)*L_H);	//K_Q
		params->ve[13] = 1e3 * H_b * A_H;	//V_T
		params->ve[14] = 3.6e3 / n_vh * pow(S_H,.5);	//c_1
		//params->ve[14] = 3.6e6 / n_vh * pow(S_H,.5);	//c_1
		//params->ve[15] = (2e-6/.6) * (L_H/A_H);	//c_2
		params->ve[15] = (2e-6) * (L_H/A_H);	//c_2
		params->ve[16] = 1e3 / params->ve[13] * K_SAT * S_H;	//c_3
		params->ve[17] = (1.0/3.6) * A_H;	//c_4
		params->ve[18] = (1e6/60.0) * A_H / H_b;	//c_5
		params->ve[19] = (1e3 / 60.0) * A_H;	//c_6
		params->ve[20] = (H_relmax > 1e-12) ? 1e6 * A_H / H_relmax : 0.0;	//c_7
	}
	else if(type == 105)
	{
		//Order of parameters: L_i,A_h,A_i,h_b,h_H,max_inf_rate,K_sat,S_h,eta,b_H,c_H,d_H,invtau,c_1,c_2,c_3
		//The numbering is:     0   1   2   3   4       5         6    7   8   9   10  11  12    13  14  15 
		//Order of global_params: v_r,lambda_1,lambda_2,Q_r,A_r,RC
		//The numbering is:        0      1        2     3   4   5
		double* vals = params->ve;
		double v_r = global_params->ve[0];
		double lambda_1 = global_params->ve[1];
		double lambda_2 = global_params->ve[2];
		double A_r = global_params->ve[4];
		double RC = global_params->ve[5];
		double rootS_h = pow(vals[7],.5);

		//invtau
		vals[12] = 60.0*v_r*pow(vals[2]/A_r,lambda_2)/((1.0-lambda_1)*vals[0]);
		//c_1
		vals[13] = 2.0 / .6 * vals[0] * rootS_h / vals[8];
		//c_2
		vals[14] = 1e-3 / 60.0 * RC;
		//c_3
		vals[15] = 60.0 * 1e-6 * 2.0 / .6 * vals[0] / vals[1] * rootS_h / vals[8];
	}
	else if(type == 219)
	{
		//Order of parameters: A_i,L_i,A_h,invtau,c_1
		//The numbering is:	0   1   2    3     4
		//Order of global_params: v_0,lambda_1,lambda_2
		//The numbering is:        0      1        2
		double* vals = params->ve;
		double A_i = params->ve[0];
		double L_i = params->ve[1];
		double A_h = params->ve[2];
		double v_0 = global_params->ve[0];
		double lambda_1 = global_params->ve[1];
		double lambda_2 = global_params->ve[2];

		vals[3] = 60.0*v_0*pow(A_i,lambda_2) / ((1.0-lambda_1)*L_i);	//[1/min]  invtau
		vals[4] = (1e-3/3600.0) * A_h;		//c_1
	}
	if(type == 250)
	{
		//Order of parameters: A_i,L_i,A_h,h_r,invtau,k_2,k_I,c_1,c_2
		//The numbering is:	0   1   2   3    4     5   6   7   8
		//Order of global_params: v_0,lambda_1,lambda_2,v_h,k_3,k_I_factor,gamma,h_b,e_pot
		//The numbering is:        0      1        2     3   4     5         6    7	8
		double* vals = params->ve;
		double A_i = params->ve[0];
		double L_i = params->ve[1];
		double A_h = params->ve[2];

		double v_0 = global_params->ve[0];
		double lambda_1 = global_params->ve[1];
		double lambda_2 = global_params->ve[2];
		double v_h = global_params->ve[3];
		double k_I_factor = global_params->ve[5];

		vals[4] = 60.0*v_0*pow(A_i,lambda_2) / ((1.0-lambda_1)*L_i);	//[1/min]  invtau
		vals[5] = v_h * L_i/A_h * 60.0;	//[1/min] k_2
		vals[6] = vals[5] * k_I_factor;	//[1/min] k_I
		vals[7] = (0.001/60.0);		//(mm/hr->m/min)  c_1
		vals[8] = A_h/60.0;	//  c_2
	}
	else if(type == 252 || type == 253)
	{
		//Order of parameters: A_i,L_i,A_h,invtau,k_2,k_i,c_1,c_2
		//The numbering is:	0   1   2    3     4   5   6   7 
		//Order of global_params: v_0,lambda_1,lambda_2,v_h,k_3,k_I_factor,h_b,S_L,A,B,exponent
		//The numbering is:        0      1        2     3   4     5        6   7  8 9  10
		double* vals = params->ve;
		double A_i = params->ve[0];
		double L_i = params->ve[1];
		double A_h = params->ve[2];

		double v_0 = global_params->ve[0];
		double lambda_1 = global_params->ve[1];
		double lambda_2 = global_params->ve[2];
		double v_h = global_params->ve[3];
		double k_i_factor = global_params->ve[5];

		vals[3] = 60.0*v_0*pow(A_i,lambda_2) / ((1.0-lambda_1)*L_i);	//[1/min]  invtau
		vals[4] = v_h * L_i/A_h * 60.0;	//[1/min] k_2
		vals[5] = vals[4] * k_i_factor;	//[1/min] k_i
		vals[6] = (0.001/60.0);		//(mm/hr->m/min)  c_1
		vals[7] = A_h/60.0;	//  c_2
	}
	else if(type == 260)
	{
		//Order of parameters: A_i,L_i,A_h | invtau,c_1,c_2
		//The numbering is:	0   1   2  |    3    4   5 
		//Order of global_params: v_0,lambda_1,lambda_2,h_b,k_D,k_2,k_dry,k_i,T_L,N,phi
		//The numbering is:        0      1        2     3   4   5   6     7   8  9  10
		double* vals = params->ve;
		double A_i = params->ve[0];
		double L_i = params->ve[1];
		double A_h = params->ve[2];

		double v_0 = global_params->ve[0];
		double lambda_1 = global_params->ve[1];
		double lambda_2 = global_params->ve[2];

		vals[3] = 60.0*v_0*pow(A_i,lambda_2) / ((1.0-lambda_1)*L_i);	//[1/min]  invtau
		vals[4] = (0.001/60.0);		//(mm/hr->m/min)  c_1
		vals[5] = A_h/60.0;	//  c_2
	}
	else if(type == 300 || type == 301)
	{
		//Order of parameters: L_i,A_h,A_i,h_b,h_H,max_inf_rate,K_sat,S_h,eta,b_H,c_H,d_H,invtau,epsilon,c_1,c_2,c_3,c_4,c_5,c_6
		//The numbering is:     0   1   2   3   4       5         6    7   8   9   10  11  12    13      14  15  16  17  18  19
		//Order of global_params: v_r,lambda_1,lambda_2,Q_r,A_r,RC
		//The numbering is:        0      1        2     3   4   5
		//Need to set entries 12-19 of params.
		double* vals = params->ve;
		double K_T = 1.0;
		double s_r = 1.0;
		double rootS_h = pow(vals[7],.5);
		double v_r = global_params->ve[0];
		double lambda_1 = global_params->ve[1];
		double lambda_2 = global_params->ve[2];
		double Q_r = global_params->ve[3];
		double A_r = global_params->ve[4];
		double RC = global_params->ve[5];

		vals[12] = 60.0*v_r*pow(vals[2]/A_r,lambda_2)/((1.0-lambda_1)*vals[0]);
		vals[13] = vals[3] / s_r;
		vals[14] = 2.0 / .6 * vals[0] * s_r / Q_r * 1.0 / vals[8] * rootS_h * pow(vals[3],2.0/3.0);
		vals[15] = vals[6] * vals[0] * vals[3] / (3600.0 * Q_r);
		vals[16] = 1e-3/(60*s_r) * RC;
		vals[17] = 60 * 1e-6 * 2.0 / .6 * vals[0] / vals[1] * 1.0 / vals[8] * rootS_h;
		vals[18] = K_T/60.0;
		vals[19] = vals[6]/(60.0*s_r);

		iparams->ve[0] = link_i->location;
	}
	else if(type == 315)
	{
		//Order of parameters: L_i,A_h,A_i,h_b,h_H,max_inf_rate,K_sat,S_h,eta,b_H,c_H,d_H,invtau,epsilon,c_1,c_2,c_3,c_4,c_5,c_6
		//The numbering is:     0   1   2   3   4       5         6    7   8   9   10  11  12    13      14  15  16  17  18  19
		//Order of global_params: v_r,lambda_1,lambda_2,v_h,A_r,RC
		//The numbering is:        0      1        2     3   4   5
		//Need to set entries 12-19 of params.
		double* vals = params->ve;
		double K_T = 1.0;
		double s_r = 1.0;
		double rootS_h = pow(vals[7],.5);
		double L = params->ve[0];
		double A_h = params->ve[1] * 1e6;	//Put into m^2
		double eta = params->ve[8];
		double v_r = global_params->ve[0];
		double lambda_1 = global_params->ve[1];
		double lambda_2 = global_params->ve[2];
		double v_h = global_params->ve[3];
		double A_r = global_params->ve[4];
		double RC = global_params->ve[5];


		vals[12] = 60.0*v_r*pow(vals[2]/A_r,lambda_2)/((1.0-lambda_1)*vals[0]);	//invtau [1/min]
		vals[13] = vals[3] / s_r; //epsilon
		vals[14] = v_h*L;	//c_1 [m^2/s]
		vals[15] = vals[6] * vals[0] * vals[3] / 3600.0; //c_2
		vals[16] = (1e-3/60.0) * RC;	//c_3
		vals[17] = 60.0*v_h*L/A_h;	//c_4 [1/min], A_h converted above
		vals[18] = K_T/60.0;
		vals[19] = vals[6]/(60.0*s_r);


		iparams->ve[0] = link_i->location;
	}
	else if(type == 2000)
	{
		//Order of parameters: L_i,A_h,A_i,h_b,h_H,max_inf_rate,K_sat,S_h,eta,b_H,c_H,d_H,invtau,epsilon,c_1,c_2,c_3,c_4,c_5,c_6
		//The numbering is:     0   1   2   3   4       5         6    7   8   9   10  11  12    13      14  15  16  17  18  19
		//Order of global_params: v_r,lambda_1,lambda_2,Q_r,A_r,RC
		//The numbering is:        0      1        2     3   4   5
		//Need to set entries 12-19 of params.
		double* vals = params->ve;
		double K_T = 1.0;
		double s_r = 1.0;
		double rootS_h = pow(vals[7],.5);
		double v_r = global_params->ve[0];
		double lambda_1 = global_params->ve[1];
		double lambda_2 = global_params->ve[2];
		double Q_r = global_params->ve[3];
		double A_r = global_params->ve[4];
		double RC = global_params->ve[5];

		vals[12] = 60.0*v_r*pow(vals[2]/A_r,lambda_2)/((1.0-lambda_1)*vals[0]);
		vals[13] = vals[3] / s_r;
		vals[14] = 2.0 / .6 * vals[0] * s_r / Q_r * 1.0 / vals[8] * rootS_h * pow(vals[3],2.0/3.0);
		vals[15] = vals[6] * vals[0] * vals[3] / (3600.0 * Q_r);
		vals[16] = 1e-3/(60*s_r) * RC;
		vals[17] = 60 * 1e-6 * 2.0 / .6 * vals[0] / vals[1] * 1.0 / vals[8] * rootS_h;
		vals[18] = K_T/60.0;
		vals[19] = vals[6]/(60.0*s_r);
	}
}

//Set the initial condition for the differential equations a link. This method will be called once for each link. The differential
//	variables from any .ini, .uini, or database are already set here. Precalculations have also been made previously. All
//	algebraic variables MUST be set here. Other initial conditions may be set here as well, depending on the model.
//VEC* global_params:	The vector of global parameters for the system.
//VEC* params:		The parameters for this link. DEM, dam parameters, and precalculations are available.
//unsigned int dam:	1 if a dam is present at this link, 0 if no dam is present.
//VEC* y_0:		The initial condition vector. Store the initial data here.
//unsigned int type:	The index of the model.
//Returns the state of the solution (use 0 if state discontinuities are not a concern).
unsigned int ReadInitData(VEC* global_params,VEC* params,IVEC* iparams,QVSData* qvs,unsigned short int dam,VEC* y_0,unsigned int type)
{
	unsigned int state;

	if(type == 19)
	{
		//For type 21, just set the state
		return LinearHillslope_Evap_Check(y_0,params,global_params,qvs,0);
	}
	else if(type == 21)
	{
		//For type 21, only the storage (S or y_0->ve[1]) has been set.
		//Order of parameters: A_i,L_i,A_h,k2,k3,invtau,orifice_area,H_spill,H_max,S_max,alpha,orifice_diam,c_1,c_2,L_spill
		//The numbering is:	0   1   2  3  4    5	       6      7       8     9	  10	    11       12  13  14
		//Order of global_params: v_r,lambda_1,lambda_2,RC,S_0,v_h
		//The numbering is:        0      1        2     3  4   5
		double RC = global_params->ve[3];
		double S_0 = global_params->ve[4];
		double A_h = params->ve[2];
		y_0->ve[2] = RC * S_0 * A_h;
		y_0->ve[3] = (1.0 - RC) * S_0 * A_h;

		state = dam_check(y_0,global_params,params,qvs,dam); 
		dam_q(y_0,global_params,params,qvs,state,y_0);
		return state;
	}
	else if(type == 22)
	{
		//For type 22, only the storage (S or y_0->ve[1]) has been set.
		//Order of parameters: A_i,L_i,A_h,RC,v_h,v_r,k2,k3,invtau,orifice_area,H_spill,H_max,S_max,alpha,orifice_diam,c_1,c_2,L_spill
		//The numbering is:	0   1   2  3   4   5   6  7   8          9	  10	  11   12     13      14        15  16   17
		//Order of global_params: lambda_1,lambda_2,S_0
		//The numbering is:         0        1       2
		double S_0 = global_params->ve[2];
		double A_h = params->ve[2];
		double RC = params->ve[3];

		y_0->ve[2] = RC * S_0 * A_h;
		y_0->ve[3] = (1.0 - RC) * S_0 * A_h;

		state = dam_check2(y_0,global_params,params,qvs,dam); 
		dam_q2(y_0,global_params,params,qvs,state,y_0);
		return state;
	}
	else if(type == 23)
	{
		//For type 22, only the storage (S or y_0->ve[1]) has been set.
		//Order of parameters: A_i,L_i,A_h,RC,v_h,v_r,k2,k3,invtau,orifice_area,H_spill,H_max,S_max,alpha,orifice_diam,c_1,c_2,L_spill
		//The numbering is:	0   1   2  3   4   5   6  7   8          9	  10	  11   12     13      14        15  16   17
		//Order of global_params: lambda_1,lambda_2,S_0
		//The numbering is:         0        1       2
		double S_0 = global_params->ve[2];
		double A_h = params->ve[2];
		double RC = params->ve[3];

		y_0->ve[2] = RC * S_0 * A_h;
		y_0->ve[3] = (1.0 - RC) * S_0 * A_h;

		state = dam_check3(y_0,global_params,params,qvs,dam); 
		dam_q3(y_0,global_params,params,qvs,state,y_0);
		return state;
	}
	else if(type == 40)
	{
		//For type 22, only the storage (S or y_0->ve[1]) has been set.
		//Order of parameters: A_i,L_i,A_h,RC,v_h,v_r,k2,k3,invtau
		//The numbering is:	0   1   2  3   4   5   6  7   8
		//Order of global_params: lambda_1,lambda_2,S_0
		//The numbering is:         0        1       2
		double S_0 = global_params->ve[2];
		double A_h = params->ve[2];
		double RC = params->ve[3];

		y_0->ve[2] = RC * S_0 * A_h;
		y_0->ve[3] = (1.0 - RC) * S_0 * A_h;

		state = dam_check_qvs(y_0,global_params,params,qvs,dam);
		dam_q_qvs(y_0,global_params,params,qvs,state,y_0);
		return state;
	}
	else if(type == 30)
	{
		//y_i has q (y_i[0]), s_p (y_i[1]), h_w (y_i[2]), theta (y_i[3])
		//Order of parameters: L_h,A_h,A_up,H_b,H_h,max_inf_rate,K_SAT,S_H,n_vh, b, c, d,K_Q,V_T,c_1,c_2,c_3,c_4,c_5,c_6,c_7
		//The numbering is:     0   1   2    3   4       5         6    7   8    9 10 11  12  13  14  15  16  17  18  19  20
		//Order of global_params: v_0,lambda_1,lambda_2,Q_r,A_r,K_T,C_r,e_pot
		//The numbering is:        0      1        2     3   4   5   6   7
/*
		double H_h = params->ve[4];

params->ve[9] = 1.0;
params->ve[10] = 0.0;
params->ve[11] = 0.0;

		if(H_h < 1e-12)		//Flat surface
		{
			y_0->ve[2] = 0.0;
			y_0->ve[3] = 0.0;
		}
*/
		return 0;
	}
	else if(type == 200)
	{
		//For type 200, only the discharge has been set. Need to set the storage.
		//Order of parameters: L_i,A_h,A_i,h_b,h_H,max_inf_rate,K_sat,S_h,eta,b_H,c_H,d_H,invtau,epsilon,c_1,c_2,c_3,c_4,c_5,c_6
		//The numbering is:     0   1   2   3   4       5         6    7   8   9   10  11  12    13      14  15  16  17  18  19
		//Order of global_params: v_0,lambda_1,lambda_2,Q_r,A_r,RC,u_0
		//The numbering is:        0      1        2     3   4   5  6
		y_0->ve[1] = params->ve[0] / (global_params->ve[6] + global_params->ve[0]) * y_0->ve[0];
		return 0;
	}
	else if(type == 300)
	{
		//For this type, all initial conditions for variational equation must be set here.
		//Order of parameters: L_i,A_h,A_i,h_b,h_H,max_inf_rate,K_sat,S_h,eta,b_H,c_H,d_H,invtau,epsilon,c_1,c_2,c_3,c_4,c_5,c_6
		//The numbering is:     0   1   2   3   4       5         6    7   8   9   10  11  12    13      14  15  16  17  18  19
		//Order of global_params: v_r,lambda_1,lambda_2,Q_r,A_r,RC
		//The numbering is:        0      1        2     3   4   5
		unsigned int i;
		unsigned int offset = type - 299;
		for(i=offset;i<y_0->dim;i++)	y_0->ve[i] = 0.0;
		y_0->ve[iparams->ve[0] + offset] = 1.0;
		return 0;
	}
	else if(type == 301)
	{
		//For this type, all initial conditions for variational equation must be set here.
		//Order of parameters: L_i,A_h,A_i,h_b,h_H,max_inf_rate,K_sat,S_h,eta,b_H,c_H,d_H,invtau,epsilon,c_1,c_2,c_3,c_4,c_5,c_6
		//The numbering is:     0   1   2   3   4       5         6    7   8   9   10  11  12    13      14  15  16  17  18  19
		//Order of global_params: v_r,lambda_1,lambda_2,Q_r,A_r,RC
		//The numbering is:        0      1        2     3   4   5
		unsigned int i;
		unsigned int offset = type - 299;

		//New
		y_0->ve[offset] = 1.0;
		y_0->ve[offset + 1] = 1.0;
		y_0->ve[offset + 2] = 0.0;
		for(i=offset+3;i<y_0->dim;i++)	y_0->ve[i] = 0.0;

/*
		//Old
		unsigned int all_states = (y_0->dim - offset)/2;
		for(i=offset;i<y_0->dim;i++)	y_0->ve[i] = 0.0;
		y_0->ve[iparams->ve[0]*2 + offset] = 1.0;
		y_0->ve[iparams->ve[0]*2 + 1 + all_states + offset] = 1.0;
*/

		return 0;
	}
	else if(type == 315)
	{
		//For this type, all initial conditions for variational equation must be set here.
		//Order of parameters: L_i,A_h,A_i,h_b,h_H,max_inf_rate,K_sat,S_h,eta,b_H,c_H,d_H,invtau,epsilon,c_1,c_2,c_3,c_4,c_5,c_6
		//The numbering is:     0   1   2   3   4       5         6    7   8   9   10  11  12    13      14  15  16  17  18  19
		//Order of global_params: v_r,lambda_1,lambda_2,Q_r,A_r,RC
		//The numbering is:        0      1        2     3   4   5
		unsigned int i;
		unsigned int offset = 2;

		//New
		y_0->ve[offset] = 1.0;
		y_0->ve[offset + 1] = 1.0;
		y_0->ve[offset + 2] = 0.0;
		for(i=offset+3;i<y_0->dim;i++)	y_0->ve[i] = 0.0;

		return 0;
	}
	else
	{
		//If not using algebraic variables, then everything is already set
		return 0;
	}
}

/*
//!!!! Note: This routine will be phased out in the future. !!!!
//Sets the dimensions for the type of .ini or .rkd file. These should be set as dim in SetParamsSizes.
//unsigned int rktype:		The type of the .rkd file.
//unsigned int inittype:	The type of the .ini file.
//unsigned int* rkdim:		The number of unknowns the .rkd file assumes.
//unsigned int* initdim:	The number of unknowns the .ini file assumes.
void SetRKInitDim(unsigned int rktype,unsigned int inittype,unsigned int* rkdim,unsigned int* initdim)
{
	if(rktype == 6)				*rkdim = 4;
	else if(rktype == 4 || rktype == 5)	*rkdim = 4;
	else if(rktype == 1 || rktype == 2)	*rkdim = 2;
	else if(rktype == 3)			*rkdim = 3;
	else if(rktype == 0)			*rkdim = 2;
	else if(rktype == 101)			*rkdim = 3;
	else if(rktype == 15)			*rkdim = 1;
	else if(rktype == 19)			*rkdim = 3;
	else if(rktype == 20)			*rkdim = 3;
	else if(rktype == 21)			*rkdim = 3;
	else if(rktype == 22 || rktype == 23)	*rkdim = 3;
	else if(rktype == 30)			*rkdim = 4;
	else if(rktype == 40)			*rkdim = 3;
	else if(rktype == 60)			*rkdim = 3;
	else if(rktype == 105)			*rkdim = 2;
	else if(rktype == 190)			*rkdim = 3;
	else if(rktype == 200)			*rkdim = 2;
	else if(rktype == 219)			*rkdim = 1;
	else if(rktype == 250)			*rkdim = 3;
	else if(rktype == 252)			*rkdim = 4;
	else if(rktype == 2000)			*rkdim = 2;
	else if(rktype == 300 || rktype == 301 || rktype == 315)
	{
		printf("Error: Type %u does not work with rkd files.\n",rktype);
		abort();
	}
	else					*rkdim = -1;

	if(inittype == 6)			*initdim = 4;
	else if(inittype == 4 || inittype == 5)	*initdim = 4;
	else if(inittype == 1 || inittype == 2)	*initdim = 2;
	else if(inittype == 3)			*initdim = 3;
	else if(inittype == 0)			*initdim = 1;
	else if(inittype == 101)		*initdim = 3;
	else if(inittype == 15)			*initdim = 2;
	else if(inittype == 19)			*initdim = 3;
	else if(inittype == 20)			*initdim = 3;
	else if(inittype == 21)			*initdim = 3;	//!!!! Should be 1? !!!!
	else if(inittype == 22 || inittype == 23)*initdim = 3;	//!!!! Should be 1? !!!!
	else if(inittype == 30)			*initdim = 4;
	else if(inittype == 40)			*initdim = 3;	//!!!! Should be 1? !!!!
	else if(inittype == 60)			*rkdim = 3;
	else if(inittype == 105)		*initdim = 2;
	else if(inittype == 190)		*initdim = 3;
	else if(inittype == 200)		*initdim = 1;
	else if(inittype == 219)		*initdim = 1;
	else if(inittype == 250)		*initdim = 3;
	else if(inittype == 252)		*initdim = 4;
	else if(inittype == 300)		*initdim = 1;
	else if(inittype == 301)		*initdim = 2;
	else if(inittype == 315)		*initdim = 2;
	else if(inittype == 2000)		*initdim = 2;
	else 					*initdim = -1;
}
*/

//If using data assimilation, sets the global errors and dim correctly
void AssimError(unsigned int N,UnivVars* GlobalVars,ErrorData* GlobalErrors)
{
	unsigned int i,j;
	unsigned int old_num_dense;

	//Remove any variaional indices from dense_indices
	for(i=0;i<GlobalVars->num_dense;i++)
	{
		if(GlobalVars->dense_indices[i] >= GlobalVars->problem_dim)
		{
			for(j=i+1;j<GlobalVars->num_dense;j++)	GlobalVars->dense_indices[j-1] = GlobalVars->dense_indices[j];
			(GlobalVars->num_dense)--;
			i--;
		}
	}

	old_num_dense = GlobalVars->num_dense;

	GlobalVars->dim = GlobalVars->problem_dim + N*GlobalVars->problem_dim*GlobalVars->problem_dim;
	GlobalVars->num_dense += N*GlobalVars->problem_dim*GlobalVars->problem_dim;
	GlobalErrors->abstol->ve = realloc(GlobalErrors->abstol->ve,GlobalVars->dim*sizeof(double));
	GlobalErrors->reltol->ve = realloc(GlobalErrors->reltol->ve,GlobalVars->dim*sizeof(double));
	GlobalErrors->abstol_dense->ve = realloc(GlobalErrors->abstol_dense->ve,GlobalVars->dim*sizeof(double));
	GlobalErrors->reltol_dense->ve = realloc(GlobalErrors->reltol_dense->ve,GlobalVars->dim*sizeof(double));
	GlobalErrors->abstol->dim = GlobalErrors->reltol->dim = GlobalErrors->reltol_dense->dim = GlobalErrors->reltol_dense->dim = GlobalVars->dim;

	//Setup error
	for(i=GlobalVars->problem_dim+1;i<GlobalVars->dim;i++)
	{
		GlobalErrors->abstol->ve[i] = GlobalErrors->abstol->ve[GlobalVars->problem_dim];
		GlobalErrors->reltol->ve[i] = GlobalErrors->reltol->ve[GlobalVars->problem_dim];
		GlobalErrors->abstol_dense->ve[i] = GlobalErrors->abstol_dense->ve[GlobalVars->problem_dim];
		GlobalErrors->reltol_dense->ve[i] = GlobalErrors->reltol_dense->ve[GlobalVars->problem_dim];
	}

	//Setup dense indices
	//Add in variational indices
	GlobalVars->dense_indices = realloc(GlobalVars->dense_indices,GlobalVars->num_dense * sizeof(unsigned int));
	for(i=old_num_dense;i<GlobalVars->num_dense;i++)
	{
		GlobalVars->dense_indices[i] = i;
	}
}

/*
//!!!! Consider putting this into problems and setting a function pointer !!!!
//Checks that the variable y is consistent with constraints
//VEC* y:		The current state of the system of ODEs.
//unsigned int dim:	The dimension of the ODEs.
//unsigned int type:	The index of the model.
void CheckConsistency(VEC* y,unsigned int dim,unsigned int type,VEC* params,VEC* global_params)
{
	if(type <= 5)	//!!!! Break this up, remove dim !!!!
	{
		if(y->ve[0] < 1e-14)		y->ve[0] = 1e-14;
		if(dim > 1 && y->ve[1] < 0.0)	y->ve[1] = 0.0;
		if(dim > 2 && y->ve[2] < 0.0)	y->ve[2] = 0.0;
		if(dim > 2 && y->ve[2] > 1.0)	y->ve[2] = 1.0;
		if(dim > 3 && y->ve[3] < 0.0)	y->ve[3] = 0.0;
		if(dim > 3 && y->ve[3] > 1.0 - y->ve[2])	y->ve[3] = 1.0 - y->ve[2];
	}
	else if(type == 6)
	{
		if(y->ve[0] < 1e-14)	y->ve[0] = 1e-14;
		if(y->ve[1] < 0.0)	y->ve[1] = 0.0;
		if(y->ve[2] < 0.0)	y->ve[2] = 0.0;
		if(y->ve[3] < 0.0)	y->ve[3] = 0.0;
	}
	else if(type == 15)
	{
		if(y->ve[0] < 1e-12)	y->ve[0] = 1e-12;
		if(y->ve[1] < 0.0)	y->ve[1] = 0.0;
	}
	else if(type == 19 || type == 60 || type == 190)
	{
		if(y->ve[0] < 1e-14)	y->ve[0] = 1e-14;
		if(y->ve[1] < 0.0)	y->ve[1] = 0.0;
		if(y->ve[2] < 0.0)	y->ve[2] = 0.0;
	}
	else if(type == 20)
	{
		if(y->ve[0] < 1e-14)	y->ve[0] = 1e-14;
		if(y->ve[1] < 0.0)	y->ve[1] = 0.0;
		if(y->ve[2] < 0.0)	y->ve[2] = 0.0;
	}
	else if(type == 21 || type == 22 || type == 23 || type == 40)
	{
		if(y->ve[0] < 1e-14)	y->ve[0] = 1e-14;
		if(y->ve[1] < 0.0)	y->ve[1] = 0.0;
		if(y->ve[2] < 0.0)	y->ve[2] = 0.0;
		if(y->ve[3] < 0.0)	y->ve[3] = 0.0;
	}
	else if(type == 30)
	{
		if(y->ve[0] < 1e-14)	y->ve[0] = 1e-14;
		if(y->ve[1] < 0.0)	y->ve[1] = 0.0;
		if(y->ve[2] < 0.0)	y->ve[2] = 0.0;
		if(y->ve[2] > params->ve[4])	y->ve[2] = params->ve[4];
		if(y->ve[3] < 0.0)	y->ve[3] = 0.0;
		else if(y->ve[3] > 1.0)	y->ve[3] = 1.0;
	}
	else if(type == 105)
	{
		if(y->ve[0] < 0.0)	y->ve[0] = 0.0;
		if(y->ve[1] < 0.0)	y->ve[1] = 0.0;
	}
	else if(type == 200)
	{
		if(y->ve[0] < 0.0)	y->ve[0] = 0.0;
		if(y->ve[1] < 0.0)	y->ve[1] = 0.0;
	}
	else if(type == 219)
	{
		if(y->ve[0] < 1e-14)	y->ve[0] = 1e-14;
	}
	else if(type == 250)
	{
		if(y->ve[0] < 1e-14)	y->ve[0] = 1e-14;
		if(y->ve[1] < 0.0)	y->ve[1] = 0.0;
		if(y->ve[2] < 0.0)	y->ve[2] = 0.0;
		if(y->ve[2] > global_params->ve[7])	y->ve[2] = global_params->ve[7];
	}
	else if(type == 252)
	{
		if(y->ve[0] < 1e-14)	y->ve[0] = 1e-14;
		if(y->ve[1] < 0.0)	y->ve[1] = 0.0;
		if(y->ve[2] < 0.0)	y->ve[2] = 0.0;
		if(y->ve[3] < 0.0)	y->ve[3] = 0.0;
	}
	else if(type == 300)
	{
		if(y->ve[0] < 0.0)	y->ve[0] = 1e-12;
	}
	else if(type == 301)
	{
		if(y->ve[0] < 0.0)	y->ve[0] = 1e-12;
		if(y->ve[1] < 0.0)	y->ve[1] = 0.0;
	}
	else if(type == 315)
	{
		if(y->ve[0] < 1e-12)	y->ve[0] = 1e-12;
		if(y->ve[1] < 0.0)	y->ve[1] = 0.0;
	}
	else if(type == 2000)
	{
		if(y->ve[0] < 1e-14)		y->ve[0] = 1e-14;
		if(y->ve[1] < 0.0)		y->ve[1] = 0.0;
	}
}
*/


