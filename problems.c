#include "problems.h"

extern VEC* dump;

//extern int flaggy;

//These are the right-hand side functions for the ODEs. For each of these functions:
//double t: The current time
//VEC* y_i: The approximate value of the solution to the ODEs at the current link at time t
//VEC** y_p: y_p[j] has the approximate values for the immediately upstream link j to the current link at time t
//unsigned short int numparents: The number of upstream links (parents) to link i
//VEC* global_params: The global parameters
//RainData* rain: The rain fall values for link i
//VEC* params: The parameters for link i
//int state: The current state of the system
//VEC* ans (set by method, assumed that space is allocated): The value returned by the right-hand side function

/*
//Type 2000
//Test for equations from a parser. Same as Type 1.
//Calculates the flow using simple parameters, using only the flow q.
//Order of parameters: L_i,A_h,A_i,h_b,h_H,max_inf_rate,K_sat,S_h,eta,b_H,c_H,d_H,invtau,epsilon,c_1,c_2,c_3,c_4,c_5,c_6
//The numbering is:     0   1   2   3   4       5         6    7   8   9   10  11  12    13      14  15  16  17  18  19
//Order of global_params: v_r,lambda_1,lambda_2,Q_r,A_r,RC
//The numbering is:        0      1        2     3   4   5
//This uses the units and functions from September 18, 2011 document
//y_i[0] = q, y_i[1] = s
void parser_test(double t,VEC* y_i,VEC** y_p,unsigned short int numparents,VEC* global_params,double* forcing_values,QVSData* qvs,VEC* params,IVEC* iparams,int state,unsigned int** upstream,unsigned int* numupstream,VEC* ans)
{
	unsigned int i;
	unsigned int dim = y_i->dim;
	double* tempresult;
	double inflow = 0.0;

	Formula* equations = (Formula*) qvs;	//!!!! Pass in somehow !!!!
	muParserHandle_t* parser = equations->parser;	
	VEC* variable_values = equations->variable_values;
	//int numvars = mupGetVarNum(parser);


	//Set variables
	for(i=0;i<dim;i++)	variable_values->ve[i] = y_i->ve[i];	//States
	variable_values->ve[dim] = t;					//Time
	for(i=0;i<numparents;i++)	inflow += y_p[i]->ve[0];
	variable_values->ve[dim+1] = inflow;				//Inflow
	variable_values->ve[dim+2] = forcing_values[0];			//Rainfall

	//Evaluate equations
	tempresult = mupEvalMulti(parser,&i);
	for(i=0;i<dim;i++)	ans->ve[i] = tempresult[i];
}
*/



//Type 250
//Order of parameters: A_i,L_i,A_h,h_r,invtau,k_2,k_I,c_1,c_2
//The numbering is:	0   1   2   3    4     5   6   7   8
//Order of global_params: v_0,lambda_1,lambda_2,v_h,k_3,k_I_factor,gamma,h_b,e_pot
//The numbering is:        0      1        2     3   4     5         6    7	8
void NonLinearHillslope(double t,VEC* y_i,VEC** y_p,unsigned short int numparents,VEC* global_params,double* forcing_values,QVSData* qvs,VEC* params,IVEC* iparams,int state,unsigned int** upstream,unsigned int* numupstream,VEC* ans)
{
	unsigned short int i;

	double lambda_1 = global_params->ve[1];
	double k_3 = global_params->ve[4];	//[1/min]
	double gamma = global_params->ve[6];
	double h_b = global_params->ve[7];	//[m]
	//double e_pot = global_params->ve[8] * (1e-3/60.0);	//[mm/hr]->[m/min]
	double e_pot = global_params->ve[8];	//[m/min]


	double L = params->ve[1];	//[m]
	double A_h = params->ve[2];	//[m^2]
	double h_r = params->ve[3];	//[m]
	double invtau = params->ve[4];	//[1/min]
	double k_2 = params->ve[5];	//[1/min]
	double k_I = params->ve[6];	//[1/min]
	double c_1 = params->ve[7];
	double c_2 = params->ve[8];

	double q = y_i->ve[0];		//[m^3/s]
	double s_p = y_i->ve[1];	//[m]
	double s_s = y_i->ve[2];	//[m]

	//Evaporation
	double C_p,C_s,C_T;
	if(e_pot > 0.0)
	{
		C_p = s_p / e_pot;
		C_s = s_s / e_pot;
		C_T = C_p + C_s;
	}
	else
	{
		C_p = 0.0;
		C_s = 0.0;
		C_T = 0.0;
	}

	double Corr_evap = (C_T > 1.0) ? 1.0/C_T : 1.0;
	double e_p = Corr_evap * C_p * e_pot;
	double e_s = Corr_evap * C_s * e_pot;

	//Fluxes
	double q_pl = k_2 * s_p;
	double q_ps = k_I * s_p * pow(1.0 - s_s/h_b,gamma);
	double q_sl = k_3 * sq(s_s * L / A_h) * (h_r + h_b);

	//Discharge
	ans->ve[0] = -q + (q_pl + q_sl) * c_2;
	for(i=0;i<numparents;i++)
		ans->ve[0] += y_p[i]->ve[0];
	ans->ve[0] = invtau * pow(q,lambda_1) * ans->ve[0];

	//Hillslope
	ans->ve[1] = forcing_values[0]*c_1 - q_pl - q_ps - e_p;
	ans->ve[2] = q_ps - q_sl - e_s;
}


//Type 252
//Contains 3 layers on hillslope: ponded, top layer, soil
//Order of parameters: A_i,L_i,A_h,invtau,k_2,k_i,c_1,c_2
//The numbering is:	0   1   2     3    4   5   6   7
//Order of global_params: v_0,lambda_1,lambda_2,v_h,k_3,k_I_factor,h_b,S_L,A,B,exponent
//The numbering is:        0      1        2     3   4     5        6   7  8 9  10
void TopLayerHillslope(double t,VEC* y_i,VEC** y_p,unsigned short int numparents,VEC* global_params,double* forcing_values,QVSData* qvs,VEC* params,IVEC* iparams,int state,unsigned int** upstream,unsigned int* numupstream,VEC* ans)
{
	unsigned short int i;

	double lambda_1 = global_params->ve[1];
	double k_3 = global_params->ve[4];	//[1/min]
	double h_b = global_params->ve[6];	//[m]
	double S_L = global_params->ve[7];	//[m]
	double A = global_params->ve[8];
	double B = global_params->ve[9];
	double exponent = global_params->ve[10];
	double e_pot = forcing_values[1] * (1e-3/(30.0*24.0*60.0));	//[mm/month] -> [m/min]
	//double e_pot = global_params->ve[11];	//[m/min]
	//double e_pot = global_params->ve[11] * (1e-3*60.0);	//[m/min]
	//double e_pot = 0.0;

	double L = params->ve[1];	//[m]
	double A_h = params->ve[2];	//[m^2]
	//double h_r = params->ve[3];	//[m]
	double invtau = params->ve[3];	//[1/min]
	double k_2 = params->ve[4];	//[1/min]
	double k_i = params->ve[5];	//[1/min]
	double c_1 = params->ve[6];
	double c_2 = params->ve[7];

	double q = y_i->ve[0];		//[m^3/s]
	double s_p = y_i->ve[1];	//[m]
	double s_t = y_i->ve[2];	//[m]
	double s_s = y_i->ve[3];	//[m]

	//Evaporation
	double e_p,e_t,e_s;
	double Corr = s_p + s_t / S_L + s_s / (h_b-S_L);
	if(e_pot > 0.0 && Corr > 1e-12)
	{
		e_p = s_p * e_pot / Corr;
		e_t = s_t/S_L * e_pot / Corr;
		e_s = s_s/(h_b-S_L) * e_pot / Corr;
	}
	else
	{
		e_p = 0.0;
		e_t = 0.0;
		e_s = 0.0;
	}

	double pow_term = (1.0-s_t/S_L > 0.0) ? pow(1.0-s_t/S_L,exponent) : 0.0;
	double k_t = (A + B * pow_term) * k_2;

	//Fluxes
	double q_pl = k_2 * s_p;
	double q_pt = k_t * s_p;
	double q_ts = k_i * s_t;
	double q_sl = k_3 * s_s;

	//Discharge
	ans->ve[0] = -q + (q_pl + q_sl) * c_2;
	for(i=0;i<numparents;i++)
		ans->ve[0] += y_p[i]->ve[0];
	ans->ve[0] = invtau * pow(q,lambda_1) * ans->ve[0];

	//Hillslope
	ans->ve[1] = forcing_values[0]*c_1 - q_pl - q_pt - e_p;
	ans->ve[2] = q_pt - q_ts - e_t;
	ans->ve[3] = q_ts - q_sl - e_s;
}


//Type 253
//Contains 3 layers on hillslope: ponded, top layer, soil
//Order of parameters: A_i,L_i,A_h,invtau,k_2,k_i,c_1,c_2
//The numbering is:	0   1   2     3    4   5   6   7
//Order of global_params: v_0,lambda_1,lambda_2,v_h,k_3,k_I_factor,h_b,S_L,A,B,exponent
//The numbering is:        0      1        2     3   4     5        6   7  8 9  10
void TopLayerHillslope_Reservoirs(double t,VEC* y_i,VEC** y_p,unsigned short int numparents,VEC* global_params,double* forcing_values,QVSData* qvs,VEC* params,IVEC* iparams,int state,unsigned int** upstream,unsigned int* numupstream,VEC* ans)
{
	ans->ve[0] = forcing_values[2];
	ans->ve[1] = 0.0;
	ans->ve[2] = 0.0;
	ans->ve[3] = 0.0;
}


//Type 254
//Contains 3 layers on hillslope: ponded, top layer, soil. Also has 3 extra states: total precip, total runoff, base flow
//Order of parameters: A_i,L_i,A_h,invtau,k_2,k_i,c_1,c_2
//The numbering is:	0   1   2     3    4   5   6   7
//Order of global_params: v_0,lambda_1,lambda_2,v_h,k_3,k_I_factor,h_b,S_L,A,B,exponent,v_B
//The numbering is:        0      1        2     3   4     5        6   7  8 9  10       11
void TopLayerHillslope_extras(double t,VEC* y_i,VEC** y_p,unsigned short int numparents,VEC* global_params,double* forcing_values,QVSData* qvs,VEC* params,IVEC* iparams,int state,unsigned int** upstream,unsigned int* numupstream,VEC* ans)
{
	unsigned short int i;

	double lambda_1 = global_params->ve[1];
	double k_3 = global_params->ve[4];	//[1/min]
	double h_b = global_params->ve[6];	//[m]
	double S_L = global_params->ve[7];	//[m]
	double A = global_params->ve[8];
	double B = global_params->ve[9];
	double exponent = global_params->ve[10];
	double v_B = global_params->ve[11];
	double e_pot = forcing_values[1] * (1e-3/(30.0*24.0*60.0));	//[mm/month] -> [m/min]

	double L = params->ve[1];	//[m]
	double A_h = params->ve[2];	//[m^2]
	//double h_r = params->ve[3];	//[m]
	double invtau = params->ve[3];	//[1/min]
	double k_2 = params->ve[4];	//[1/min]
	double k_i = params->ve[5];	//[1/min]
	double c_1 = params->ve[6];
	double c_2 = params->ve[7];

	double q = y_i->ve[0];		//[m^3/s]
	double s_p = y_i->ve[1];	//[m]
	double s_t = y_i->ve[2];	//[m]
	double s_s = y_i->ve[3];	//[m]
	//double s_precip = y_i->ve[4];	//[m]
	//double V_r = y_i->ve[5];	//[m^3]
	double q_b = y_i->ve[6];	//[m^3/s]

	//Evaporation
	double e_p,e_t,e_s;
	double Corr = s_p + s_t / S_L + s_s / (h_b-S_L);
	if(e_pot > 0.0 && Corr > 1e-12)
	{
		e_p = s_p * e_pot / Corr;
		e_t = s_t/S_L * e_pot / Corr;
		e_s = s_s/(h_b-S_L) * e_pot / Corr;
	}
	else
	{
		e_p = 0.0;
		e_t = 0.0;
		e_s = 0.0;
	}

	double pow_term = (1.0-s_t/S_L > 0.0) ? pow(1.0-s_t/S_L,exponent) : 0.0;
	double k_t = (A + B * pow_term) * k_2;

	//Fluxes
	double q_pl = k_2 * s_p;
	double q_pt = k_t * s_p;
	double q_ts = k_i * s_t;
	double q_sl = k_3 * s_s;	//[m/min]

	//Discharge
	ans->ve[0] = -q + (q_pl + q_sl) * c_2;
	for(i=0;i<numparents;i++)
		ans->ve[0] += y_p[i]->ve[0];
	ans->ve[0] = invtau * pow(q,lambda_1) * ans->ve[0];

	//Hillslope
	ans->ve[1] = forcing_values[0]*c_1 - q_pl - q_pt - e_p;
	ans->ve[2] = q_pt - q_ts - e_t;
	ans->ve[3] = q_ts - q_sl - e_s;

	//Additional states
	ans->ve[4] = forcing_values[0]*c_1;
	ans->ve[5] = q_pl;
	ans->ve[6] = q_sl * A_h - q_b*60.0;
	for(i=0;i<numparents;i++)
		ans->ve[6] += y_p[i]->ve[6] * 60.0;
		//ans->ve[6] += k_3*y_p[i]->ve[3]*A_h;
	ans->ve[6] *= v_B/L;
}


//Type 260
//Contains 3 layers on hillslope: ponded, top layer, soil
//Order of parameters: A_i,L_i,A_h | invtau,c_1,c_2
//The numbering is:	0   1   2  |    3    4   5 
//Order of global_params: v_0,lambda_1,lambda_2,h_b,k_D,k_2,k_dry,k_i,T_L,N,phi
//The numbering is:        0      1        2     3   4   5   6     7   8  9  10
void TopLayerNonlinearExp(double t,VEC* y_i,VEC** y_p,unsigned short int numparents,VEC* global_params,double* forcing_values,QVSData* qvs,VEC* params,IVEC* iparams,int state,unsigned int** upstream,unsigned int* numupstream,VEC* ans)
{
	unsigned short int i;
	double pers_to_permin = 60.0;

	//Global params
	double lambda_1 = global_params->ve[1];
	double h_b = global_params->ve[3];	//[m]
	double k_D = global_params->ve[4];	//[1/s]
	double k_2 = global_params->ve[5];	//[1/s]
	double k_dry = global_params->ve[6];	//[1/s]
	double k_i = global_params->ve[7];	//[1/s]
	double T_L = global_params->ve[8];	//[m]
	double N = global_params->ve[9];	//[]
	double phi = global_params->ve[10];	//[]

	//Forcings
	double e_pot = forcing_values[1] * (1e-3/(30.0*24.0*60.0));	//[mm/month] -> [m/min]

	//Precalculations
	double invtau = params->ve[3];	//[1/min]
	double c_1 = params->ve[4];
	double c_2 = params->ve[5];

	//System states
	double q = y_i->ve[0];		//[m^3/s]
	double s_p = y_i->ve[1];	//[m]
	double s_t = y_i->ve[2];	//[m]
	double s_s = y_i->ve[3];	//[m]

	//Evaporation
	double e_p,e_t,e_s;
	double Corr = s_p + s_t / T_L + s_s / (h_b-T_L);
	if(e_pot > 0.0 && Corr > 1e-12)
	{
		e_p = s_p * e_pot / Corr;
		e_t = s_t/T_L * e_pot / Corr;
		e_s = s_s/(h_b-T_L) * e_pot / Corr;
	}
	else
	{
		e_p = 0.0;
		e_t = 0.0;
		e_s = 0.0;
	}

	//Fluxes
	double q_pl = k_2 * pers_to_permin * pow(s_p,phi);
	double q_pt = k_dry * pers_to_permin * ((1.0-s_t/T_L > 0.0) ? pow(1.0-s_t/T_L,N) : 0.0) * s_p;
	double q_ts = k_i * pers_to_permin * s_t;
	double q_sl = k_D * pers_to_permin * s_s;

	//Discharge
	ans->ve[0] = -q + (q_pl + q_sl) * c_2;
	for(i=0;i<numparents;i++)
		ans->ve[0] += y_p[i]->ve[0];
	ans->ve[0] = invtau * pow(q,lambda_1) * ans->ve[0];

	//Hillslope
	ans->ve[1] = forcing_values[0]*c_1 - q_pl - q_pt - e_p;
	ans->ve[2] = q_pt - q_ts - e_t;
	ans->ve[3] = q_ts - q_sl - e_s;
}


//Type 261
//Contains 2 states in the channel: dischage, storage
//Contains 3 layers on hillslope: ponded, top layer, soil
//Contains 3 extra states: 
//Order of parameters: A_i,L_i,A_h,S_h,T_L,h_b,k_D,k_dry,k_i | invtau,c_1,c_2,c_3
//The numbering is:	0   1   2   3   4   5   6   7     8  |   9     10  11  12
//Order of global_params: v_0,lambda_1,lambda_2,N,phi,v_B
//The numbering is:        0      1        2    3  4   5 
void TopLayerNonlinearExpSoilvel(double t,VEC* y_i,VEC** y_p,unsigned short int numparents,VEC* global_params,double* forcing_values,QVSData* qvs,VEC* params,IVEC* iparams,int state,unsigned int** upstream,unsigned int* numupstream,VEC* ans)
{
	unsigned short int i;
	double pers_to_permin = 60.0;

	//Global params
	double lambda_1 = global_params->ve[1];
	double N = global_params->ve[3];	//[]
	double phi = global_params->ve[4];	//[]
	double v_B = global_params->ve[5];	//[m/s]

	//Local params
	double L_i = params->ve[1];	//[m]
	double A_h = params->ve[2];	//[m^2]
	double S_h = params->ve[3];	//[]
	double T_L = params->ve[4];	//[m]
	double h_b = global_params->ve[5];	//[m]
	double k_D = global_params->ve[6];	//[1/s]
	double k_dry = global_params->ve[7];	//[1/s]
	double k_i = global_params->ve[8];	//[1/s]

	//Precalculations
	double invtau = params->ve[9];	//[1/min]
	double c_1 = params->ve[10];
	double c_2 = params->ve[11];
	double c_3 = params->ve[12];

	//Forcings
	double rainfall = forcing_values[0] * c_1;			//[mm/hr] -> [m/min]
	double e_pot = forcing_values[1] * (1e-3/(30.0*24.0*60.0));	//[mm/month] -> [m/min]
	double eta = forcing_values[2];					//[]

	//System states
	double q = y_i->ve[0];		//[m^3/s]
	double S = y_i->ve[1];		//[m^3]
	double s_p = y_i->ve[2];	//[m]
	double s_t = y_i->ve[3];	//[m]
	double s_s = y_i->ve[4];	//[m]
	//double s_precip = y_i->ve[5];	//[m]
	//double V_r = y_i->ve[6];	//[m^3]
	double q_b = y_i->ve[7];	//[m^3/s]

	//Evaporation
	double e_p,e_t,e_s;
	double Corr = s_p + s_t / T_L + s_s / (h_b-T_L);
	if(e_pot > 0.0 && Corr > 1e-12)
	{
		e_p = s_p * e_pot / Corr;
		e_t = s_t/T_L * e_pot / Corr;
		e_s = s_s/(h_b-T_L) * e_pot / Corr;
	}
	else
	{
		e_p = 0.0;
		e_t = 0.0;
		e_s = 0.0;
	}

	//Fluxes
	double k_2 = c_3 / eta;
	double q_pl = k_2 * pers_to_permin * pow(s_p,phi);	//[m/min]
	double q_pt = k_dry * pers_to_permin * ((1.0-s_t/T_L > 0.0) ? pow(1.0-s_t/T_L,N) : 0.0) * s_p;	//[m/min]
	double q_ts = k_i * pers_to_permin * s_t;	//[m/min]
	double q_sl = k_D * pers_to_permin * s_s;	//[m/min]

	//Discharge
	dam_TopLayerNonlinearExpSoilvel(y_i,global_params,params,qvs,state,ans);	//ans is used for convenience
	double qm = ans->ve[0] * 60.0;

	//Storage
	ans->ve[1] = (q_pl + q_sl) * A_h - qm;
	for(i=0;i<numparents;i++)
		ans->ve[1] += y_p[i]->ve[0] * 60.0;

	//Hillslope
	ans->ve[2] = rainfall - q_pl - q_pt - e_p;
	ans->ve[3] = q_pt - q_ts - e_t;
	ans->ve[4] = q_ts - q_sl - e_s;

	//Additional states
	ans->ve[5] = rainfall;
	ans->ve[6] = q_pl;
	ans->ve[7] = q_sl * A_h - q_b*60.0;
	for(i=0;i<numparents;i++)
		ans->ve[7] += y_p[i]->ve[7] * 60.0;
	ans->ve[7] *= v_B/L_i;
}

//Type 261
//Contains 2 states in the channel: dischage, storage
//Contains 3 layers on hillslope: ponded, top layer, soil
//Order of parameters: A_i,L_i,A_h,S_h,T_L,h_b,k_D,k_dry,k_i | invtau,c_1,c_2,c_3
//The numbering is:	0   1   2   3   4   5   6   7     8  |   9     10  11  12
//Order of global_params: v_0,lambda_1,lambda_2,N,phi,v_B
//The numbering is:        0      1        2    3  4   5 
void dam_TopLayerNonlinearExpSoilvel(VEC* y,VEC* global_params,VEC* params,QVSData* qvs,int state,VEC* ans)
{
	double q1,q2,S1,S2,S_max,q_max,S;

	//Parameters
	double lambda_1 = global_params->ve[1];
	double invtau = params->ve[9];

	//Find the discharge in [m^3/s]
	if(state == -1)
	{
		S = (y->ve[1] < 0.0) ? 0.0 : y->ve[1];
		//ans->ve[0] = invtau/60.0*pow(S,1.0/(1.0-lambda_1));
		ans->ve[0] = pow((1.0-lambda_1)*invtau/60.0 * S,1.0/(1.0-lambda_1));
	}
	else if(state == (int) qvs->n_values - 1)
	{
		S_max = qvs->points[qvs->n_values - 1][0];
		q_max = qvs->points[qvs->n_values - 1][1];
		ans->ve[0] = q_max;
	}
	else
	{
		S = (y->ve[1] < 0.0) ? 0.0 : y->ve[1];
		q2 = qvs->points[state+1][1];
		q1 = qvs->points[state][1];
		S2 = qvs->points[state+1][0];
		S1 = qvs->points[state][0];
		ans->ve[0] = (q2-q1)/(S2-S1) * (S-S1) + q1;
	}
}

//Type 261
//Contains 2 states in the channel: dischage, storage
//Contains 3 layers on hillslope: ponded, top layer, soil
//Order of parameters: A_i,L_i,A_h,S_h,T_L,h_b,k_D,k_dry,k_i | invtau,c_1,c_2,c_3
//The numbering is:	0   1   2   3   4   5   6   7     8  |   9     10  11  12
//Order of global_params: v_0,lambda_1,lambda_2,N,phi,v_B
//The numbering is:        0      1        2    3  4   5 
void TopLayerNonlinearExpSoilvel_Reservoirs(double t,VEC* y_i,VEC** y_p,unsigned short int numparents,VEC* global_params,double* forcing_values,QVSData* qvs,VEC* params,IVEC* iparams,int state,unsigned int** upstream,unsigned int* numupstream,VEC* ans)
{
	ans->ve[0] = forcing_values[3];
	ans->ve[1] = 0.0;
	ans->ve[2] = 0.0;
	ans->ve[3] = 0.0;
	ans->ve[4] = 0.0;
	ans->ve[5] = 0.0;
	ans->ve[6] = 0.0;
	ans->ve[7] = 0.0;
}

//Type 262
//Contains 2 states in the channel: dischage, storage
//Contains 3 layers on hillslope: ponded, top layer, soil
//Order of parameters: A_i,L_i,A_h,S_h,T_L,eta,h_b,k_D,k_dry,k_i | invtau,c_1,c_2,k_2
//The numbering is:	0   1   2   3   4   5   6   7   8     9  |   10    11  12  13
//Order of global_params: v_0,lambda_1,lambda_2,N,phi,v_B
//The numbering is:        0      1        2    3  4   5
void TopLayerNonlinearExpSoilvel_ConstEta(double t,VEC* y_i,VEC** y_p,unsigned short int numparents,VEC* global_params,double* forcing_values,QVSData* qvs,VEC* params,IVEC* iparams,int state,unsigned int** upstream,unsigned int* numupstream,VEC* ans)
{
	unsigned short int i;
	double pers_to_permin = 60.0;

	//Global params
	double lambda_1 = global_params->ve[1];
	double N = global_params->ve[3];	//[]
	double phi = global_params->ve[4];	//[]
	double v_B = global_params->ve[5];	//[m/s]

	//Local params
	double L_i = params->ve[1];	//[m]
	double A_h = params->ve[2];	//[m^2]
	double S_h = params->ve[3];	//[]
	double T_L = params->ve[4];	//[m]
	double eta = params->ve[5];	//[]
	double h_b = global_params->ve[6];	//[m]
	double k_D = global_params->ve[7];	//[1/s]
	double k_dry = global_params->ve[8];	//[1/s]
	double k_i = global_params->ve[9];	//[1/s]

	//Precalculations
	double invtau = params->ve[10];	//[1/min]
	double c_1 = params->ve[11];
	double c_2 = params->ve[12];
	double k_2 = params->ve[13];

	//Forcings
	double rainfall = forcing_values[0] * c_1;			//[mm/hr] -> [m/min]
	double e_pot = forcing_values[1] * (1e-3/(30.0*24.0*60.0));	//[mm/month] -> [m/min]

	//System states
	double q = y_i->ve[0];		//[m^3/s]
	double S = y_i->ve[1];		//[m^3]
	double s_p = y_i->ve[2];	//[m]
	double s_t = y_i->ve[3];	//[m]
	double s_s = y_i->ve[4];	//[m]
	//double s_precip = y_i->ve[5];	//[m]
	//double V_r = y_i->ve[6];	//[m^3]
	double q_b = y_i->ve[7];	//[m^3/s]

	//Evaporation
	double e_p,e_t,e_s;
	double Corr = s_p + s_t / T_L + s_s / (h_b-T_L);
	if(e_pot > 0.0 && Corr > 1e-12)
	{
		e_p = s_p * e_pot / Corr;
		e_t = s_t/T_L * e_pot / Corr;
		e_s = s_s/(h_b-T_L) * e_pot / Corr;
	}
	else
	{
		e_p = 0.0;
		e_t = 0.0;
		e_s = 0.0;
	}

	//Fluxes
	//double k_2 = c_3 / eta;
	double q_pl = k_2 * pers_to_permin * pow(s_p,phi);	//[m/min]
	double q_pt = k_dry * pers_to_permin * ((1.0-s_t/T_L > 0.0) ? pow(1.0-s_t/T_L,N) : 0.0) * s_p;	//[m/min]
	double q_ts = k_i * pers_to_permin * s_t;	//[m/min]
	double q_sl = k_D * pers_to_permin * s_s;	//[m/min]

	//Discharge
	dam_TopLayerNonlinearExpSoilvel(y_i,global_params,params,qvs,state,ans);	//ans is used for convenience
	double qm = ans->ve[0] * 60.0;

	//Storage
	ans->ve[1] = (q_pl + q_sl) * A_h - qm;
	for(i=0;i<numparents;i++)
		ans->ve[1] += y_p[i]->ve[0] * 60.0;

	//Hillslope
	ans->ve[2] = rainfall - q_pl - q_pt - e_p;
	ans->ve[3] = q_pt - q_ts - e_t;
	ans->ve[4] = q_ts - q_sl - e_s;

	//Additional states
	ans->ve[5] = rainfall;
	ans->ve[6] = q_pl;
	ans->ve[7] = q_sl * A_h - q_b*60.0;
	for(i=0;i<numparents;i++)
		ans->ve[7] += y_p[i]->ve[7] * 60.0;
	ans->ve[7] *= v_B/L_i;
}

//Type 262
//Contains 2 states in the channel: dischage, storage
//Contains 3 layers on hillslope: ponded, top layer, soil
//Order of parameters: A_i,L_i,A_h,S_h,T_L,eta,h_b,k_D,k_dry,k_i | invtau,c_1,c_2,k_2
//The numbering is:	0   1   2   3   4   5   6   7   8     9  |   10    11  12  13
//Order of global_params: v_0,lambda_1,lambda_2,N,phi,v_B
//The numbering is:        0      1        2    3  4   5
void dam_TopLayerNonlinearExpSoilvel_ConstEta(VEC* y,VEC* global_params,VEC* params,QVSData* qvs,int state,VEC* ans)
{
	double q1,q2,S1,S2,S_max,q_max,S;

	//Parameters
	double lambda_1 = global_params->ve[1];
	double invtau = params->ve[10];

	//Find the discharge in [m^3/s]
	if(state == -1)
	{
		S = (y->ve[1] < 0.0) ? 0.0 : y->ve[1];
		//ans->ve[0] = invtau/60.0*pow(S,1.0/(1.0-lambda_1));
		ans->ve[0] = pow((1.0-lambda_1)*invtau/60.0 * S,1.0/(1.0-lambda_1));
	}
	else if(state == (int) qvs->n_values - 1)
	{
		S_max = qvs->points[qvs->n_values - 1][0];
		q_max = qvs->points[qvs->n_values - 1][1];
		ans->ve[0] = q_max;
	}
	else
	{
		S = (y->ve[1] < 0.0) ? 0.0 : y->ve[1];
		q2 = qvs->points[state+1][1];
		q1 = qvs->points[state][1];
		S2 = qvs->points[state+1][0];
		S1 = qvs->points[state][0];
		ans->ve[0] = (q2-q1)/(S2-S1) * (S-S1) + q1;
	}
}

//Type 262
//Contains 2 states in the channel: dischage, storage
//Contains 3 layers on hillslope: ponded, top layer, soil
//Order of parameters: A_i,L_i,A_h,S_h,T_L,eta,h_b,k_D,k_dry,k_i | invtau,c_1,c_2,k_2
//The numbering is:	0   1   2   3   4   5   6   7   8     9  |   10    11  12  13
//Order of global_params: v_0,lambda_1,lambda_2,N,phi,v_B
//The numbering is:        0      1        2    3  4   5
void TopLayerNonlinearExpSoilvel_ConstEta_Reservoirs(double t,VEC* y_i,VEC** y_p,unsigned short int numparents,VEC* global_params,double* forcing_values,QVSData* qvs,VEC* params,IVEC* iparams,int state,unsigned int** upstream,unsigned int* numupstream,VEC* ans)
{
	ans->ve[0] = forcing_values[2];
	ans->ve[1] = 0.0;
	ans->ve[2] = 0.0;
	ans->ve[3] = 0.0;
	ans->ve[4] = 0.0;
	ans->ve[5] = 0.0;
	ans->ve[6] = 0.0;
	ans->ve[7] = 0.0;
}

//Type 19
//Order of parameters: A_i,L_i,A_h,k2,k3,invtau,c_1,c_2
//The numbering is:	0   1   2   3  4    5    6   7
//Order of global_params: v_r,lambda_1,lambda_2,RC,v_h,v_g,e_pot
//The numbering is:        0      1        2     3  4   5    6
void LinearHillslope_Evap(double t,VEC* y_i,VEC** y_p,unsigned short int numparents,VEC* global_params,double* forcing_values,QVSData* qvs,VEC* params,IVEC* iparams,int state,unsigned int** upstream,unsigned int* numupstream,VEC* ans)
{
	unsigned short int i;

	double lambda_1 = global_params->ve[1];
	double e_pot = global_params->ve[6] * (1e-3/60.0);	//[mm/hr]->[m/min]

	double A_h = params->ve[2];
	double k2 = params->ve[3];
	double k3 = params->ve[4];
	double invtau = params->ve[5];
	double c_1 = params->ve[6];
	double c_2 = params->ve[7];

	double q = y_i->ve[0];		//[m^3/s]
	double s_p = y_i->ve[1];	//[m]
	double s_a = y_i->ve[2];	//[m]

	double q_pl = k2 * s_p;
	double q_al = k3 * s_a;

	//Evaporation
	double C_p,C_a,C_T,Corr_evap;
	if(e_pot > 0.0)
	{
		C_p = s_p / e_pot;
		C_a = s_a / e_pot;
		C_T = C_p + C_a;
	}
	else
	{
		C_p = 0.0;
		C_a = 0.0;
		C_T = 0.0;
	}

	//Corr_evap = (!state && C_T > 0.0) ? 1.0/C_T : 1.0;
	Corr_evap = (C_T > 1.0) ? 1.0/C_T : 1.0;

	double e_p = Corr_evap * C_p * e_pot;
	double e_a = Corr_evap * C_a * e_pot;

	//Discharge
	ans->ve[0] = -q + (q_pl + q_al) * A_h/60.0;
	for(i=0;i<numparents;i++)
		ans->ve[0] += y_p[i]->ve[0];
	ans->ve[0] = invtau * pow(q,lambda_1) * ans->ve[0];

	//Hillslope
	ans->ve[1] = forcing_values[0]*c_1 - q_pl - e_p;
	ans->ve[2] = forcing_values[0]*c_2 - q_al - e_a;
}


//Type 19
//State = 0: C_T >= 1.0
//State = 1: C_T < 1.0
int LinearHillslope_Evap_Check(VEC* y,VEC* params,VEC* global_params,QVSData* qvs,unsigned int dam)
{
	double s_p = y->ve[1];	//[m]
	double s_a = y->ve[2];	//[m]

	double e_pot = global_params->ve[6] * (1e-3/60.0);	//[mm/hr]->[m/min]

	double C_p,C_a,C_T;
	if(e_pot > 0.0)
	{
		C_p = s_p / e_pot;
		C_a = s_a / e_pot;
		C_T = C_p + C_a;
	}
	else
		return 0;

	if(C_T >= 1.0)	return 0;
	else 		return 1;
}



//Type 190
//Order of parameters: A_i,L_i,A_h,k2,k3,invtau,c_1,c_2
//The numbering is:	0   1   2   3  4    5    6   7
//Order of global_params: v_r,lambda_1,lambda_2,RC,v_h,v_g
//The numbering is:        0      1        2     3  4   5 
void LinearHillslope_MonthlyEvap(double t,VEC* y_i,VEC** y_p,unsigned short int numparents,VEC* global_params,double* forcing_values,QVSData* qvs,VEC* params,IVEC* iparams,int state,unsigned int** upstream,unsigned int* numupstream,VEC* ans)
{
	unsigned short int i;

	double lambda_1 = global_params->ve[1];
	//double e_pot = global_params->ve[6] * (1e-3/60.0);	//[mm/hr]->[m/min]

	double A_h = params->ve[2];
	double k2 = params->ve[3];
	double k3 = params->ve[4];
	double invtau = params->ve[5];
	double c_1 = params->ve[6];
	double c_2 = params->ve[7];

	double q = y_i->ve[0];		//[m^3/s]
	double s_p = y_i->ve[1];	//[m]
	double s_a = y_i->ve[2];	//[m]

	double q_pl = k2 * s_p;
	double q_al = k3 * s_a;

	//Evaporation
	double C_p,C_a,C_T,Corr_evap;
	//double e_pot = forcing_values[1] * (1e-3/60.0);
	double e_pot = forcing_values[1] * (1e-3/(30.0*24.0*60.0));	//[mm/month] -> [m/min]

	if(e_pot > 0.0)
	{
		C_p = s_p / e_pot;
		C_a = s_a / e_pot;
		C_T = C_p + C_a;
	}
	else
	{
		C_p = 0.0;
		C_a = 0.0;
		C_T = 0.0;
	}

	//Corr_evap = (!state && C_T > 0.0) ? 1.0/C_T : 1.0;
	Corr_evap = (C_T > 1.0) ? 1.0/C_T : 1.0;

	double e_p = Corr_evap * C_p * e_pot;
	double e_a = Corr_evap * C_a * e_pot;

	//Discharge
	ans->ve[0] = -q + (q_pl + q_al) * A_h/60.0;
	for(i=0;i<numparents;i++)
		ans->ve[0] += y_p[i]->ve[0];
	ans->ve[0] = invtau * pow(q,lambda_1) * ans->ve[0];

	//Hillslope
	ans->ve[1] = forcing_values[0]*c_1 - q_pl - e_p;
	ans->ve[2] = forcing_values[0]*c_2 - q_al - e_a;
}


//Type 191
//Order of parameters: A_i,L_i,A_h,k2,k3,invtau,c_1,c_2
//The numbering is:	0   1   2   3  4    5    6   7
//Order of global_params: v_r,lambda_1,lambda_2,RC,v_h,v_g,v_B
//The numbering is:        0      1        2     3  4   5   6
void LinearHillslope_MonthlyEvap_extras(double t,VEC* y_i,VEC** y_p,unsigned short int numparents,VEC* global_params,double* forcing_values,QVSData* qvs,VEC* params,IVEC* iparams,int state,unsigned int** upstream,unsigned int* numupstream,VEC* ans)
{
	unsigned short int i;

	double lambda_1 = global_params->ve[1];
	double v_B = global_params->ve[6];

	double L = params->ve[1];
	double A_h = params->ve[2];
	double k2 = params->ve[3];
	double k3 = params->ve[4];
	double invtau = params->ve[5];
	double c_1 = params->ve[6];
	double c_2 = params->ve[7];

	double q = y_i->ve[0];		//[m^3/s]
	double s_p = y_i->ve[1];	//[m]
	double s_a = y_i->ve[2];	//[m]
	double q_b = y_i->ve[5];	//[m^3/s]

	double q_pl = k2 * s_p;
	double q_al = k3 * s_a;

	//Evaporation
	double C_p,C_a,C_T,Corr_evap;
	//double e_pot = forcing_values[1] * (1e-3/60.0);
	double e_pot = forcing_values[1] * (1e-3/(30.0*24.0*60.0));	//[mm/month] -> [m/min]

	if(e_pot > 0.0)
	{
		C_p = s_p / e_pot;
		C_a = s_a / e_pot;
		C_T = C_p + C_a;
	}
	else
	{
		C_p = 0.0;
		C_a = 0.0;
		C_T = 0.0;
	}

	//Corr_evap = (!state && C_T > 0.0) ? 1.0/C_T : 1.0;
	Corr_evap = (C_T > 1.0) ? 1.0/C_T : 1.0;

	double e_p = Corr_evap * C_p * e_pot;
	double e_a = Corr_evap * C_a * e_pot;

	//Discharge
	ans->ve[0] = -q + (q_pl + q_al) * A_h/60.0;
	for(i=0;i<numparents;i++)
		ans->ve[0] += y_p[i]->ve[0];
	ans->ve[0] = invtau * pow(q,lambda_1) * ans->ve[0];

	//Hillslope
	ans->ve[1] = forcing_values[0]*c_1 - q_pl - e_p;
	ans->ve[2] = forcing_values[0]*c_2 - q_al - e_a;

	//Additional states
	ans->ve[3] = forcing_values[0]*c_1;
	ans->ve[4] = q_pl;
	ans->ve[5] = q_al * A_h - q_b * 60.0;
	for(i=0;i<numparents;i++)
		ans->ve[5] += y_p[i]->ve[5] * 60.0;
		//ans->ve[5] += k3*y_p[i]->ve[2]*A_h;
	ans->ve[5] *= v_B/L;
}


//Type 191
//Order of parameters: A_i,L_i,A_h,k2,k3,invtau,c_1,c_2
//The numbering is:	0   1   2   3  4    5    6   7
//Order of global_params: v_r,lambda_1,lambda_2,RC,v_h,v_g,v_B
//The numbering is:        0      1        2     3  4   5   6
void LinearHillslope_Reservoirs_extras(double t,VEC* y_i,VEC** y_p,unsigned short int numparents,VEC* global_params,double* forcing_values,QVSData* qvs,VEC* params,IVEC* iparams,int state,unsigned int** upstream,unsigned int* numupstream,VEC* ans)
{
	ans->ve[0] = forcing_values[2];
	ans->ve[1] = 0.0;
	ans->ve[2] = 0.0;
	ans->ve[3] = 0.0;
	ans->ve[4] = 0.0;
	ans->ve[5] = 0.0;
	ans->ve[6] = 0.0;
}


//Type 20
//Order of parameters: A_i,L_i,A_h,invtau,c_1,c_2
//The numbering is:	0   1   2    3     4   5
//Order of global_params: v_r,lambda_1,lambda_2,beta,k_p,k_a,theta_p,theta_a,scale_p,scale_a
//The numbering is:        0      1        2     3    4   5     6       7      8        9

//Order of parameters: A_i,L_i,A_h,k2,k3,invtau,c_1,c_2
//The numbering is:	0   1   2   3  4    5    6   7
//Order of global_params: v_r,lambda_1,lambda_2,RC,v_h,v_g
//The numbering is:        0      1        2     3  4   5 
void Hillslope_Toy(double t,VEC* y_i,VEC** y_p,unsigned short int numparents,VEC* global_params,double* forcing_values,QVSData* qvs,VEC* params,IVEC* iparams,int state,unsigned int** upstream,unsigned int* numupstream,VEC* ans)
{

	unsigned short int i;
	double lambda_1 = global_params->ve[1];


	//For complex
	double omega = global_params->ve[3];
	double alpha = global_params->ve[4];
	double beta = global_params->ve[5];
	double A_p = global_params->ve[6];
	double A_a = global_params->ve[7];

/*
	//For gamma
	double beta = global_params->ve[3];
	double k_p = global_params->ve[4];
	double k_a = global_params->ve[5];
	double theta_p = global_params->ve[6];
	double theta_a = global_params->ve[7];
	double scale_p = global_params->ve[8];
	double scale_a = global_params->ve[9];
*/
/*
	//For lognormal
	double mu = global_params->ve[3];
	double sigma2 = global_params->ve[4];
*/
/*
	//For exponential	
	double lambda = global_params->ve[3];
	double scale = global_params->ve[4];
*/
/*
	//For line and drop
	double t_p = global_params->ve[3];
	double p = global_params->ve[4];
	double alpha = global_params->ve[5];
*/
/*
	//For GEV
	double mu = global_params->ve[3];
	double sigma = global_params->ve[4];
	double ksi = global_params->ve[5];
	double shift = global_params->ve[6];
*/
/*
	//For linear
	double k2 = global_params->ve[3];
	double k3 = global_params->ve[4];
	double scale = global_params->ve[5];
*/
	//General use
	double A_h = params->ve[2];
	double invtau = params->ve[3];
	double c_1 = params->ve[4];
	double c_2 = params->ve[5];

	//States
	double q = y_i->ve[0];
	double s_p = y_i->ve[1];
	double s_a = y_i->ve[2];

	//For Complex Model
	double Q_pl = A_p * pow(sin(omega*s_p),2.0);
	double Q_al = A_a * (exp(alpha*s_a) - 1.0);
	double Q_pa = beta * s_p;

/*
	//For Gamma Forcing Model
	//double Q_pl = A_p * pow(sin(omega*s_p),2.0);
	//double Q_al = A_a * (exp(alpha*s_a) - 1.0);
	double Q_pl = scale_p * c_1 * pow(t,k_p-1.0) * exp(-t/theta_p); //[m/min]
	//double Q_al = 0.0;
	double Q_al = scale_a * c_2 * pow(t,k_a-1.0) * exp(-t/theta_a); //[m/min]
	double Q_pa = beta * s_p; //[m/min]
*/
/*
	//For Log Normal Forcing Model
	double Q_pl = (t <= 0.0) ? 0.0 : c_1/t * exp( -sq(log(t)-mu)/(2.0*sigma2) );
	double Q_al = 0.0;
*/
/*
	//For Exponential Forcing Model
	double Q_pl = scale * lambda * exp(-lambda*t);
	double Q_al = 0.0;
*/
/*
	//For Line and Drop Model
	double Q_pl = (t <= t_p) ? p/t_p * t : p * exp(-alpha*(t-t_p));
	double Q_al = 0.0;
*/
/*
	//For GEV
	double holder = pow(1.0 + (t-shift-mu)/sigma*ksi,-1.0/ksi);
	//double holder = (t < 25.0) ? pow(1.0 + (t-shift-mu)/sigma*ksi,-1.0/ksi) : 0.0;
	//double holder = (t < 6.0 + shift) ? pow(1.0 + (t-shift-mu)/sigma*ksi,-1.0/ksi) : 0.0;
	//double holder = (t < 1000.0 + shift) ? pow(1.0 + (t-shift-mu)/sigma*ksi,-1.0/ksi) : 0.0;
	double Q_pl = 1.0/sigma * pow(holder,ksi+1) * exp(-holder);
	double Q_al = 0.0;
*/
/*
	//For Linear Model
	double Q_pl = scale * (k2 * 0.015 * exp(-k2*t) + k3 * 0.015 * exp(-k3*t));
	double Q_al = 0.0;
*/
	//Discharge
	ans->ve[0] = -q + (Q_pl + Q_al) * A_h/60.0;
	for(i=0;i<numparents;i++)
		ans->ve[0] += y_p[i]->ve[0];
	ans->ve[0] = invtau * pow(q,lambda_1) * ans->ve[0];

	//Hillslope
	//ans->ve[1] = 0.0;
	//ans->ve[2] = 0.0;
	ans->ve[1] = forcing_values[0]*(0.001/60.0) - Q_pl - Q_pa;
	ans->ve[2] = -Q_al + Q_pa;

}

//Type 60
//Order of parameters: A_i,L_i,A_h,h_b,k2,k3,invtau,c_1
//The numbering is:	0   1   2   3   4  5    6    7 
//Order of global_params: v_r,lambda_1,lambda_2,v_h,v_g,e_pot
//The numbering is:        0      1        2     3   4   5  
void LinearHillslope_Evap_RC(double t,VEC* y_i,VEC** y_p,unsigned short int numparents,VEC* global_params,double* forcing_values,QVSData* qvs,VEC* params,IVEC* iparams,int state,unsigned int** upstream,unsigned int* numupstream,VEC* ans)
{
	unsigned short int i;

	double lambda_1 = global_params->ve[1];
	double e_pot = global_params->ve[5] * (1e-3/60.0);	//[mm/hr]->[m/min]

	double A_h = params->ve[2];	//[m^2]
	double h_b = params->ve[3];	//[m]
	double k2 = params->ve[4];
	double k3 = params->ve[5];
	double invtau = params->ve[6];
	double c_1 = params->ve[7];

	double q = y_i->ve[0];		//[m^3/s]
	double s_p = y_i->ve[1];	//[m]
	double s_a = y_i->ve[2];	//[m]

	double q_pl = k2 * s_p;
	double q_al = k3 * s_a;

	//Run off
	double d_soil = (s_a > h_b) ? 0.0 : h_b - s_a;
	double RC = (s_p > 0.0) ? (s_p*(s_p+2.0*d_soil))/sq(s_p+d_soil) : 0.0;

	//Evaporation
	double C_p,C_a,C_T,Corr_evap;
	if(e_pot > 0.0)
	{
		C_p = s_p / e_pot;
		C_a = s_a / e_pot;
		C_T = C_p + C_a;
	}
	else
	{
		C_p = 0.0;
		C_a = 0.0;
		C_T = 0.0;
	}
	Corr_evap = (C_T > 1.0) ? 1.0/C_T : 1.0;

	double e_p = Corr_evap * C_p * e_pot;
	double e_a = Corr_evap * C_a * e_pot;

	//Discharge
	ans->ve[0] = -q + (q_pl + q_al) * A_h/60.0;
	for(i=0;i<numparents;i++)
		ans->ve[0] += y_p[i]->ve[0];
	ans->ve[0] = invtau * pow(q,lambda_1) * ans->ve[0];

	//Hillslope
	ans->ve[1] = forcing_values[0]*c_1*RC - q_pl - e_p;
	ans->ve[2] = forcing_values[0]*c_1*(1.0-RC) - q_al - e_a;
}


// Type 21 ************************************************

//Type 21
//Order of parameters: A_i,L_i,A_h,k2,k3,invtau
//The numbering is:	0   1   2  3  4    5
//Order of global_params: v_r,lambda_1,lambda_2,RC,S_0,v_h,v_g
//The numbering is:        0      1        2     3  4   5   6
void nodam_rain_hillslope(double t,VEC* y_i,VEC** y_p,unsigned short int numparents,VEC* global_params,double* forcing_values,QVSData* qvs,VEC* params,IVEC* iparams,int state,unsigned int** upstream,unsigned int* numupstream,VEC* ans)
{
	unsigned short int i;
	double lambda_1 = global_params->ve[1];
	double RC = global_params->ve[3];
	double A_h = params->ve[2];
	double k2 = params->ve[3];
	double k3 = params->ve[4];
	double invtau = params->ve[5];
	//double F_et = 1.0;
	double F_et = 0.05;

	double S = y_i->ve[1];
	double Ss = y_i->ve[2];
	double Sg = y_i->ve[3];

	//double q = invtau*S;
	double qm = invtau*pow(S,1.0/(1.0-lambda_1));

	ans->ve[1] = k2*Ss + k3*Sg - qm;
	for(i=0;i<numparents;i++)	ans->ve[1] += y_p[i]->ve[0] * 60.0;
	//ans->ve[1] = ans->ve[1];

	ans->ve[2] = RC*forcing_values[0]*A_h*(.001/60.0) - k2*Ss;
	ans->ve[3] = (1.0-RC)*F_et*forcing_values[0]*A_h*(.001/60.0) - k3*Sg;
	//ans->ve[3] = 0.0;
}


//Type 21
//Order of parameters: A_i,L_i,A_h,k2,k3,invtau,orifice_area,H_spill,H_max,S_max,alpha,orifice_diam,c_1,c_2,L_spill
//The numbering is:	0   1   2  3  4    5	       6      7       8     9	  10	    11       12  13  14
//Order of global_params: v_r,lambda_1,lambda_2,RC,S_0,v_h,v_g
//The numbering is:        0      1        2     3  4   5   6
void dam_rain_hillslope(double t,VEC* y_i,VEC** y_p,unsigned short int numparents,VEC* global_params,double* forcing_values,QVSData* qvs,VEC* params,IVEC* iparams,int state,unsigned int** upstream,unsigned int* numupstream,VEC* ans)
{
	unsigned short int i;
	double lambda_1 = global_params->ve[1];
	double RC = global_params->ve[3];
	double A_h = params->ve[2];
	double k2 = params->ve[3];
	double k3 = params->ve[4];
	double invtau = params->ve[5];
	double orifice_area = params->ve[6];
	double H_spill = params->ve[7];
	double H_max = params->ve[8];
	double S_max = params->ve[9];
	double alpha = params->ve[10];
	double diam = params->ve[11];
	double c_1 = params->ve[12];
	double c_2 = params->ve[13];
	double L_spill = params->ve[14];
	double g = 9.81;
	double F_et = 1.0;

	double qm;
	double S = y_i->ve[1];
	double Ss = y_i->ve[2];
	double Sg = y_i->ve[3];
	double h = H_max * pow(S/S_max,alpha);
	double diff = (h - H_spill >= 0) ? h - H_spill : 0.0;

	if(state == 0)
		//q = invtau*S;
		qm = invtau*pow(S,1.0/(1.0-lambda_1));
	else if(state == 1)
		qm = c_1*orifice_area*pow(2*g*h,.5) * 60.0;
		//qm = 0.0;
	else if(state == 2)
		qm = c_1*orifice_area*pow(2*g*h,.5) * 60.0 + c_2*L_spill*pow(diff,1.5) * 60.0;
	else if(state == 3)
		qm = c_1*orifice_area*pow(2*g*h,.5) * 60.0 + c_2*L_spill*pow(diff,1.5) * 60.0 + invtau*pow(S-S_max,1.0/(1.0-lambda_1));
	else if(state == 4)
	{
		double r = diam/2.0;
		double frac = (h < 2*r) ? (h-r)/r : 1.0;
		double A = -r*r*(acos(frac) - pow(1-frac*frac,.5)*frac - 3.1415926535897932);
		qm = c_1*A*pow(2*g*h,.5) * 60.0;
	}
	else printf("Error in function evaluation. Invalid state %u.\n",state);

	ans->ve[1] = k2*Ss + k3*Sg - qm;
	for(i=0;i<numparents;i++)	ans->ve[1] += y_p[i]->ve[0] * 60.0;
	//ans->ve[1] = ans->ve[1];

	ans->ve[2] = RC*forcing_values[0]*A_h*(.001/60.0) - k2*Ss;
	ans->ve[3] = (1.0-RC)*F_et*forcing_values[0]*A_h*(.001/60.0) - k3*Sg;
}

//Type 21
//Order of parameters: A_i,L_i,A_h,k2,k3,invtau,orifice_area,H_spill,H_max,S_max,alpha,orifice_diam,c_1,c_2,L_spill
//The numbering is:	0   1   2  3  4    5	       6      7       8     9	  10	    11       12  13  14
//Order of global_params: v_r,lambda_1,lambda_2,RC,S_0,v_h,v_g
//The numbering is:        0      1        2     3  4   5   6
void dam_q(VEC* y,VEC* global_params,VEC* params,QVSData* qvs,int state,VEC* ans)
{
	double lambda_1 = global_params->ve[1];
	double invtau = params->ve[5];
	double S = (y->ve[1] < 0.0) ? 0.0 : y->ve[1];

	if(state == 0)
		//ans->ve[0] = invtau*S;
		ans->ve[0] = invtau/60.0*pow(S,1.0/(1.0-lambda_1));
	else
	{
		double orifice_area = params->ve[6];
		double H_spill = params->ve[7];
		double H_max = params->ve[8];
		double S_max = params->ve[9];
		double alpha = params->ve[10];
		double diam = params->ve[11];
		double c_1 = params->ve[12];
		double c_2 = params->ve[13];
		double L_spill = params->ve[14];
		double g = 9.81;

		double h = H_max * pow(S/S_max,alpha);
		double diff = (h - H_spill >= 0) ? h - H_spill : 0.0;

		if(state == 1)
			ans->ve[0] = c_1*orifice_area*pow(2*g*h,.5);
			//ans->ve[0] = 0.0;
		else if(state == 2)
			ans->ve[0] = c_1*orifice_area*pow(2*g*h,.5) + c_2*L_spill*pow(diff,1.5);
		else if(state == 3)
			ans->ve[0] = c_1*orifice_area*pow(2*g*h,.5) + c_2*L_spill*pow(diff,1.5) + invtau/60.0*pow(S-S_max,1.0/(1.0-lambda_1));
		else if(state == 4)
		{
			double r = diam/2.0;
			double frac = (h < 2*r) ? (h-r)/r : 1.0;
			double A = -r*r*(acos(frac) - pow(1.0-frac*frac,.5)*frac - 3.141592653589);
			ans->ve[0] = c_1*A*pow(2*g*h,.5);
		}
		else printf("Error in dam evaluation. Invalid state %u.\n",state);
	}
}


//Type 21
//Order of parameters: A_i,L_i,A_h,k2,k3,invtau,orifice_area,H_spill,H_max,S_max,alpha,orifice_diam,c_1,c_2,L_spill
//The numbering is:	0   1   2  3  4    5	       6      7       8     9	  10	    11       12  13  14
//Order of global_params: v_r,lambda_1,lambda_2,RC,S_0,v_h,v_g
//The numbering is:        0      1        2     3  4   5   6
int dam_check(VEC* y,VEC* global_params,VEC* params,QVSData* qvs,unsigned int dam)
{
	if(dam == 0)	return 0;

	double H_spill = params->ve[7];
	double H_max = params->ve[8];
	double S_max = params->ve[9];
	double alpha = params->ve[10];
	double diam = params->ve[11];
	double S = y->ve[1];
	double h = H_max * pow(S/S_max,alpha);

	if(h < diam)		return 4;
	if(h <= H_spill)	return 1;
	if(h <= H_max)		return 2;
	return 3;
}

// Type 22 ************************************************

//Type 22
//Order of parameters: A_i,L_i,A_h,RC,v_h,v_r,k2,k3,invtau
//The numbering is:	0   1   2   3  4   5   6  7    8
//Order of global_params: lambda_1,lambda_2,S_0,v_g
//The numbering is:         0        1       2   3
void nodam_rain_hillslope2(double t,VEC* y_i,VEC** y_p,unsigned short int numparents,VEC* global_params,double* forcing_values,QVSData* qvs,VEC* params,IVEC* iparams,int state,unsigned int** upstream,unsigned int* numupstream,VEC* ans)
{
	unsigned short int i;
	double lambda_1 = global_params->ve[0];
	double A_h = params->ve[2];
	double RC = params->ve[3];
	double k2 = params->ve[6];
	double k3 = params->ve[7];
	double invtau = params->ve[8];

	double S = y_i->ve[1];
	double Ss = y_i->ve[2];
	double Sg = y_i->ve[3];

	//double q = invtau*S;
	double qm = invtau*pow(S,1.0/(1.0-lambda_1));

	ans->ve[1] = k2*Ss + k3*Sg - qm;
	for(i=0;i<numparents;i++)	ans->ve[1] += y_p[i]->ve[0] * 60.0;
	//ans->ve[1] = ans->ve[1];

	ans->ve[2] = RC*forcing_values[0]*A_h*(.001/60.0) - k2*Ss;
	ans->ve[3] = (1.0-RC)*forcing_values[0]*A_h*(.001/60.0) - k3*Sg;
	//ans->ve[3] = 0.0;
}


//Type 22
//Order of parameters: A_i,L_i,A_h,RC,v_h,v_r,k2,k3,invtau,orifice_area,H_spill,H_max,S_max,alpha,orifice_diam,c_1,c_2,L_spill
//The numbering is:	0   1   2  3   4   5   6  7   8          9	  10	  11   12     13      14        15  16   17
//Order of global_params: lambda_1,lambda_2,S_0,v_g
//The numbering is:         0        1       2   3
void dam_rain_hillslope2(double t,VEC* y_i,VEC** y_p,unsigned short int numparents,VEC* global_params,double* forcing_values,QVSData* qvs,VEC* params,IVEC* iparams,int state,unsigned int** upstream,unsigned int* numupstream,VEC* ans)
{
	unsigned short int i;
	double lambda_1 = global_params->ve[0];
	double A_h = params->ve[2];
	double RC = params->ve[3];
	double k2 = params->ve[6];
	double k3 = params->ve[7];
	double invtau = params->ve[8];
	double orifice_area = params->ve[9];
	double H_spill = params->ve[10];
	double H_max = params->ve[11];
	double S_max = params->ve[12];
	double alpha = params->ve[13];
	double diam = params->ve[14];
	double c_1 = params->ve[15];
	double c_2 = params->ve[16];
	double L_spill = params->ve[17];
	double g = 9.81;

	double qm;
	double S = y_i->ve[1];
	double Ss = y_i->ve[2];
	double Sg = y_i->ve[3];
	double h = H_max * pow(S/S_max,alpha);
	double diff = (h - H_spill >= 0) ? h - H_spill : 0.0;

	if(state == 0)
		//q = invtau*S;
		qm = invtau*pow(S,1.0/(1.0-lambda_1));
	else if(state == 1)
		qm = c_1*orifice_area*pow(2*g*h,.5) * 60.0;
		//qm = 0.0;
	else if(state == 2)
		qm = c_1*orifice_area*pow(2*g*h,.5) * 60.0 + c_2*L_spill*pow(diff,1.5) * 60.0;
	else if(state == 3)
		qm = c_1*orifice_area*pow(2*g*h,.5) * 60.0 + c_2*L_spill*pow(diff,1.5) * 60.0 + invtau*pow(S-S_max,1.0/(1.0-lambda_1));
	else if(state == 4)
	{
		double r = diam/2.0;
		double frac = (h < 2*r) ? (h-r)/r : 1.0;
		double A = -r*r*(acos(frac) - pow(1-frac*frac,.5)*frac - 3.1415926535897932);
		qm = c_1*A*pow(2*g*h,.5) * 60.0;
	}
	else printf("Error in function evaluation. Invalid state %u.\n",state);

	ans->ve[1] = k2*Ss + k3*Sg - qm;
	for(i=0;i<numparents;i++)	ans->ve[1] += y_p[i]->ve[0] * 60.0;
	//ans->ve[1] = ans->ve[1];

	ans->ve[2] = RC*forcing_values[0]*A_h*(.001/60.0) - k2*Ss;
	ans->ve[3] = (1.0-RC)*forcing_values[0]*A_h*(.001/60.0) - k3*Sg;
}


//Type 22
//Order of parameters: A_i,L_i,A_h,RC,v_h,v_r,k2,k3,invtau,orifice_area,H_spill,H_max,S_max,alpha,orifice_diam,c_1,c_2,L_spill
//The numbering is:	0   1   2  3   4   5   6  7   8          9	  10	  11   12     13      14        15  16   17
//Order of global_params: lambda_1,lambda_2,S_0,v_g
//The numbering is:         0        1       2   3
void dam_q2(VEC* y,VEC* global_params,VEC* params,QVSData* qvs,int state,VEC* ans)
{
	double lambda_1 = global_params->ve[0];
	double invtau = params->ve[8];
	double S = (y->ve[1] < 0.0) ? 0.0 : y->ve[1];

	if(state == 0)
		//ans->ve[0] = invtau*S;
		ans->ve[0] = invtau/60.0*pow(S,1.0/(1.0-lambda_1));
	else
	{
		double orifice_area = params->ve[9];
		double H_spill = params->ve[10];
		double H_max = params->ve[11];
		double S_max = params->ve[12];
		double alpha = params->ve[13];
		double diam = params->ve[14];
		double c_1 = params->ve[15];
		double c_2 = params->ve[16];
		double L_spill = params->ve[17];
		double g = 9.81;

		double h = H_max * pow(S/S_max,alpha);
		double diff = (h - H_spill >= 0) ? h - H_spill : 0.0;

		if(state == 1)
			ans->ve[0] = c_1*orifice_area*pow(2*g*h,.5);
			//ans->ve[0] = 0.0;
		else if(state == 2)
			ans->ve[0] = c_1*orifice_area*pow(2*g*h,.5) + c_2*L_spill*pow(diff,1.5);
		else if(state == 3)
			ans->ve[0] = c_1*orifice_area*pow(2*g*h,.5) + c_2*L_spill*pow(diff,1.5) + invtau/60.0*pow(S-S_max,1.0/(1.0-lambda_1));
		else if(state == 4)
		{
			double r = diam/2.0;
			double frac = (h < 2*r) ? (h-r)/r : 1.0;
			double A = -r*r*(acos(frac) - pow(1.0-frac*frac,.5)*frac - 3.141592653589);
			ans->ve[0] = c_1*A*pow(2*g*h,.5);
		}
		else printf("Error in dam evaluation. Invalid state %u.\n",state);
	}
}


//Type 22
//Order of parameters: A_i,L_i,A_h,RC,v_h,v_r,k2,k3,invtau,orifice_area,H_spill,H_max,S_max,alpha,orifice_diam,c_1,c_2,L_spill
//The numbering is:	0   1   2  3   4   5   6  7   8          9	  10	  11   12     13      14        15  16   17
//Order of global_params: lambda_1,lambda_2,S_0,v_g
//The numbering is:         0        1       2   3
int dam_check2(VEC* y,VEC* global_params,VEC* params,QVSData* qvs,unsigned int dam)
{
	if(dam == 0)	return 0;

	double H_spill = params->ve[10];
	double H_max = params->ve[11];
	double S_max = params->ve[12];
	double alpha = params->ve[13];
	double diam = params->ve[14];
	double S = y->ve[1];
	double h = H_max * pow(S/S_max,alpha);

	if(h < diam)		return 4;
	if(h <= H_spill)	return 1;
	if(h <= H_max)		return 2;
	return 3;
}

// Type 23 ************************************************

//Type 23
//Order of parameters: A_i,L_i,A_h,RC,v_h,v_r,k2,k3,invtau
//The numbering is:	0   1   2   3  4   5   6  7    8
//Order of global_params: lambda_1,lambda_2,S_0,v_g
//The numbering is:         0        1       2   3
void nodam_rain_hillslope3(double t,VEC* y_i,VEC** y_p,unsigned short int numparents,VEC* global_params,double* forcing_values,QVSData* qvs,VEC* params,IVEC* iparams,int state,unsigned int** upstream,unsigned int* numupstream,VEC* ans)
{
	unsigned short int i;
	double lambda_1 = global_params->ve[0];
	double A_h = params->ve[2];
	double RC = params->ve[3];
	double k2 = params->ve[6];
	double k3 = params->ve[7];
	double invtau = params->ve[8];

	double S = y_i->ve[1];
	double Ss = y_i->ve[2];
	double Sg = y_i->ve[3];

	//double q = invtau*S;
	double qm = invtau*pow(S,1.0/(1.0-lambda_1));

	ans->ve[1] = k2*Ss + k3*Sg - qm;
	for(i=0;i<numparents;i++)	ans->ve[1] += y_p[i]->ve[0] * 60.0;
	//ans->ve[1] = ans->ve[1];

	ans->ve[2] = RC*forcing_values[0]*A_h*(.001/60.0) - k2*Ss;
	ans->ve[3] = (1.0-RC)*forcing_values[0]*A_h*(.001/60.0) - k3*Sg;
	//ans->ve[3] = 0.0;
}


//Type 23
//Order of parameters: A_i,L_i,A_h,RC,v_h,v_r,k2,k3,invtau,orifice_area,H_spill,H_max,S_max,alpha,orifice_diam,c_1,c_2,L_spill
//The numbering is:	0   1   2  3   4   5   6  7   8          9	  10	  11   12     13      14        15  16   17
//Order of global_params: lambda_1,lambda_2,S_0,v_g
//The numbering is:         0        1       2   3
void dam_rain_hillslope3(double t,VEC* y_i,VEC** y_p,unsigned short int numparents,VEC* global_params,double* forcing_values,QVSData* qvs,VEC* params,IVEC* iparams,int state,unsigned int** upstream,unsigned int* numupstream,VEC* ans)
{
	unsigned short int i;
	double lambda_1 = global_params->ve[0];
	double A_h = params->ve[2];
	double RC = params->ve[3];
	double k2 = params->ve[6];
	double k3 = params->ve[7];
	double invtau = params->ve[8];
	double orifice_area = params->ve[9];
	double H_spill = params->ve[10];
	double H_max = params->ve[11];
	double S_max = params->ve[12];
	double alpha = params->ve[13];
	double diam = params->ve[14];
	double c_1 = params->ve[15];
	double c_2 = params->ve[16];
	double L_spill = params->ve[17];
	double g = 9.81;

	double qm;
	double S = y_i->ve[1];
	double Ss = y_i->ve[2];
	double Sg = y_i->ve[3];
	double h = H_max * pow(S/S_max,alpha);
	double diff = (h - H_spill >= 0) ? h - H_spill : 0.0;

	if(state == 0)
		//q = invtau*S;
		qm = invtau*pow(S,1.0/(1.0-lambda_1));
	else if(state == 1)
		qm = c_1*orifice_area*pow(2*g*h,.5) * 60.0;
		//qm = 0.0;
	else if(state == 2)
		qm = c_1*orifice_area*pow(2*g*h,.5) * 60.0 + c_2*L_spill*pow(diff,1.5) * 60.0;
	else if(state == 3)
		qm = c_1*orifice_area*pow(2*g*h,.5) * 60.0 + c_2*L_spill*pow(diff,1.5) * 60.0 + invtau*pow(S-S_max,1.0/(1.0-lambda_1));
	else if(state == 4)
	{
//if(h < 1e-6)	qm = 0.0;
		qm = c_2 * 2*sqrt(h*diam-h*h) * pow(h,1.5) * 60.0;
	}
	else if(state == 5)
		qm = c_2 * diam * pow(h,1.5) * 60.0;
	else printf("Error in function evaluation. Invalid state %u.\n",state);

	ans->ve[1] = k2*Ss + k3*Sg - qm;
	for(i=0;i<numparents;i++)	ans->ve[1] += y_p[i]->ve[0] * 60.0;
	//ans->ve[1] = ans->ve[1];

	ans->ve[2] = RC*forcing_values[0]*A_h*(.001/60.0) - k2*Ss;
	ans->ve[3] = (1.0-RC)*forcing_values[0]*A_h*(.001/60.0) - k3*Sg;
}


//Type 23
//Order of parameters: A_i,L_i,A_h,RC,v_h,v_r,k2,k3,invtau,orifice_area,H_spill,H_max,S_max,alpha,orifice_diam,c_1,c_2,L_spill
//The numbering is:	0   1   2  3   4   5   6  7   8          9	  10	  11   12     13      14        15  16   17
//Order of global_params: lambda_1,lambda_2,S_0,v_g
//The numbering is:         0        1       2   3
void dam_q3(VEC* y,VEC* global_params,VEC* params,QVSData* qvs,int state,VEC* ans)
{
	double lambda_1 = global_params->ve[0];
	double invtau = params->ve[8];
	double S = (y->ve[1] < 0.0) ? 0.0 : y->ve[1];

	if(state == 0)
		//ans->ve[0] = invtau*S;
		ans->ve[0] = invtau/60.0*pow(S,1.0/(1.0-lambda_1));
	else
	{
		double orifice_area = params->ve[9];
		double H_spill = params->ve[10];
		double H_max = params->ve[11];
		double S_max = params->ve[12];
		double alpha = params->ve[13];
		double diam = params->ve[14];
		double c_1 = params->ve[15];
		double c_2 = params->ve[16];
		double L_spill = params->ve[17];
		double g = 9.81;

		double h = H_max * pow(S/S_max,alpha);
		double diff = (h - H_spill >= 0) ? h - H_spill : 0.0;

		if(state == 1)
			ans->ve[0] = c_1*orifice_area*pow(2*g*h,.5);
			//ans->ve[0] = 0.0;
		else if(state == 2)
			ans->ve[0] = c_1*orifice_area*pow(2*g*h,.5) + c_2*L_spill*pow(diff,1.5);
		else if(state == 3)
			ans->ve[0] = c_1*orifice_area*pow(2*g*h,.5) + c_2*L_spill*pow(diff,1.5) + invtau/60.0*pow(S-S_max,1.0/(1.0-lambda_1));
		else if(state == 4)
{
//if(h < 1e-6)	ans->ve[0] = 0.0;
			ans->ve[0] = c_2 * 2*sqrt(h*diam-h*h) * pow(h,1.5);
}
		else if(state == 5)
			ans->ve[0] = c_2 * diam * pow(h,1.5);
		else printf("Error in dam evaluation. Invalid state %u.\n",state);
	}
}


//Type 23
//Order of parameters: A_i,L_i,A_h,RC,v_h,v_r,k2,k3,invtau,orifice_area,H_spill,H_max,S_max,alpha,orifice_diam,c_1,c_2,L_spill
//The numbering is:	0   1   2  3   4   5   6  7   8          9	  10	  11   12     13      14        15  16   17
//Order of global_params: lambda_1,lambda_2,S_0,v_g
//The numbering is:         0        1       2   3
int dam_check3(VEC* y,VEC* global_params,VEC* params,QVSData* qvs,unsigned int dam)
{
	if(dam == 0)	return 0;

	double H_spill = params->ve[10];
	double H_max = params->ve[11];
	double S_max = params->ve[12];
	double alpha = params->ve[13];
	double diam = params->ve[14];
	double S = y->ve[1];
	double h = H_max * pow(S/S_max,alpha);

	if(h < 0.5*diam)	return 4;
	if(h < diam)		return 5;
	if(h <= H_spill)	return 1;
	if(h <= H_max)		return 2;
	return 3;
}


// Type 40 ************************************************

//Type 40
//Order of parameters: A_i,L_i,A_h,RC,v_h,v_r,k2,k3,invtau
//The numbering is:	0   1   2   3  4   5   6  7    8
//Order of global_params: lambda_1,lambda_2,S_0,v_g
//The numbering is:         0        1       2   3
void nodam_rain_hillslope_qsv(double t,VEC* y_i,VEC** y_p,unsigned short int numparents,VEC* global_params,double* forcing_values,QVSData* qvs,VEC* params,IVEC* iparams,int state,unsigned int** upstream,unsigned int* numupstream,VEC* ans)
{
	unsigned short int i;
	double lambda_1 = global_params->ve[0];
	double A_h = params->ve[2];
	double RC = params->ve[3];
	double k2 = params->ve[6];
	double k3 = params->ve[7];
	double invtau = params->ve[8];

	double S = y_i->ve[1];
	double Ss = y_i->ve[2];
	double Sg = y_i->ve[3];

	double qm = invtau*pow(S,1.0/(1.0-lambda_1));

	ans->ve[1] = k2*Ss + k3*Sg - qm;
	for(i=0;i<numparents;i++)	ans->ve[1] += y_p[i]->ve[0] * 60.0;

	ans->ve[2] = RC*forcing_values[0]*A_h*(.001/60.0) - k2*Ss;
	ans->ve[3] = (1.0-RC)*forcing_values[0]*A_h*(.001/60.0) - k3*Sg;
}


//Type 40
//Order of parameters: A_i,L_i,A_h,RC,v_h,v_r,k2,k3,invtau
//The numbering is:	0   1   2  3   4   5   6  7   8
//Order of global_params: lambda_1,lambda_2,S_0,v_g
//The numbering is:         0        1       2   3
void dam_rain_hillslope_qsv(double t,VEC* y_i,VEC** y_p,unsigned short int numparents,VEC* global_params,double* forcing_values,QVSData* qvs,VEC* params,IVEC* iparams,int state,unsigned int** upstream,unsigned int* numupstream,VEC* ans)
{
	unsigned short int i;
	double lambda_1 = global_params->ve[0];
	double A_h = params->ve[2];
	double RC = params->ve[3];
	double k2 = params->ve[6];
	double k3 = params->ve[7];
	double invtau = params->ve[8];

	double qm,q1,q2,S1,S2,S_max,q_max;
	double S = y_i->ve[1];
	double Ss = y_i->ve[2];
	double Sg = y_i->ve[3];

	if(state == -1)
		qm = invtau*pow(S,1.0/(1.0-lambda_1));
	else if(state == (int) qvs->n_values - 1)
	{
		S_max = qvs->points[qvs->n_values - 1][0];
		q_max = qvs->points[qvs->n_values - 1][1];
		qm = q_max*60.0 + invtau*pow(max(S-S_max,0.0),1.0/(1.0-lambda_1));
		//qm = q_max * 60.0;
	}
	else
	{
		q2 = qvs->points[state+1][1];
		q1 = qvs->points[state][1];
		S2 = qvs->points[state+1][0];
		S1 = qvs->points[state][0];
		qm = ((q2-q1)/(S2-S1) * (S-S1) + q1) * 60.0;
	}

	ans->ve[1] = k2*Ss + k3*Sg - qm;
	for(i=0;i<numparents;i++)	ans->ve[1] += y_p[i]->ve[0] * 60.0;

	ans->ve[2] = RC*forcing_values[0]*A_h*(.001/60.0) - k2*Ss;
	ans->ve[3] = (1.0-RC)*forcing_values[0]*A_h*(.001/60.0) - k3*Sg;
}


//Type 40
//Order of parameters: A_i,L_i,A_h,RC,v_h,v_r,k2,k3,invtau
//The numbering is:	0   1   2  3   4   5   6  7   8
//Order of global_params: lambda_1,lambda_2,S_0,v_g
//The numbering is:         0        1       2   3
void dam_q_qvs(VEC* y,VEC* global_params,VEC* params,QVSData* qvs,int state,VEC* ans)
{
	double S_max,q_max;
	double lambda_1 = global_params->ve[0];
	double invtau = params->ve[8];
	double S = (y->ve[1] < 0.0) ? 0.0 : y->ve[1];
	double q1,q2,S1,S2;

	if(state == -1)
		ans->ve[0] = invtau/60.0*pow(S,1.0/(1.0-lambda_1));
	else if(state == (int) qvs->n_values - 1)
	{
		S_max = qvs->points[qvs->n_values - 1][0];
		q_max = qvs->points[qvs->n_values - 1][1];
		ans->ve[0] = q_max + invtau/60.0*pow(max(S-S_max,0.0),1.0/(1.0-lambda_1));
		//ans->ve[0] = q_max;
	}
	else
	{
		q2 = qvs->points[state+1][1];
		q1 = qvs->points[state][1];
		S2 = qvs->points[state+1][0];
		S1 = qvs->points[state][0];
		ans->ve[0] = (q2-q1)/(S2-S1) * (S-S1) + q1;
	}
}


//Type 40 / 261 / 262
int dam_check_qvs(VEC* y,VEC* global_params,VEC* params,QVSData* qvs,unsigned int dam)
{
	unsigned int i,iterations,max,min,k;
	double S = y->ve[1];

	if(dam == 0)	return -1;

	iterations = qvs->n_values - 1;
	for(i=0;i<iterations;i++)
	{
		if(qvs->points[i][0] <= S && S < qvs->points[i+1][0])
			return i;
	}

	return i;
}


// Rodica's Models **************************

//Type 0
//Function for simple river system
//Calculates the flow using simple parameters, using only the flow q.
//Order of parameters: L_i,A_h,A_i,h_b,h_H,max_inf_rate,K_sat,S_h,eta,b_H,c_H,d_H,invtau,epsilon,c_1,c_2,c_3,c_4,c_5,c_6
//The numbering is:     0   1   2   3   4       5         6    7   8   9   10  11  12    13      14  15  16  17  18  19
//Order of global_params: v_r,lambda_1,lambda_2,Q_r,A_r,RC
//The numbering is:        0      1        2     3   4   5
//This uses the units and functions from September 18, 2011 document
//y_i[0] = q
//void simple_river(double t,VEC* y_i,VEC** y_p,unsigned short int numparents,VEC* global_params,double* forcing_values,VEC* params,int state,VEC* ans)
void simple_river(double t,VEC* y_i,VEC** y_p,unsigned short int numparents,VEC* global_params,double* forcing_values,QVSData* qvs,VEC* params,IVEC* iparams,int state,unsigned int** upstream,unsigned int* numupstream,VEC* ans)
{
	unsigned short int i;

	double q = y_i->ve[0];
	double lambda_1 = global_params->ve[1];
	double invtau = params->ve[12];

	//Flux equation (y_i[0])
	ans->ve[0] = -q;
	for(i=0;i<numparents;i++)
		ans->ve[0] += y_p[i]->ve[0];
	ans->ve[0] = invtau * pow(q,lambda_1) * ans->ve[0];
}

//Type 0
//This assumes there is only 1 variable being passed link to link.
//double NormJx_simple_river(VEC* y_i,VEC* global_params,VEC* params)
void Jx_simple_river(VEC* y_i,VEC* global_params,VEC* params,VEC* ans)
{
	ans->ve[0] = params->ve[12] * pow(y_i->ve[0],global_params->ve[1]);
	//return params->ve[12] * pow(y_i->ve[0],global_params->ve[1]);
}

//Type 0
//Jacobian of the simple river system******************
//y_i has q_i (y_i[0])
//Order of parameters: L_i,A_h,A_i,h_b,h_H,max_inf_rate,K_sat,S_h,eta,b_H,c_H,d_H,invtau,epsilon,c_1,c_2,c_3,c_4,c_5,c_6
//The numbering is:     0   1   2   3   4       5         6    7   8   9   10  11   12     13     14  15  16  17  18  19
//Order of global_params: v_r,lambda_1,lambda_2,Q_r,A_r,RC
//The numbering is:        0      1        2     3   4   5
//This uses the units and functions from September 18, 2011 document
void Jsimple(double t,VEC* y_i,VEC** y_p,unsigned short int numparents,VEC* global_params,double* forcing_values,VEC* params,MAT* ans)
{
	unsigned short int i;
	double q = y_i->ve[0];
	double lambda_1 = global_params->ve[1];
	double invtau = params->ve[12];

	ans->me[0][0] = 0.0;
	for(i=0;i<numparents;i++)
		ans->me[0][0] += y_p[i]->ve[0];
	ans->me[0][0] = -(lambda_1 + 1.0) * invtau * pow(q,lambda_1) + lambda_1 * invtau * pow(q,lambda_1-1.0) * ans->me[0][0];
}

//Type 1/3
//Function for the simple river system with rainfall
//Calculates the flow using simple parameters, using only the flow q.
//Order of parameters: L_i,A_h,A_i,h_b,h_H,max_inf_rate,K_sat,S_h,eta,b_H,c_H,d_H,invtau,epsilon,c_1,c_2,c_3,c_4,c_5,c_6
//The numbering is:     0   1   2   3   4       5         6    7   8   9   10  11  12    13      14  15  16  17  18  19
//Order of global_params: v_r,lambda_1,lambda_2,Q_r,A_r,RC
//The numbering is:        0      1        2     3   4   5
//This uses the units and functions from September 18, 2011 document
//y_i[0] = q, y_i[1] = s
void river_rainfall(double t,VEC* y_i,VEC** y_p,unsigned short int numparents,VEC* global_params,double* forcing_values,QVSData* qvs,VEC* params,IVEC* iparams,int state,unsigned int** upstream,unsigned int* numupstream,VEC* ans)
{
	unsigned short int i;

	double q = y_i->ve[0];
	double s = y_i->ve[1];

	double invtau = params->ve[12];
	double c_1 = params->ve[14];
	double c_3 = params->ve[16];
	double c_4 = params->ve[17];
	double lambda_1 = global_params->ve[1];

	double s53 = pow(s,5.0/3.0);

	//Flux equation (y_i[0])
	ans->ve[0] = -q + c_1 * s53;
	for(i=0;i<numparents;i++)
		ans->ve[0] += y_p[i]->ve[0];
	ans->ve[0] = invtau * pow(q,lambda_1) * ans->ve[0];

	//Ponded water equation (y_i[1])
	ans->ve[1] = c_3 * forcing_values[0] - c_4 * s53;
}

//Type 2
//Function for simple river system with hillslope runoff
//Calculates the flow using simple parameters, using only the flow q.
//Order of parameters: L_i,A_h,A_i,h_b,h_H,max_inf_rate,K_sat,S_h,eta,b_H,c_H,d_H,invtau,epsilon,c_1,c_2,c_3,c_4,c_5,c_6
//The numbering is:     0   1   2   3   4       5         6    7   8   9   10  11  12    13      14  15  16  17  18  19
//Order of global_params: v_r,lambda_1,lambda_2,Q_r,A_r,RC
//The numbering is:        0      1        2     3   4   5
//This uses the units and functions from September 18, 2011 document
//y_i[0] = q, y_i[1] = s
void simple_hillslope(double t,VEC* y_i,VEC** y_p,unsigned short int numparents,VEC* global_params,double* forcing_values,QVSData* qvs,VEC* params,IVEC* iparams,int state,unsigned int** upstream,unsigned int* numupstream,VEC* ans)
{
	unsigned short int i;

	double q = y_i->ve[0];
	double s = y_i->ve[1];

	double invtau = params->ve[12];
	double c_1 = params->ve[14];
	double c_4 = params->ve[17];
	double lambda_1 = global_params->ve[1];

	double s53 = pow(s,5.0/3.0);

	//Flux equation (y_i[0])
	ans->ve[0] = -q + c_1 * s53;
	for(i=0;i<numparents;i++)
		ans->ve[0] += y_p[i]->ve[0];
	ans->ve[0] = invtau * pow(q,lambda_1) * ans->ve[0];

	//Ponded water equation (y_i[1])
	ans->ve[1] = - c_4 * s53;
}

//Type 4
//Function for simple river system with soil properties******************
//Calculates the flow using the soil parameters.
//y_i has q_i (y_i[0]), s (y_i[1]), a (y_i[2]), v (y_i[3])
//Order of parameters: L_i,A_h,A_i,h_b,h_H,max_inf_rate,K_sat,S_h,eta,b_H,c_H,d_H,invtau,epsilon,c_1,c_2,c_3,c_4,c_5,c_6
//The numbering is:     0   1   2   3   4       5         6    7   8   9   10  11  12    13      14  15  16  17  18  19
//Order of global_params: v_r,lambda_1,lambda_2,Q_r,A_r,RC
//The numbering is:        0      1        2     3   4   5
//This uses the units and functions from September 18, 2011 document  
void simple_soil(double t,VEC* y_i,VEC** y_p,unsigned short int numparents,VEC* global_params,double* forcing_values,QVSData* qvs,VEC* params,IVEC* iparams,int state,unsigned int** upstream,unsigned int* numupstream,VEC* ans)
{
	unsigned short int i;

	double q = y_i->ve[0];
	double s = y_i->ve[1];
	double a = y_i->ve[2];
	double v = y_i->ve[3];

	double invtau = params->ve[12];
	double epsilon = params->ve[13];
	double c_1 = params->ve[14];
	double c_4 = params->ve[17];
	double c_5 = params->ve[18];
	double c_6 = params->ve[19];
	double lambda_1 = global_params->ve[1];

	double s53 = pow(s,5.0/3.0);
	double exp_term,big_frac;

	if(a < 1.0)
	{
		exp_term = c_6 * (1.0 - a) * v * exp(-(1.0-a-v)/(1.0-a));
		//big_frac = c_5 * s * epsilon*epsilon * (1.0-a) * (1.0-a-v)*(1.0-a-v)/((s*(1.0-a) + epsilon*(1.0 - a - v))*(s*(1.0-a) + 1.0 - a - v));
		big_frac = c_5 * s * sq(epsilon) * (1.0-a) * sq(1.0-a-v) / sq( s*(1.0-a) + epsilon*(1.0 - a - v) );
	}
	else
	{
		exp_term = 0.0;
		big_frac = 0.0;
	}

	//Flux equation (y_i[0])
	ans->ve[0] = -q + c_1 * s53;
	for(i=0;i<numparents;i++)
		ans->ve[0] += y_p[i]->ve[0];
	ans->ve[0] = invtau * pow(q,lambda_1)*ans->ve[0];

	//Ponded water equation (y_i[1])
	ans->ve[1] = -c_4 * s53 - big_frac + c_6 * a * a;

	//Impermeable area (y_i[2])
	//ans->ve[2] = -c_6 * a * a + c_6 * (1.0 - a) + exp_term;
	ans->ve[2] = (-c_6 * a * a + exp_term) / epsilon;

	//Soil Volumetric Water Content (y_i[3])
	//ans->ve[3] = big_frac - exp_term;
	ans->ve[3] = (big_frac - exp_term) / epsilon;
}

//Type 5
//Function for simple river system with soil properties and rainfall******************
//Calculates the flow using the soil parameters.
//y_i has q_i (y_i[0]), s (y_i[1]), a (y_i[2]), v (y_i[3])
//Order of parameters: L_i,A_h,A_i,h_b,h_H,max_inf_rate,K_sat,S_h,eta,b_H,c_H,d_H,invtau,epsilon,c_1,c_2,c_3,c_4,c_5,c_6
//The numbering is:     0   1   2   3   4       5         6    7   8   9   10  11   12     13     14  15  16  17  18  19
//Order of global_params: v_r,lambda_1,lambda_2,Q_r,A_r,RC
//The numbering is:        0      1        2     3   4   5
//This uses the units and functions from September 18, 2011 document
void soil_rainfall(double t,VEC* y_i,VEC** y_p,unsigned short int numparents,VEC* global_params,double* forcing_values,QVSData* qvs,VEC* params,IVEC* iparams,int state,unsigned int** upstream,unsigned int* numupstream,VEC* ans)
{
	unsigned short int i;

	double q = y_i->ve[0];
	double s = y_i->ve[1];
	double a = y_i->ve[2];
	double v = y_i->ve[3];

	double invtau = params->ve[12];
	double epsilon = params->ve[13];
	double c_1 = params->ve[14];
	double c_3 = params->ve[16];
	double c_4 = params->ve[17];
	double c_5 = params->ve[18];
	double c_6 = params->ve[19];
	double lambda_1 = global_params->ve[1];

	double s53 = pow(s,5.0/3.0);
	double exp_term,big_frac;

	if(a < 1.0)
	{
		exp_term = c_6 * (1.0 - a) * v * exp(-(1.0-a-v)/(1.0-a));
		//big_frac = c_5 * s * epsilon*epsilon * (1.0-a) * (1.0-a-v)*(1.0-a-v)/((s*(1.0-a) + epsilon*(1.0 - a - v))*(s*(1.0-a) + 1.0 - a - v));
		big_frac = c_5 * s * sq(epsilon) * (1.0-a) * sq(1.0-a-v) / sq( s*(1.0-a) + epsilon*(1.0 - a - v) );
	}
	else
	{
		exp_term = 0.0;
		big_frac = 0.0;
	}

	//Flux equation (y_i[0])
	ans->ve[0] = -q + c_1 * s53;
	for(i=0;i<numparents;i++)
		ans->ve[0] += y_p[i]->ve[0];
	ans->ve[0] = invtau * pow(q,lambda_1)*ans->ve[0];

	//Ponded water equation (y_i[1])
	ans->ve[1] = c_3 * forcing_values[0] - c_4 * s53 - big_frac + c_6 * a * a;

	//Impermeable area (y_i[2])
	//ans->ve[2] = -c_6 * a * a + c_6 * (1.0 - a) + exp_term;
	ans->ve[2] = (-c_6 * a * a + exp_term) / epsilon;

	//Soil Volumetric Water Content (y_i[3])
	//ans->ve[3] = big_frac - exp_term;
	ans->ve[3] = (big_frac - exp_term) / epsilon;
}

//Type 6
//RHS function for discharge, surface ponding, saturated zone, unsaturated zone, and vegetation
//Includes rainfall
//This uses units and functions from Rodica's 1st (published) paper
//Order of parameters: L_i,A_h,A_i,h_b,K_sat,K_sp,d_0,d_1,d_2,invbeta,alpha_N,alpha_soil,a_res,v_res,invtau,gamma,c_1,c_2,c_3,c_4
//The numbering is:     0   1   2   3   4     5    6   7   8     9      10       11       12    13     14    15   16  17  18  19
//Order of global_params: v_r,lambda_1,lambda_2,Q_r,A_r
//The numbering is:        0      1        2     3   4
void qsav_rainfall(double t,VEC* y_i,VEC** y_p,unsigned short int numparents,VEC* global_params,double* forcing_values,QVSData* qvs,VEC* params,IVEC* iparams,int state,unsigned int** upstream,unsigned int* numupstream,VEC* ans)
{
	unsigned short int i;

	//State variables
	double q = y_i->ve[0];
	double s_p = y_i->ve[1];
	double a = y_i->ve[2];
	double v = y_i->ve[3];

	//Parameters and constants
	double h_b = params->ve[3];
	double d_0 = params->ve[6];
	double d_1 = params->ve[7];
	double d_2 = params->ve[8];
	double invbeta = params->ve[9];
	double a_res = params->ve[12];
	double v_res = params->ve[13];
	double invtau = params->ve[14];
	double gamma = params->ve[15];
	double c_1 = params->ve[16];
	double c_2 = params->ve[17];
	double c_3 = params->ve[18];
	double c_4 = params->ve[19];
	double lambda_1 = global_params->ve[1];

	//Calculations
	double a_adj = a - a_res;
	double v_adj = v - v_res;
	double av_adjsum = a_adj + v_adj;
	double polyterm = d_0 * v_adj + d_1 * v_adj * sq(a_adj) + d_2 * sq(a_adj);
	double expterm = exp(a_adj/c_4);

	//Evaporation
/*
	double c_evap = 1.0/60.0;
	if(t < 8000.0)		c_evap *= 0.002;
	else if(t < 15000.0)	c_evap *= 0.0005;
	else if(t < 30000.0)	c_evap *= 0.0001;
	else			c_evap = 0.0;

	double evap;
	if(t < 5000.0)		evap = 0.62;
	else if(t < 18000.0)	evap = 0.48;
	else if(t < 25000.0)	evap = 0.5;
	else if(t < 30000.0)	evap = 0.48;
	else if(t < 35000.0)	evap = 0.51;
	else if(t < 40000.0)	evap = 0.44;
	else			evap = 0.3;
*/
	double c_evap = 0.0;
	double evap = 0.0;

	//Flux equation
	ans->ve[0] = -q + gamma * c_1 * s_p * av_adjsum + gamma * c_2 * a_adj * expterm;
	for(i=0;i<numparents;i++)
		ans->ve[0] += y_p[i]->ve[0];
	ans->ve[0] = invtau * pow(q,lambda_1) * ans->ve[0];

	//Ponded water equation
	//ans->ve[1] = (1.0 - evap) * c_3 * forcing_values[0] - c_1 * s_p * av_adjsum - c_1 * s_p * (h_b - av_adjsum);
	ans->ve[1] = (1.0 - evap) * c_3 * forcing_values[0] - c_1 * h_b * s_p;

	//Saturated zone equation (a)
	ans->ve[2] = invbeta * (polyterm - c_2 * a_adj * expterm - c_evap * a_adj);

	//Unsaturated zone equation (v)
	ans->ve[3] = invbeta * ( c_1 * s_p * (h_b - av_adjsum) - polyterm );
}



//Type 15
//Function for the simple river system with rainfall
//Calculates the flow using simple parameters, using only the flow q.
//Order of parameters: L_i,A_h,A_i,h_b,h_H,max_inf_rate,K_sat,S_h,eta,b_H,c_H,d_H,invtau,epsilon,c_1,c_2,c_3,c_4,c_5,c_6
//The numbering is:     0   1   2   3   4       5         6    7   8   9   10  11  12    13      14  15  16  17  18  19
//Order of global_params: v_r,lambda_1,lambda_2,v_h,A_r,RC
//The numbering is:        0      1        2     3   4   5
//This uses the units and functions from September 18, 2011 document
//y_i[0] = q, y_i[1] = s_p
void river_rainfall_adjusted(double t,VEC* y_i,VEC** y_p,unsigned short int numparents,VEC* global_params,double* forcing_values,QVSData* qvs,VEC* params,IVEC* iparams,int state,unsigned int** upstream,unsigned int* numupstream,VEC* ans)
{
	unsigned short int i;

	double q = y_i->ve[0];
	double s_p = y_i->ve[1];

	double L = params->ve[0];
	double invtau = params->ve[12];
	double c_1 = params->ve[14];
	double c_3 = params->ve[16];
	double c_4 = params->ve[17];
	double lambda_1 = global_params->ve[1];

	//double q_pl = pow(0.5*L*s_p/(L+s_p),2.0/3.0);
	double q_pl = s_p;

	//Flux equation (y_i[0])
	ans->ve[0] = -q + c_1 * q_pl;
	for(i=0;i<numparents;i++)
		ans->ve[0] += y_p[i]->ve[0];
	ans->ve[0] = invtau * pow(q,lambda_1) * ans->ve[0];

	//Ponded water equation (y_i[1])
	ans->ve[1] = c_3 * forcing_values[0] - c_4 * q_pl;
	//ans->ve[1] = c_3 * (forcing_values[0] + 100.0*sin(t/5.0)) - c_4 * q_pl;
}

//Type 315
//Function for simple river system with data assimilation.
//Calculates the flow using simple parameters, using only the flow q.
//Order of parameters: L_i,A_h,A_i,h_b,h_H,max_inf_rate,K_sat,S_h,eta,b_H,c_H,d_H,invtau,epsilon,c_1,c_2,c_3,c_4,c_5,c_6
//The numbering is:     0   1   2   3   4       5         6    7   8   9   10  11  12    13      14  15  16  17  18  19
//Order of global_params: v_r,lambda_1,lambda_2,v_h,A_r,RC
//The numbering is:        0      1        2     3   4   5
//This uses the units and functions from September 18, 2011 document
//y_i[0] = q, y_i[1] = s, followed by N entries for the variational equation
void assim_river_rainfall_adjusted(double t,VEC* y_i,VEC** y_p,unsigned short int numparents,VEC* global_params,double* forcing_values,QVSData* qvs,VEC* params,IVEC* iparams,int state,unsigned int** upstream,unsigned int* numupstream,VEC* ans)
{
	unsigned int i,j;
	unsigned int dim = ans->dim;
	unsigned int offset = 2;		//!!!! This needs to be num_dense, but without variational eqs !!!!
	unsigned int parent_offset;
	unsigned int problem_dim = 2;
	unsigned int all_states = (dim-offset)/problem_dim;
	unsigned int loc = iparams->ve[0];
	double inflow = 0.0;

	double q = y_i->ve[0];
	double s_p = y_i->ve[1];

	double L = params->ve[0];
	double invtau = params->ve[12];
	double c_1 = params->ve[14];
	double c_3 = params->ve[16];
	double c_4 = params->ve[17];
	double lambda_1 = global_params->ve[1];

	double q_to_lambda_1 = pow(q,lambda_1);
	//double q_to_lambda_1_m1 = (q > 1e-12) ? q_to_lambda_1 / q : 1e-12;
	double q_to_lambda_1_m1 = (q > 1e-12) ? q_to_lambda_1 / q : pow(1e-12,lambda_1 - 1.0);
	double deriv_qpl = 1.0;

	double q_pl = s_p;

	//Flux equation (y_i[0])
	ans->ve[0] = -q + c_1 * q_pl;
	for(i=0;i<numparents;i++)
		inflow += y_p[i]->ve[0];
	ans->ve[0] = invtau * q_to_lambda_1 * (inflow + ans->ve[0]);

	//Ponded water equation (y_i[1])
	ans->ve[1] = c_3 * forcing_values[0] - c_4 * q_pl;
	//ans->ve[1] = c_3 * ( max(forcing_values[0] + 20.0*sin(t/5.0),0.0)) - c_4 * q_pl;

	//!!!! Pull if statements out of loops (should just need two cases total) !!!!
	//!!!! A lot of terms get repeated !!!!

	//Eqs for variational equations
	for(i=offset;i<dim;i++)	ans->ve[i] = 0.0;

	//s variable from this link
	ans->ve[offset] = -c_4*deriv_qpl*y_i->ve[offset];

	//q variables from this link
//	if(lambda_1 > 1e-12 && (inflow) > 1e-12)
		ans->ve[offset + 1] = (lambda_1 * invtau * q_to_lambda_1_m1 * (inflow + c_1*s_p) - (lambda_1 + 1) * invtau * q_to_lambda_1) * y_i->ve[offset + 1];
//	else
//		ans->ve[offset + 1] = -(lambda_1 + 1.0) * invtau * q_to_lambda_1 * y_i->ve[offset + 1];

//	if(lambda_1 > 1e-12 && (inflow) > 1e-12)
		ans->ve[offset + 2] = (lambda_1 * invtau * q_to_lambda_1_m1 * (inflow + c_1*s_p) - (lambda_1 + 1) * invtau * q_to_lambda_1) * y_i->ve[offset + 2] + invtau*c_1*q_to_lambda_1*deriv_qpl * y_i->ve[offset];
//	else
//		ans->ve[offset + 2] = -(lambda_1 + 1.0) * invtau * q_to_lambda_1 * y_i->ve[offset + 2] + invtau*c_1*deriv_qpl*y_i->ve[offset];

	//Adjust offset
	offset += 3;

	//Variables from parents
	for(i=0;i<numparents;i++)
	{
		parent_offset = 1 + problem_dim;

		for(j=0;j<numupstream[i];j++)
		{
			ans->ve[offset] = invtau * q_to_lambda_1 * y_p[i]->ve[parent_offset];
//			if(lambda_1 > 1e-12 && (inflow) > 1e-12)
				ans->ve[offset] += (lambda_1 * invtau * q_to_lambda_1_m1 * (inflow + c_1*s_p) - (lambda_1 + 1) * invtau * q_to_lambda_1) * y_i->ve[offset];
//			else
//				ans->ve[offset] += -(lambda_1 + 1.0) * invtau * q_to_lambda_1 * y_i->ve[offset];

			ans->ve[offset + 1] = invtau * q_to_lambda_1 * y_p[i]->ve[parent_offset + 1];
//			if(lambda_1 > 1e-12 && (inflow) > 1e-12)
				ans->ve[offset + 1] += (lambda_1 * invtau * q_to_lambda_1_m1 * (inflow + c_1*s_p) - (lambda_1 + 1) * invtau * q_to_lambda_1) * y_i->ve[offset + 1];
//			else
//				ans->ve[offset + 1] += -(lambda_1 + 1.0) * invtau * q_to_lambda_1 * y_i->ve[offset + 1];

			offset += 2;
			parent_offset += 2;
		}
	}
}



//Type 30
//Function for river system with soil properties and rainfall******************
//This uses Luciana's model
//Calculates the flow using the soil parameters.
//y_i has q (y_i[0]), s_p (y_i[1]), h_w (y_i[2]), theta (y_i[3])
//Order of parameters: L_h,A_h,A_up,H_b,H_h,max_inf_rate,K_SAT,S_H,n_vh, b, c, d,K_Q,V_T,c_1,c_2,c_3,c_4,c_5,c_6,c_7
//The numbering is:     0   1   2    3   4       5         6    7   8    9 10 11  12  13  14  15  16  17  18  19  20
//Order of global_params: v_0,lambda_1,lambda_2,Q_r,A_r,K_T,C_r,e_pot
//The numbering is:        0      1        2     3   4   5   6   7
/*
void lcuencas_soilrain(double t,VEC* y_i,VEC** y_p,unsigned short int numparents,VEC* global_params,double* forcing_values,QVSData* qvs,VEC* params,IVEC* iparams,int state,unsigned int** upstream,unsigned int* numupstream,VEC* ans)
{
	unsigned int i;

	//States
	double q = y_i->ve[0];		//[m^3/s]
	double s_p = y_i->ve[1];	//[mm]
	double h_w = y_i->ve[2];	//[m]
	double theta = y_i->ve[3];	//[%]

	//Global params
	//double v_0 = global_params->ve[0];
	double lambda_1 = global_params->ve[1];
	//double lambda_2 = global_params->ve[2];
	//double A_r = global_params->ve[4];
	double K_T = global_params->ve[5];
	double e_pot = global_params->ve[7];

	//Local params
	double A_h = params->ve[1];
	double H_b = params->ve[3];
	double H_h = params->ve[4];
	double maxinf = params->ve[5];
	//double K_SAT = params->ve[6];
	double K_SAT = (H_h > 1e-12) ? params->ve[6] : 0.0;
	double b = params->ve[9];
	double c = params->ve[10];
	double d = params->ve[11];
	double K_Q = params->ve[12];
	double V_T = params->ve[13];
	double c_1 = params->ve[14];
	double c_2 = params->ve[15];
	double c_3 = params->ve[16];
	double c_4 = params->ve[17];
	double c_5 = params->ve[18];
	double c_6 = params->ve[19];
	double c_7 = params->ve[20];

	//Other
	//double H_T = H_b + H_h;
	double h_rel = h_w - H_b * 1e-3;
	//double H_relmax = H_T - H_b;
	double H_relmax = (H_h > 1e-12) ? H_h : 1e-6;
	double hoverH = h_rel / H_relmax;
	double a_IoverA_H = b * hoverH + c * sq(hoverH) + d * sq(hoverH) * hoverH;
	if(a_IoverA_H < 0.0)	a_IoverA_H = 1e-6;
	if(a_IoverA_H > 0.9999)	a_IoverA_H = 0.9999;
	double a_PoverA_H = 1.0 - a_IoverA_H;
	double v_H = c_1 * pow(s_p * 1e-3,2.0/3.0);
	if(v_H > 350.0)	v_H = 350.0;
	double v_ssat = V_T * a_IoverA_H;
	double v_sunsat = V_T - v_ssat;
	double d_soil = (1.0 - theta) * H_b;
	double RC = (s_p <= 0.0) ? 0.0 : s_p*(s_p + 2*d_soil) / sq(s_p + d_soil);
	double deriv_a_I = c_7 * (b + 2*c*hoverH + 3*d*sq(hoverH));
	if(deriv_a_I <= (1e-5*1e6)) deriv_a_I = (100.0*1e6);		//Should use this
	//if(deriv_a_I <= (1e-0)) deriv_a_I = 1e-0;
	double K_UNSAT = K_SAT * exp(5.0 * (theta - 1.0));
	double infrate = (s_p > maxinf) ? infrate = maxinf : s_p * K_T;
	//double infrate = s_p*K_T;
	//if(infrate > maxinf)	infrate = maxinf;


	//Evaporation
	//double C_p = (e_pot > 0.0) ? s_p / e_pot : 0.0;
	double D_unsat = 1e-3 * v_sunsat * theta / A_h;
	double D_sat = 1e-3 * v_ssat / A_h;

	double C_p,C_unsat,C_sat,C_T;
	if(e_pot > 0.0)
	{
		C_p = s_p / e_pot;
		C_unsat = D_unsat / e_pot;
		C_sat = D_sat / e_pot;
		C_T = C_p + C_unsat + C_sat;
	}
	else
	{
		C_p = 0.0;
		C_unsat = 0.0;
		C_sat = 0.0;
		C_T = 0.0;
	}


	//double C_unsat = D_unsat / H_b * a_PoverA_H;
	//double C_sat = D_sat/H_b * a_IoverA_H;
	//double C_s = 0.0;
	//C_T = C_s + C_unsat + C_sat;

	double Corr_evap = (C_T > 1.0) ? 1.0/C_T : 1.0;
	double e_p = Corr_evap * C_p * e_pot;
	double e_sat = Corr_evap * C_sat * e_pot;
	double e_unsat = Corr_evap * C_unsat * e_pot;


	//Fluxes
	double q_pl = c_2 * (a_PoverA_H * RC + a_IoverA_H) * v_H * s_p;
	double q_sl = c_3 * v_ssat * a_IoverA_H;
	double snow_M = 0.8 * 0.0;	//Assume no snow for now
	double q_pu = a_PoverA_H * infrate * (1.0 - RC);
	double q_us = 1e3 * K_UNSAT * theta * a_PoverA_H;

	//!!!! This might be cleaner if H_b is converted into m !!!!

	//Channel transport
	ans->ve[0] = -q;
	for(i=0;i<numparents;i++)	ans->ve[0] += y_p[i]->ve[0];
	ans->ve[0] = K_Q * pow(q,lambda_1) * (ans->ve[0] + c_4 * (q_pl + q_sl + snow_M));

	//Land surface storage
	ans->ve[1] = (1.0/60.0) * (forcing_values[0] - q_pl - q_pu - e_p);

	//Water table
	ans->ve[2] = (H_b < 1e-6) ? 0.0 : c_5 / deriv_a_I * (q_us - q_sl - e_sat);

	//Soil volumetric water content
	//ans->ve[3] = (v_sunsat < 1e-6) ? 0.0 : c_6 / v_sunsat * (q_pu - q_us - e_unsat + theta*(q_us - q_sl - e_sat));
	//ans->ve[3] = (v_sunsat < 1e-12) ? -100.0 : c_6 / v_sunsat * (q_pu - q_us - e_unsat + theta*(q_us - q_sl - e_sat));
	ans->ve[3] = c_6 / v_sunsat * (q_pu - q_us - e_unsat + theta*(q_us - q_sl - e_sat));
}
*/

void lcuencas_soilrain(double t,VEC* y_i,VEC** y_p,unsigned short int numparents,VEC* global_params,double* forcing_values,QVSData* qvs,VEC* params,IVEC* iparams,int state,unsigned int** upstream,unsigned int* numupstream,VEC* ans)
{
	unsigned int i;

	//States
	double q = y_i->ve[0];
	double s_p = y_i->ve[1];
	double h_w = y_i->ve[2];
	double theta = y_i->ve[3];

	//Global params
	//double v_0 = global_params->ve[0];
	double lambda_1 = global_params->ve[1];
	//double lambda_2 = global_params->ve[2];
	//double A_r = global_params->ve[4];
	double K_T = global_params->ve[5];
	double e_pot = forcing_values[1];
	//double e_pot = global_params->ve[7];

	//Local params
	//double L_h = params->ve[0];
	double A_h = params->ve[1];
	//double A_up = params->ve[2];
	double H_b = params->ve[3];
	double H_h = params->ve[4];
	double maxinf = params->ve[5];
	double K_SAT = params->ve[6];
	//double S_H = params->ve[7];
	//double n_vh = params->ve[8];
	double b = params->ve[9];
	double c = params->ve[10];
	double d = params->ve[11];
	double K_Q = params->ve[12];
	double V_T = params->ve[13];
	double c_1 = params->ve[14];
	double c_2 = params->ve[15];
	double c_3 = params->ve[16];
	double c_4 = params->ve[17];
	double c_5 = params->ve[18];
	double c_6 = params->ve[19];
	double c_7 = params->ve[20];

	//Other
	//double H_T = H_b + H_h;
	double h_rel = h_w - H_b * 1e-3;
//if(h_rel <= 0)	h_rel = 0.01;
	double H_relmax = (H_h > 1e-12) ? H_h : 1e-6;
	double hoverH = h_rel / H_relmax;
	double a_IoverA_H = b * hoverH + c * sq(hoverH) + d * sq(hoverH) * hoverH;
	if(a_IoverA_H < 0.0)	a_IoverA_H = 1e-6;
//if(a_IoverA_H <= 0)	a_IoverA_H = 1e-8 / A_h;
	double a_PoverA_H = 1.0 - a_IoverA_H;
	double v_H = c_1 * pow(s_p * 1e-3,2.0/3.0);
	if(v_H > 350.0)	v_H = 350.0;
	double v_ssat = V_T * a_IoverA_H;
	double v_sunsat = V_T - v_ssat;
	double d_soil = (1.0 - theta) * H_b;
	double RC = (s_p <= 0.0) ? 0.0 : s_p*(s_p + 2*d_soil) / sq(s_p + d_soil);
	double deriv_a_I = c_7 * (b + 2*c*hoverH + 3*d*sq(hoverH));
	//if(deriv_a_I <= 0.0)	deriv_a_I = 100.0;
	if(deriv_a_I <= 1e-5) deriv_a_I = 100.0;
	double K_UNSAT = K_SAT * exp(5.0 * (theta - 1.0));
	double infrate = (s_p > maxinf) ? infrate = maxinf : s_p * K_T;

	//Evaporation
	//double e_pot = 0.0;	//!!!! This should be an input function !!!!
	double C_p = (e_pot > 0.0) ? s_p / e_pot : 0.0;
	double D_unsat = 1e-3 * v_sunsat * theta / A_h;
	double D_sat = 1e-3 * v_ssat / A_h;
	double C_unsat = D_unsat / H_b * a_PoverA_H;
	double C_sat = D_sat/H_b * a_IoverA_H;
	//double C_s = 0.0;		//!!!! Don't know what to use here !!!!
	//double C_T = C_s + C_unsat + C_sat;
	double C_T = C_p + C_unsat + C_sat;
	double Corr_evap = (C_T > 1.0) ? 1.0/C_T : 1.0;
	double e_p = Corr_evap * C_p * e_pot;
	double e_sat = Corr_evap * C_sat * e_pot;
	double e_unsat = Corr_evap * C_unsat * e_pot;

	//Fluxes
	double q_pl = c_2 * (a_PoverA_H * RC + a_IoverA_H) * v_H * s_p;
	double q_sl = c_3 * v_ssat * a_IoverA_H;
	double snow_M = 0.8 * 0.0;	//Assume no snow for now
	double q_pu = a_PoverA_H * infrate * (1.0 - RC);
	double q_us = 1e3 * K_UNSAT * theta * a_PoverA_H;

	//!!!! This might be cleaner if H_b is converted into m !!!!

	//Channel transport
	ans->ve[0] = -q;
	for(i=0;i<numparents;i++)	ans->ve[0] += y_p[i]->ve[0];
	ans->ve[0] = K_Q * pow(q,lambda_1) * (ans->ve[0] + c_4 * (q_pl + q_sl + snow_M));

	//Land surface storage
	ans->ve[1] = (1.0/60.0) * (forcing_values[0] - q_pl - q_pu - e_p);

	//Water table
	ans->ve[2] = (H_b < 1e-6) ? 0.0 : c_5 / deriv_a_I * (q_us - q_sl - e_sat);

	//Soil volumetric water content
	ans->ve[3] = (v_sunsat < 1e-6) ? 0.0 : c_6 / v_sunsat * (q_pu - q_us - e_unsat + theta*(q_us - q_sl - e_sat));
/*
dump->ve[0] = e_p + e_sat + e_unsat;
dump->ve[1] = q_pu;
dump->ve[2] = q_pl;
dump->ve[3] = q_sl;
*/
}


//Type 4
//Jacobian of the simple river system with soil properties function******************
//y_i has q_i (y_i[0]), s (y_i[1]), a (y_i[2]), v (y_i[3])
//Order of parameters: L_i,A_h,A_i,h_b,h_H,max_inf_rate,K_sat,S_h,eta,b_H,c_H,d_H,invtau,epsilon,c_1,c_2,c_3,c_4,c_5,c_6
//The numbering is:     0   1   2   3   4       5         6    7   8   9   10  11   12     13     14  15  16  17  18  19
//Order of global_params: v_r,lambda_1,lambda_2,Q_r,A_r,RC
//The numbering is:        0      1        2     3   4   5
//This uses the units and functions from September 18, 2011 document
void Jsimple_soil(double t,VEC* y_i,VEC** y_p,unsigned short int numparents,VEC* global_params,double* forcing_values,VEC* params,MAT* ans)
{
	unsigned short int i;

	double q = y_i->ve[0];
	double s = y_i->ve[1];
	double a = y_i->ve[2];
	double v = y_i->ve[3];
	double** result = ans->me;

	double invtau = params->ve[12];
	double epsilon = params->ve[13];
	double c_1 = params->ve[14];
	double c_4 = params->ve[17];
	double c_5 = params->ve[18];
	double c_6 = params->ve[19];
	double lambda_1 = global_params->ve[1];

	double qlambda_1 = pow(q,lambda_1);
	double s23 = pow(s,2.0/3.0);
	double expterm = (v==1.0-a) ? 0.0 : exp(-(1-a-v)/(1-a));
	double denom = (s*(1-a) + epsilon*(1-a-v))*(s*(1-a) + epsilon*(1-a-v))*(s*(1-a) + epsilon*(1-a-v));

	//Calculate total inflow
	result[0][0] = 0.0;
	for(i=0;i<numparents;i++)
		result[0][0] += y_p[i]->ve[0];

	//Calculate the jacobian
	result[0][0] = qlambda_1 * ( lambda_1 * invtau * (result[0][0] - q + c_1*s23*s) / q - invtau );
	//result[0][0] = qlambda_1 * ( lambda_1 * invtau * (result[0][0] - q + c_1*s23*s + 100.0) / q - invtau );
	result[0][1] = 5/3 * c_1 * invtau * qlambda_1 * s23;
	result[0][2] = 0.0;
	result[0][3] = 0.0;

	result[1][0] = 0.0;
	if(denom == 0.0)
	{
		result[1][1] = c_5 * (1.0 - a);
		result[1][2] = 0.0;
		result[1][3] = 0.0;
	}
	else
	{
		result[1][1] = -5.0/3.0*c_4*s23 - c_5*(1-a)*(epsilon*(1-a-v))*(epsilon*(1-a-v))*(-s*(1-a)+epsilon*(1-a-v)) / denom;
		result[1][2] = c_5*epsilon*epsilon*s*(1-a-v)*((1-a)*s*(1+v-a)+epsilon*(1-a-v)*(1-a-v)) / denom;
		result[1][3] = 2*c_5*epsilon*epsilon*s*s*(1-a)*(1-a)*(1-a-v) / denom;
	}

	result[2][0] = 0.0;
	result[2][1] = 0.0;
	result[2][2] = (v==1.0-a) ? -2*c_6*a/epsilon : c_6/((1-a)*epsilon) * (v*(v-1)*expterm + a*v*expterm) + 2*c_6*a/epsilon;
	result[2][3] = c_6/epsilon*(1-a+v)*expterm;

	result[3][0] = 0.0;
	if(denom == 0)
	{
		result[3][1] = -c_5/epsilon * (1-a);
		result[3][2] = (a == 1.0) ? 0.0 : c_6 * v/epsilon * (1-a-v)/(1-a) * expterm;
		result[3][3] = -c_6/epsilon * (1+v-a) * expterm;
	}
	else
	{
		result[3][1] = c_5*epsilon*(1-a)*(1-a-v)*(1-a-v)*(s-a*s-epsilon*(1-a-v)) / denom;
		result[3][2] = (1-a-v)/epsilon * (c_6/(1-a)*v*expterm + c_5*epsilon*epsilon*s*((1-a)*s*(a-v-a)-epsilon*(1-a-v)*(1-a-v)) / denom);
		result[3][3] = -c_6/epsilon*(1+v-a)*expterm + 2*c_5*epsilon*s*s*(1-a)*(1-a)*(-1+a+v) / denom;
	}
}

//Type 101
void Robertson(double t,VEC* y_i,VEC** y_p,unsigned short int numparents,VEC* global_params,double* forcing_values,QVSData* qvs,VEC* params,IVEC* iparams,int state,unsigned int** upstream,unsigned int* numupstream,VEC* ans)
{
	double y1 = y_i->ve[0];
	double y2 = y_i->ve[1];
	double y3 = y_i->ve[2];

	ans->ve[0] = -.04 * y1 + 1e4 * y2 * y3;
	ans->ve[1] = .04 * y1 - 1e4 * y2 * y3 - 3e7 * y2 * y2;
	ans->ve[2] = 3e7 * y2 * y2;
}

//Type 101
void JRobertson(double t,VEC* y_i,VEC** y_p,unsigned short int numparents,VEC* global_params,double* forcing_values,VEC* params,MAT* ans)
{
	double y2 = y_i->ve[1];
	double y3 = y_i->ve[2];

	ans->me[0][0] = -.04;
	ans->me[0][1] = 1e4 * y3;
	ans->me[0][2] = 1e4 * y2;

	ans->me[1][0] = .04;
	ans->me[1][1] = -1e4 * y3 - 6e7 * y2;
	ans->me[1][2] = -1e4 * y2;

	ans->me[2][0] = 0.0;
	ans->me[2][1] = 6e7 * y2;
	ans->me[2][2] = 0.0;
}

//Type 105
//This is the model used in the Asynch summary paper
//Function for the simple river system with rainfall
//Calculates the flow using simple parameters, using only the flow q.
//Order of parameters: L_i,A_h,A_i,h_b,h_H,max_inf_rate,K_sat,S_h,eta,b_H,c_H,d_H,invtau,c_1,c_2,c_3
//The numbering is:     0   1   2   3   4       5         6    7   8   9   10  11  12    13  14  15 
//Order of global_params: v_r,lambda_1,lambda_2,Q_r,A_r,RC
//The numbering is:        0      1        2     3   4   5
//This uses the units and functions from September 18, 2011 document
//y_i[0] = q, y_i[1] = s
void river_rainfall_summary(double t,VEC* y_i,VEC** y_p,unsigned short int numparents,VEC* global_params,double* forcing_values,QVSData* qvs,VEC* params,IVEC* iparams,int state,unsigned int** upstream,unsigned int* numupstream,VEC* ans)
{
	unsigned short int i;

	double q = y_i->ve[0];
	double s = y_i->ve[1];

	double invtau = params->ve[12];
	double c_1 = params->ve[13];
	double c_2 = params->ve[14];
	double c_3 = params->ve[15];
	double lambda_1 = global_params->ve[1];

	double s53 = pow(s,5.0/3.0);

	//Flux equation (y_i[0])
	ans->ve[0] = -q + c_1 * s53;
	for(i=0;i<numparents;i++)
		ans->ve[0] += y_p[i]->ve[0];
	ans->ve[0] = invtau * pow(q,lambda_1) * ans->ve[0];

	//Ponded water equation (y_i[1])
	ans->ve[1] = c_2 * forcing_values[0] - c_3 * s53;
}



//Type 219
//This model is used for the tiling experiments. The rainfall is treated as the input from the hillslope.
//Order of parameters: A_i,L_i,A_h,invtau,c_1
//The numbering is:	0   1   2    3     4
//Order of global_params: v_0,lambda_1,lambda_2
//The numbering is:        0      1        2
void Tiling(double t,VEC* y_i,VEC** y_p,unsigned short int numparents,VEC* global_params,double* forcing_values,QVSData* qvs,VEC* params,IVEC* iparams,int state,unsigned int** upstream,unsigned int* numupstream,VEC* ans)
{
	unsigned short int i;

	double lambda_1 = global_params->ve[1];
	double A_h = params->ve[2];	//[m^2]
	double invtau = params->ve[3];	//[1/min]
	double c_1 = params->ve[4];

	double hillslope_flux = forcing_values[0] * c_1;	//[m^3/s]

	double q = y_i->ve[0];		//[m^3/s]

	//Discharge
	ans->ve[0] = -q + hillslope_flux;
	for(i=0;i<numparents;i++)
		ans->ve[0] += y_p[i]->ve[0];
	ans->ve[0] = invtau * pow(q,lambda_1) * ans->ve[0];
}


//Type 300
//Function for simple river system with data assimilation.
//Calculates the flow using simple parameters, using only the flow q.
//Order of parameters: L_i,A_h,A_i,h_b,h_H,max_inf_rate,K_sat,S_h,eta,b_H,c_H,d_H,invtau,epsilon,c_1,c_2,c_3,c_4,c_5,c_6
//The numbering is:     0   1   2   3   4       5         6    7   8   9   10  11  12    13      14  15  16  17  18  19
//Order of global_params: v_r,lambda_1,lambda_2,Q_r,A_r,RC
//The numbering is:        0      1        2     3   4   5
//This uses the units and functions from September 18, 2011 document
//y_i[0] = q, followed by N entries for the variational equation
void assim_simple_river(double t,VEC* y_i,VEC** y_p,unsigned short int numparents,VEC* global_params,double* forcing_values,QVSData* qvs,VEC* params,IVEC* iparams,int state,unsigned int** upstream,unsigned int* numupstream,VEC* ans)
{
	unsigned int i,j;
	unsigned int dim = ans->dim;
	unsigned int offset = 1;
	unsigned int loc = iparams->ve[0];
	double inflow = 0.0;

	double q = y_i->ve[0];
	double lambda_1 = global_params->ve[1];
	double invtau = params->ve[12];

	double q_to_lambda_1 = pow(q,lambda_1);

	//Flux equation (y_i[0])
	for(i=0;i<numparents;i++)
		inflow += y_p[i]->ve[0];
	ans->ve[0] = invtau * q_to_lambda_1 * (inflow - q);

	//Eqs for variational equations
	//!!!! I don't think this is really needed !!!!
	for(i=offset;i<dim;i++)	ans->ve[i] = 0.0;

	//Variables from parents
	for(i=0;i<numparents;i++)
	{
		for(j=0;j<numupstream[i];j++)
		{
			ans->ve[upstream[i][j] + offset] = invtau * q_to_lambda_1 * y_p[i]->ve[upstream[i][j] + offset];
			if(lambda_1 > 1e-12 && inflow > 1e-12)
				ans->ve[upstream[i][j] + offset] += (lambda_1 * invtau * q_to_lambda_1/q * inflow - (lambda_1 + 1) * invtau * q_to_lambda_1) * y_i->ve[upstream[i][j] + offset];
			else
				ans->ve[upstream[i][j] + offset] += -(lambda_1 + 1.0) * invtau * q_to_lambda_1 * y_i->ve[upstream[i][j] + offset];
		}
	}

	//Variable from this link
	if(lambda_1 > 1e-12 && inflow > 1e-12)
		ans->ve[loc + offset] = (lambda_1 * invtau * q_to_lambda_1/q * inflow - (lambda_1 + 1) * invtau * q_to_lambda_1) * y_i->ve[loc + offset];
	else
		ans->ve[loc + offset] = -(lambda_1 + 1.0) * invtau * q_to_lambda_1 * y_i->ve[loc + offset];
}

//Type 301
//Function for simple river system with data assimilation.
//Calculates the flow using simple parameters, using only the flow q.
//Order of parameters: L_i,A_h,A_i,h_b,h_H,max_inf_rate,K_sat,S_h,eta,b_H,c_H,d_H,invtau,epsilon,c_1,c_2,c_3,c_4,c_5,c_6
//The numbering is:     0   1   2   3   4       5         6    7   8   9   10  11  12    13      14  15  16  17  18  19
//Order of global_params: v_r,lambda_1,lambda_2,Q_r,A_r,RC
//The numbering is:        0      1        2     3   4   5
//This uses the units and functions from September 18, 2011 document
//y_i[0] = q, y_i[1] = s, followed by N entries for the variational equation
void assim_river_rainfall(double t,VEC* y_i,VEC** y_p,unsigned short int numparents,VEC* global_params,double* forcing_values,QVSData* qvs,VEC* params,IVEC* iparams,int state,unsigned int** upstream,unsigned int* numupstream,VEC* ans)
{
	unsigned int i,j;
	unsigned int dim = ans->dim;
	unsigned int offset = 2;		//!!!! This needs to be num_dense, but without variational eqs !!!!
	unsigned int parent_offset;
	unsigned int problem_dim = 2;
	unsigned int all_states = (dim-offset)/problem_dim;
	unsigned int loc = iparams->ve[0];
	double inflow = 0.0;

	double q = y_i->ve[0];
	double s = y_i->ve[1];

	double lambda_1 = global_params->ve[1];
	double invtau = params->ve[12];
	double c_1 = params->ve[14];
	double c_3 = params->ve[16];
	double c_4 = params->ve[17];

	double q_to_lambda_1 = pow(q,lambda_1);
	double q_to_lambda_1_m1 = (q > 1e-12) ? q_to_lambda_1 / q : 0.0;
	double s53 = pow(s,5.0/3.0);
	double s23 = (s53 < 1e-12) ? 0.0 : s53/s;

	//Flux equation (y_i[0])
	for(i=0;i<numparents;i++)
		inflow += y_p[i]->ve[0];
	ans->ve[0] = invtau * q_to_lambda_1 * (inflow - q + c_1 * s53);

	//Ponded water equation (y_i[1])
	ans->ve[1] = c_3 * forcing_values[0] - c_4 * s53;

	//!!!! Pull if statements out of loops !!!!
	//!!!! A lot of terms get repeated !!!!

	//Eqs for variational equations
	for(i=offset;i<dim;i++)	ans->ve[i] = 0.0;




	//New
unsigned int stophere = 1;
	//s variable from this link
	ans->ve[offset] = (-5.0/3.0)*c_4*s23*y_i->ve[offset];
//if(loc == stophere)	printf("1 putting into %u  %e\n",offset,ans->ve[offset]);

	//q variables from this link
	if(lambda_1 > 1e-12 && (inflow + c_1*s53) > 1e-12)
		ans->ve[offset + 1] = (lambda_1 * invtau * q_to_lambda_1_m1 * (inflow + c_1*s53) - (lambda_1 + 1) * invtau * q_to_lambda_1) * y_i->ve[offset + 1];
	else
		ans->ve[offset + 1] = -(lambda_1 + 1.0) * invtau * q_to_lambda_1 * y_i->ve[offset + 1];
//if(loc == stophere)	printf("2 putting into %u  %e\n",offset+1,ans->ve[offset + 1]);

	if(lambda_1 > 1e-12 && (inflow + c_1*s53) > 1e-12)
		ans->ve[offset + 2] = (lambda_1 * invtau * q_to_lambda_1_m1 * (inflow + c_1*s53) - (lambda_1 + 1) * invtau * q_to_lambda_1) * y_i->ve[offset + 2] + (5.0/3.0)*invtau*c_1*q_to_lambda_1*s23 * y_i->ve[offset];
	else
		//ans->ve[offset + 2] = -(lambda_1 + 1.0) * invtau * q_to_lambda_1 * y_i->ve[offset + 2];
		ans->ve[offset + 2] = -(lambda_1 + 1.0) * invtau * q_to_lambda_1 * y_i->ve[offset + 2] + invtau*(5.0/3.0)*c_1*s23*y_i->ve[offset];
//if(loc == stophere && t > 60.0)	printf("3 putting into %u  %e  %e %e\n",offset+2,ans->ve[offset + 2],q_to_lambda_1,y_i->ve[offset + 2]);

	//Adjust offset
	offset += 3;

	//Variables from parents
	for(i=0;i<numparents;i++)
	{
		parent_offset = 1 + problem_dim;

		for(j=0;j<numupstream[i];j++)
		{
//if(loc == stophere)	printf("4 putting into %u using %u with parent %u\n",offset,parent_offset,i);
			ans->ve[offset] = invtau * q_to_lambda_1 * y_p[i]->ve[parent_offset];
			if(lambda_1 > 1e-12 && (inflow + c_1*s53) > 1e-12)
				ans->ve[offset] += (lambda_1 * invtau * q_to_lambda_1_m1 * (inflow + c_1*s53) - (lambda_1 + 1) * invtau * q_to_lambda_1) * y_i->ve[offset];
			else
				ans->ve[offset] += -(lambda_1 + 1.0) * invtau * q_to_lambda_1 * y_i->ve[offset];

//if(loc == stophere)	printf("5 putting into %u using %u\n",offset+1,parent_offset + 1);
			ans->ve[offset + 1] = invtau * q_to_lambda_1 * y_p[i]->ve[parent_offset + 1];
			if(lambda_1 > 1e-12 && (inflow + c_1*s53) > 1e-12)
				ans->ve[offset + 1] += (lambda_1 * invtau * q_to_lambda_1_m1 * (inflow + c_1*s53) - (lambda_1 + 1) * invtau * q_to_lambda_1) * y_i->ve[offset + 1];
			else
				ans->ve[offset + 1] += -(lambda_1 + 1.0) * invtau * q_to_lambda_1 * y_i->ve[offset + 1];

			offset += 2;
			parent_offset += 2;
		}
	}
/*
//if(loc == stophere)
if(loc == stophere && t > 60.0)
{
printf("ans is\n");
Print_Vector(ans);
printf("Waiting...\n");
getchar();
}
*/



/*
	//Old
unsigned int stophere = -1;
	//q variables
	//Variables from parents
	for(i=0;i<numparents;i++)
	{
		for(j=0;j<numupstream[i];j++)
		{
			ans->ve[upstream[i][j]*problem_dim + offset] = invtau * q_to_lambda_1 * y_p[i]->ve[upstream[i][j]*problem_dim + offset];
			if(lambda_1 > 1e-12 && (inflow + c_1*s53) > 1e-12)
				ans->ve[upstream[i][j]*problem_dim + offset] += (lambda_1 * invtau * q_to_lambda_1_m1 * (inflow + c_1*s53) - (lambda_1 + 1) * invtau * q_to_lambda_1) * y_i->ve[upstream[i][j]*problem_dim + offset];
			else
				ans->ve[upstream[i][j]*problem_dim + offset] += -(lambda_1 + 1.0) * invtau * q_to_lambda_1 * y_i->ve[upstream[i][j]*problem_dim + offset];

			ans->ve[upstream[i][j]*problem_dim + 1 + offset] = invtau * q_to_lambda_1 * y_p[i]->ve[upstream[i][j]*problem_dim + 1 + offset];
			if(lambda_1 > 1e-12 && (inflow + c_1*s53) > 1e-12)
				ans->ve[upstream[i][j]*problem_dim + 1 + offset] += (lambda_1 * invtau * q_to_lambda_1_m1 * (inflow + c_1*s53) - (lambda_1 + 1) * invtau * q_to_lambda_1) * y_i->ve[upstream[i][j]*problem_dim + 1 + offset];
			else
				ans->ve[upstream[i][j]*problem_dim + 1 + offset] += -(lambda_1 + 1.0) * invtau * q_to_lambda_1 * y_i->ve[upstream[i][j]*problem_dim + 1 + offset];
		}
	}

	//Variable from this link
	if(lambda_1 > 1e-12 && (inflow + c_1*s53) > 1e-12)
		ans->ve[loc*problem_dim + offset] = (lambda_1 * invtau * q_to_lambda_1_m1 * (inflow + c_1*s53) - (lambda_1 + 1) * invtau * q_to_lambda_1) * y_i->ve[loc*problem_dim + offset];
	else
		ans->ve[loc*problem_dim + offset] = -(lambda_1 + 1.0) * invtau * q_to_lambda_1 * y_i->ve[loc*problem_dim + offset];

	if(lambda_1 > 1e-12 && (inflow + c_1*s53) > 1e-12)
		ans->ve[loc*problem_dim + 1 + offset] = (lambda_1 * invtau * q_to_lambda_1_m1 * (inflow + c_1*s53) - (lambda_1 + 1) * invtau * q_to_lambda_1) * y_i->ve[loc*problem_dim + 1 + offset] + (5.0/3.0)*invtau*c_1*q_to_lambda_1*s23 * y_i->ve[loc*problem_dim + 1 + offset + all_states];
	else
		ans->ve[loc*problem_dim + 1 + offset] = -(lambda_1 + 1.0) * invtau * q_to_lambda_1 * y_i->ve[loc*problem_dim + 1 + offset];

	//s variables
	ans->ve[loc*problem_dim + 1 + offset + all_states] = (-5.0/3.0)*c_4*s23*y_i->ve[loc*problem_dim + 1 + offset + all_states];

if(loc == stophere)
{
printf("ans is\n");
Print_Vector(ans);
printf("Waiting...\n");
getchar();
}
*/

}





// Consistency routines **************************************************************************

void CheckConsistency_Nonzero_1States(VEC* y,VEC* params,VEC* global_params)
{
	if(y->ve[0] < 1e-14)	y->ve[0] = 1e-14;
}

void CheckConsistency_Nonzero_2States(VEC* y,VEC* params,VEC* global_params)
{
	if(y->ve[0] < 1e-14)	y->ve[0] = 1e-14;
	if(y->ve[1] < 0.0)	y->ve[1] = 0.0;
}

void CheckConsistency_Nonzero_3States(VEC* y,VEC* params,VEC* global_params)
{
	if(y->ve[0] < 1e-14)	y->ve[0] = 1e-14;
	if(y->ve[1] < 0.0)	y->ve[1] = 0.0;
	if(y->ve[2] < 0.0)	y->ve[2] = 0.0;
}

void CheckConsistency_Nonzero_4States(VEC* y,VEC* params,VEC* global_params)
{
	if(y->ve[0] < 1e-14)	y->ve[0] = 1e-14;
	if(y->ve[1] < 0.0)	y->ve[1] = 0.0;
	if(y->ve[2] < 0.0)	y->ve[2] = 0.0;
	if(y->ve[3] < 0.0)	y->ve[3] = 0.0;
}

void CheckConsistency_Model5(VEC* y,VEC* params,VEC* global_params)
{
	if(y->ve[0] < 1e-14)	y->ve[0] = 1e-14;
	if(y->ve[1] < 0.0)	y->ve[1] = 0.0;
	if(y->ve[2] < 0.0)	y->ve[2] = 0.0;
	if(y->ve[3] < 0.0)	y->ve[3] = 0.0;
	if(y->ve[3] > 1.0 - y->ve[2])	y->ve[3] = 1.0 - y->ve[2];
}

void CheckConsistency_Model30(VEC* y,VEC* params,VEC* global_params)
{
	if(y->ve[0] < 1e-14)	y->ve[0] = 1e-14;
	if(y->ve[1] < 0.0)	y->ve[1] = 0.0;
	if(y->ve[2] < 0.0)	y->ve[2] = 0.0;
	if(y->ve[2] > params->ve[4])	y->ve[2] = params->ve[4];
	if(y->ve[3] < 0.0)	y->ve[3] = 0.0;
	else if(y->ve[3] > 1.0)	y->ve[3] = 1.0;
}

void CheckConsistency_Nonzero_AllStates_q(VEC* y,VEC* params,VEC* global_params)
{
	unsigned int i;

	if(y->ve[0] < 1e-14)	y->ve[0] = 1e-14;
	for(i=1;i<y->dim;i++)
		if(y->ve[i] < 0.0)	y->ve[i] = 0.0;
}




