#include "mathmethods.h"

//Allocates space for a vector with size entries.
VEC* v_get(int size)
{
	VEC* v = (VEC*) malloc(sizeof(VEC));
	v->ve = (double*) calloc(size,sizeof(double));
	v->dim = size;
	return v;
}

//Allocates space for an ivector with size entries.
IVEC* iv_get(int size)
{
	IVEC* v = (IVEC*) malloc(sizeof(IVEC));
	v->ve = (int*) calloc(size,sizeof(int));
	v->dim = size;
	return v;
}

//Allocates space for an m x n matrix.
MAT* m_get(int m,int n)
{
	int i;
	MAT* A = (MAT*) malloc(sizeof(MAT));
	A->array = (double*) calloc(m*n,sizeof(double));
	A->me = (double**) malloc(m*sizeof(double*));
	for(i=0;i<m;i++)	A->me[i] = &(A->array[i*n]);
	A->m = m;
	A->n = n;
	return A;
}

//Deallocates the vector v.
void v_free(VEC* v)
{
	if(v == NULL)	return;
	free(v->ve);
	free(v);
}

//Deallocates the ivector v.
void iv_free(IVEC* v)
{
	if(v == NULL)	return;
	free(v->ve);
	free(v);
}

//Deallocates the matrix A.
void m_free(MAT* A)
{
	if(A == NULL)	return;
	free(A->array);
	free(A->me);
	free(A);
}

//Calculates w = u + v. start is the index of the first entry of u.
void v_add(VEC* u,VEC* v,VEC* w,unsigned int start)
{
	unsigned int i;
	double* a = &(u->ve[start]);
	double* b = v->ve;
	double* c = w->ve;
	unsigned int size = w->dim;

	for(i=0;i<size;i++)
		c[i] = a[i] + b[i];
}

//Calculates C = A + B
void m_add(MAT* A,MAT* B,MAT* C)
{
	unsigned int i,j;
	unsigned int m = C->m;
	unsigned int n = C->n;

	for(i=0;i<m;i++)
		for(j=0;j<n;j++)	C->me[i][j] = A->me[i][j] + B->me[i][j];
}

//Calculates v = A * v. Assumes dimensions work
void mv_mlt(MAT* A,VEC* u,VEC* v)
{
	unsigned int i,j;
	unsigned int m = A->m;
	unsigned int n = A->n;

	for(i=0;i<m;i++)
	{
		v->ve[i] = A->me[i][0] * u->ve[0];
		for(j=1;j<n;j++)	v->ve[i] += A->me[i][j] * u->ve[j];
	}
}

//Calculates C = A * B. Assumes dimensions work.
void mm_mlt(MAT* A,MAT* B,MAT* C)
{
	unsigned int i,j,k;
	unsigned int m = A->m;
	unsigned int n = B->n;
	unsigned int inner = A->n;

	for(i=0;i<m;i++)
	{
		for(j=0;j<n;j++)	C->me[i][j] = 0.0;
		for(k=0;k<inner;k++)
		{
			for(j=0;j<n;j++)
				C->me[i][j] += A->me[i][k] * B->me[k][j];
		}
	}

/*
	for(i=0;i<m;i++)
	{
		for(j=0;j<n;j++)
		{
			C->me[i][j] = 0.0;
			for(k=0;k<inner;k++)
				C->me[i][j] += A->me[i][k] * B->me[k][j];
		}
	}
*/
}

//Calculates C = A^T * B. Assumes dimensions work.
void mTm_mlt(MAT* A,MAT* B,MAT* C)
{
	unsigned int i,j,k;
	unsigned int m = A->n;
	unsigned int n = B->n;
	unsigned int inner = A->m;

	for(i=0;i<m;i++)
	{
		for(j=0;j<n;j++)	C->me[i][j] = 0.0;
		for(k=0;k<inner;k++)
		{
			for(j=0;j<n;j++)
				C->me[i][j] += A->me[k][i] * B->me[k][j];
		}
	}

/*
	for(i=0;i<m;i++)
	{
		for(j=0;j<n;j++)
		{
			C->me[i][j] = 0.0;
			for(k=0;k<inner;k++)
				C->me[i][j] += A->me[k][i] * B->me[k][j];
		}
	}
*/
}

//Calculates C = A * B^T. Assumes dimensions work.
void mmT_mlt(MAT* A,MAT* B,MAT* C)
{
	unsigned int i,j,k;
	unsigned int m = A->m;
	unsigned int n = B->m;
	unsigned int inner = A->n;

	for(i=0;i<m;i++)
	{
		for(j=0;j<n;j++)
		{
			C->me[i][j] = 0.0;
			for(k=0;k<inner;k++)
				C->me[i][j] += A->me[i][k] * B->me[j][k];
		}
	}

/*
	for(i=0;i<m;i++)
	{
		for(j=0;j<n;j++)
		{
			C->me[i][j] = 0.0;
			for(k=0;k<inner;k++)
				C->me[i][j] += A->me[i][k] * B->me[j][k];
		}
	}
*/
}

//Calculates w = u - v. start is the index of the first entries.
void v_sub(VEC* u,VEC* v,VEC* w,unsigned int start)
{
	unsigned int i;
	double* a = u->ve;
	double* b = v->ve;
	double* c = w->ve;
	unsigned int size = w->dim;

	for(i=start;i<size;i++)
		c[i] = a[i] - b[i];
}

//Calculates B = alpha * A
void sm_mlt(double alpha,MAT* A,MAT* B,unsigned int startrow,unsigned int startcol)
{
	unsigned int i,j;
	unsigned int m = A->m;
	unsigned int n = A->n;
	for(i=0;i<m;i++)
		for(j=0;j<n;j++)	B->me[i+startrow][j+startcol] = alpha * A->me[i][j];
}

//Calculates B = I + alpha * A
void dipaa(double alpha,MAT* A,MAT* B,unsigned int startrow,unsigned int startcol)
{
	unsigned int i,j;
	unsigned int m = A->m;
	unsigned int n = A->n;
	for(i=0;i<m;i++)
	{
		for(j=0;j<n;j++)	B->me[i+startrow][j+startcol] = alpha * A->me[i][j];
		B->me[i+startrow][i+startcol] += 1.0;
	}
}

//Calculates A * B = C, assuming A, B, C are diagonal.
void diag_mm_mlt(MAT* A,MAT* B,MAT* C)
{
	unsigned int i;
	unsigned int size = (A->m < A->n) ? A->m : A->n;
	size = (size < B->n) ? size : B->n;
	for(i=0;i<size;i++)
		C->me[i][i] = A->me[i][i] * B->me[i][i];
}

//Calculates A * B^T = C, assuming A, B, C are diagonal.
void diag_mm_mltT(MAT* A,MAT* B,MAT* C)
{
	unsigned int i;
	unsigned int size = (A->m < A->n) ? A->m : A->n;
	size = (size < B->m) ? size : B->m;
	for(i=0;i<size;i++)
		C->me[i][i] = A->me[i][i] * B->me[i][i];
}

//Calculates A + B = C, assuming A, B, C are diagonal.
void diag_m_add(MAT* A,MAT* B,MAT* C)
{
	unsigned int i;
	unsigned int size = (A->m < A->n) ? A->m : A->n;
	for(i=0;i<size;i++)
		C->me[i][i] = A->me[i][i] + B->me[i][i];
}

//Calculates B = A^-1, assuming A is diagonal.
void diag_invert(MAT* A,MAT* B)
{
	unsigned int i;
	unsigned m = A->m;
	for(i=0;i<m;i++)
		B->me[i][i] = 1.0/A->me[i][i];
}

//Calculates v = A * u, assuming A is diagonal.
void diag_mv_mlt(MAT* A,VEC* u,VEC* v)
{
	unsigned int i;
	unsigned int size = (A->m < A->n) ? A->m : A->n;
	for(i=0;i<size;i++)
		v->ve[i] = A->me[i][i] * u->ve[i];
}


//Calculates C = A^T * B, assuming B is diagonal.
void mTdiag_mlt(MAT* A,VEC* B,MAT* C)
{
	unsigned int i,j;
	unsigned int m = A->n;
	unsigned int n = B->dim;

	for(j=0;j<n;j++)
		for(i=0;i<m;i++)
			C->me[i][j] = A->me[j][i] * B->ve[j];
}

//Calculates C = A * B, assuming B is diagonal
void mdiag_mlt(MAT* A,VEC* B,MAT* C)
{
	unsigned int i,j;
	unsigned int m = A->m;
	unsigned int n = B->dim;

	for(i=0;i<m;i++)
	{
		for(j=0;j<n;j++)
			C->me[i][j] = A->me[i][j] * B->ve[j];
	}
}

//Calculates C = A * B, assuming A is diagonal
void diagm_mlt(VEC* A,MAT* B,MAT* C)
{
	unsigned int i,j;
	unsigned int m = A->dim;
	unsigned int n = B->n;

	for(i=0;i<m;i++)
	{
		for(j=0;j<n;j++)
			C->me[i][j] = A->ve[i] * B->me[i][j];
	}
}


//Calculates A = V*D*V^T, where D is diagonal and square
void VDVT_mlt(MAT* V,VEC* D,MAT* A)
{
	unsigned int i,j,k;
	unsigned int m = A->m;
	unsigned int n = A->n;
	unsigned int inner = D->dim;

	for(i=0;i<m;i++)
	{
		for(j=0;j<n;j++)
		{
			A->me[i][j] = 0.0;
			for(k=0;k<inner;k++)	A->me[i][j] += V->me[i][k] * V->me[j][k] * D->ve[k];
		}
	}
}

//Calculates A = V*D^-1*V^T, where D is diagonal and square
void VDinvVT_mlt(MAT* V,VEC* D,MAT* A,VEC* temp)
{
	unsigned int i,j,k;
	unsigned int m = A->m;
	unsigned int n = A->n;
	unsigned int inner = D->dim;

	for(i=0;i<m;i++)
		temp->ve[i] = 1.0/D->ve[i];

	for(i=0;i<m;i++)
	{
		for(j=0;j<n;j++)
		{
			A->me[i][j] = 0.0;
			for(k=0;k<inner;k++)	A->me[i][j] += V->me[i][k] * V->me[j][k] * temp->ve[k];
		}
	}
}

//Calculates A = V*sqrt(D)^-1*V^T, where D is diagonal and square
void VsqrtDinvVT_mlt(MAT* V,VEC* D,MAT* A,VEC* temp)
{
	unsigned int i,j,k;
	unsigned int m = A->m;
	unsigned int n = A->n;
	unsigned int inner = D->dim;

	for(i=0;i<m;i++)
		temp->ve[i] = 1.0/sqrt(D->ve[i]);

	for(i=0;i<m;i++)
	{
		for(j=0;j<n;j++)
		{
			A->me[i][j] = 0.0;
			for(k=0;k<inner;k++)	A->me[i][j] += V->me[i][k] * V->me[j][k] * temp->ve[k];
		}
	}
}

//Calculates A = V*sqrt(D)*V^T, where D is diagonal and square
void VsqrtDVT_mlt(MAT* V,VEC* D,MAT* A,VEC* temp)
{
	unsigned int i,j,k;
	unsigned int m = A->m;
	unsigned int n = A->n;
	unsigned int inner = D->dim;

	for(i=0;i<m;i++)
		temp->ve[i] = sqrt(D->ve[i]);

	for(i=0;i<m;i++)
	{
		for(j=0;j<n;j++)
		{
			A->me[i][j] = 0.0;
			for(k=0;k<inner;k++)	A->me[i][j] += V->me[i][k] * V->me[j][k] * temp->ve[k];
		}
	}
}

//Computes v^T*A*v
double vTAv(MAT* A,VEC* v)
{
	unsigned int i,j;
	double result = 0.0;
	unsigned size = A->m;

	for(i=0;i<size;i++)
		for(j=0;j<size;j++)
			result += v->ve[i] * A->me[i][j] * v->ve[j];

	return result;
}


//Copies the contents of v into w.
//Note: This seems to be a bit faster than memcpy.
void v_copy(VEC* v,VEC* w)
{
	unsigned int i;
	double* a = v->ve;
	double* b = w->ve;
	unsigned int size = w->dim;

	for(i=0;i<size;i++)
		b[i] = a[i];
}

//Copies the contents of A into B.
void m_copy(MAT* A,MAT* B)
{
	unsigned int i;
	double* a = A->array;
	double* b = B->array;
	unsigned int size = A->m * A->n;

	for(i=0;i<size;i++)
		b[i] = a[i];
}

//Calculate the weights for the Lagrange interpolation polynomial
//Assumes that 0 and each entry of c are node points
//This does not store the weight for 0
VEC* lagrange_weights(short unsigned int s,VEC* c)
{
	unsigned int i,j;
	VEC* w = v_get(s);

	for(i=0;i<s;i++)
	{
		//w->ve[i] = 1.0;
		w->ve[i] = 1.0/c->ve[i];	//For the 0 node
		for(j=0;j<i;j++)	w->ve[i] *= 1.0/(c->ve[i] - c->ve[j]);
		for(j=i+1;j<s;j++)	w->ve[i] *= 1.0/(c->ve[i] - c->ve[j]);
	}

	return w;
}
/*
//Evaluates the Lagrange polynomial at theta using the barycentric formula
double lagrange_bary(double theta,short unsigned int s,VEC* c,VEC** Q,VEC* w)
{
	unsigned int i;
	double sum = 0.0;
	double prod = 1.0;

	for(i=0;i<s;i++)
	{
		if( fabs(theta - c->ve[i]) < 1e-14)	return Q[i]->ve[0];
		prod *= (theta - c->ve[i]);
		sum += w->ve[i] / (theta - c->ve[i]) * Q[i]->ve[0];
	}

	return prod * sum;
}
*/

//Evaluates the Lagrange polynomial at theta using the barycentric formula.
//Assumes 0 and each component of c are the node points. At theta = 0, the poly returns 0.
//This assumes that theta will be larger than 1.
void lagrange_bary(double theta,VEC* c,VEC** Z,VEC* w,VEC* sum)
{
	unsigned short int i;
	unsigned short int s = w->dim;
	unsigned short int dim = sum->dim;
	for(i=0;i<dim;i++)	sum->ve[i] = 0.0;
	double prod = theta;	//This is for the value at theta = 0

	for(i=0;i<s;i++)
	{
		prod *= (theta - c->ve[i]);
		daxpy(w->ve[i] / (theta - c->ve[i]),Z[i],sum,0);
	}

	for(i=0;i<dim;i++)	sum->ve[i] *= prod;
}

//Evaluates the Lagrange polynomial at theta
double lagrange(double theta,short unsigned int s,VEC* c,VEC** Q)
{
	unsigned int i,j;
	double sum = 0.0;
	double product;

	for(i=0;i<s;i++)
	{
		product = 1.0;
		for(j=0;j<i;j++)	product *= (theta - c->ve[j]) / (c->ve[i] - c->ve[j]);
		for(j=i+1;j<s;j++)	product *= (theta - c->ve[j]) / (c->ve[i] - c->ve[j]);
		sum += Q[i]->ve[0] * product;
	}

	return sum;
}

//Computes the infinity norm of v. v_i is divided first by w_i.
double norm_inf(VEC* v,VEC* w,unsigned int start)
{
	unsigned int i;
	double max = fabs(v->ve[start]/w->ve[start]);
	double val;
	for(i=start+1;i<v->dim;i++)
	{
		val = fabs(v->ve[i]/w->ve[i]);
		max = (val > max) ? val : max;
	}

	return max;
}

//Computes the infinity norm of the vector v.
double vector_norminf(VEC* v,unsigned int start)
{
	unsigned int i;
	double norm = fabs(v->ve[start]);
	double val;
	for(i=start+1;i<v->dim;i++)
	{
		val = fabs(v->ve[i]);
		norm = (norm < val) ? val : norm;
	}
	return norm;
}

//Prints the vector v to stdout.
void Print_Vector(VEC* v)
{
	unsigned int i;
	for(i=0;i<v->dim;i++)
		printf("%.16f ",v->ve[i]);
	printf("\n");
}

//Prints the matrix A to stdout.
void Print_Matrix(MAT* A)
{
	unsigned int i,j;
	for(i=0;i<A->m;i++)
	{
		for(j=0;j<A->n;j++)	printf("%.16f ",A->me[i][j]);
		printf("; \n");
	}
}

//Prints the vector v in C format to stdout.
void Print_VectorC(VEC* v)
{
	unsigned int i;
	printf("{ %.16f",v->ve[0]);
	for(i=1;i<v->dim;i++)
		printf(", %.16f",v->ve[i]);
	printf(" }\n");
}

//Prints the matrix A in C format to stdout.
void Print_MatrixC(MAT* A)
{
	unsigned int i,j;
	printf("{ ");
	for(i=0;i<A->m;i++)
	{
		printf("{ %.16f",A->me[i][0]);
		for(j=1;j<A->n;j++)	printf(", %.16f",A->me[i][j]);
		printf("},\n");
	}
	printf("}\n");
}

