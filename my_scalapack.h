#ifndef MY_SCALAPACK_H
#define MY_SCALAPACK_H

typedef char* F_CHAR_T;

void pdgemm_(F_CHAR_T,F_CHAR_T,int*,int*,int*,double*,double*,int*,int*,int*,double*,int*,int*,int*,double*,double*,int*,int*,int*);
void pdtran_(int*,int*,double*,double*,int*,int*,int*,double*,double*,int*,int*,int*);
void pdgemv_(F_CHAR_T,int*,int*,double*,double*,int*,int*,int*,double*,int*,int*,int*,int*,double*,double*,int*,int*,int*,int*);
void pdsyevr_(F_CHAR_T,F_CHAR_T,F_CHAR_T,int*,double*,int*,int*,int*,double*,double*,int*,int*,int*,int*,double*,double*,int*,int*,int*,double*,int*,int*,int*,int*);
int numroc_(int*,int*,int*,int*,int*);

//CBlacs
void Cblacs_get(int ConTxt, int what, int *val);
void Cblacs_gridinit(int *ConTxt, char *order, int nprow, int npcol);

#endif

