#ifndef H_INCLUDED
#define H_INCLUDED
#include <stdio.h>

/* !!! KEY FUNCTIONS DO NOT CHANGE !!! */
void initial(int * ptl, float p, int N); /* INITIALIZE LATTICE */
void wrtl (FILE * fp, int *ptl, int printl, int L); /* OUTPUT LATTICE OCCUPATION */
int fmag(int *ptl, int N);
int fener(int *ptl, int L);
int fepx (int *ptl, int x, int L);/* calculate the nearest-neighbor interaction energy of site (row,col) periodic boundary conditions applied */
int fep3(int *ptl, int x, int L);
int fdep (int *ptl, int x, int L);/*return energy change per site */
void  insert_sort(int * a, int * b, int length);
void  insert_sort3(int * a, int * b, double *c, int length);
double min(double * x,  int length);

/* sampling related functions */
int flip(int *ptl, int * array, int x, int L);
void metropolis(int *ptl, double t, int * array, int L);
int hist(int E, int * x, int * num, double * g, int count);/* for Wang-Landau Sampling */
double diff(int * x,  int length);


/* other functions  */
double fcv1( int * E, double * P, double k, int N, int length);
void reweight(FILE * fp3,FILE * fp4, int * pointer,double t, int len, double * k,int klen, int N);/* histogram reweighting */
void correlation(FILE *fp1, int * earr, int * marr, int clen, int len);
void thermal1(FILE *fp5, int * E, int * M,double t, int N, int length);
void thermal2(FILE *fp5, double * P, int * E, double t, int N, int length); /* for Wang-Landau Sampling */

#endif
