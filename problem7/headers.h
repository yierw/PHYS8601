#ifndef H_INCLUDED
#define H_INCLUDED
#include <stdio.h>

/* for quantum spin */
void initial(int * ptl, float p, int N); /* INITIALIZE LATTICE */
void wrtl (FILE * fp, int *ptl, int printl, int L); /* OUTPUT LATTICE OCCUPATION */
int fmag(int *ptl, int N);
int fener(int *ptl, int L);
int fepx (int *ptl, int x, int L);
int fep3(int *ptl, int x, int L);
int fdep (int *ptl, int x, int L);/*return energy change per site */
int flip(int *ptl, int * array, int x, int L);
void metropolis(int *ptl, double t, int * array, int L);
void thermal1(FILE *fp5, double * E, double * M,double t, int N, int length);

/* for classical spins*/
void cmetropolis(double * p_th, double * p_ph, double t, double * array, int L, int flag);
double dep(double * p_th, double * p_ph, double th, double ph, int x, int L);
double ss(double th1, double ph1, double th2, double ph2);
void iarr(double * p_th, double * p_ph, double * array,int N);
double cfepx (double * p_th, double * p_ph, int x, int L);
double cfep3(double * p_th, double * p_ph, int x, int L);

#endif
