#include "headers.h"
#include<math.h>
#include<stdlib.h>

void thermal1(FILE *fp5, double * E, double * M,double t, int N, int length)
{
    double cv,e2bar,ebar;
    double ch,m2bar,mbar;
    double k;
    int i;
    
    if(t<0){k=0.0;/* infinite temperature*/}
    else if (t==0.0){printf("cannot calculate zero temperature\n");exit(1);}
    else{k=1.0/t;}
    
    ebar=0.0;e2bar=0.0;
    
    for(i=0;i<length;i++)
    {/* unscaled by lattice size */
        ebar+= E[i] ;
        e2bar+= E[i] *  E[i];
        mbar+= M[i];
        m2bar+= M[i] *  M[i];
    }
    ebar/= (double) length;
    e2bar/= (double) length;
    mbar/= (double) length;
    m2bar/= (double) length;
    cv= (e2bar-ebar*ebar)*k*k;
    ch= (m2bar-mbar*mbar)*k;
    
    /* scaled by lattice size */
    ebar/=(double) N;
    mbar/=(double) N;
    cv/=(double) N;
    ch/=(double) N;
    
    fprintf(fp5," %10.6f %10.6f %10.6f %10.6f %10.6f\n", t, ebar, cv,mbar,ch);
}




