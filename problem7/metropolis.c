#include "headers.h"
#include<math.h>
#include<stdlib.h>
int fdep (int *ptl, int x, int L)
{
    int epo,epf;
    int tmp,dep;
    
    /*energy before flipping*/
    epo=fep3(ptl,x,L);
    /* flip spin */
    tmp=*(ptl+x);
    *(ptl+x)=-tmp;
    /*energy after flipping*/
    epf=fep3(ptl,x,L);
    /*energy change in this flipping*/
    dep=epf-epo;
    return dep;
}

int flip(int *ptl, int * array, int x, int L)
{
    int dep,ep,mag;
    int spin;
    
    ep=*array;mag=*(array+1);
    
    /* flip spin; calculate energy change*/
    dep=fdep(ptl,x,L);
    
    spin=*(ptl+x); /* value of the current spin */
    
    ep+=dep;mag+=2.0*spin;
    
    *array=ep;*(array+1)=mag;
    
    return dep;
    
}

void metropolis(int *ptl, double t, int * array, int L)
{
    double r;
    int dep,tmp,ep,mag;
    int x,spin; /* position and value of the current spin */
    int N=L*L;
    
    /*choose a spin to flip*/
    x=rand() % N;
    /* flip spin;return energy change*/
    dep=flip(ptl,array,x,L);
    
    if(t<0){/*do nothing for infinite temperature; just flip randomly*/}
    else
    {
        r=(rand()%1000000)/1000000.0;
        
        if(r>exp(-dep/t))
        {/* flip back */
            dep=flip(ptl,array,x,L);
            /*printf("do not accept\n");*/
        }
    }
    
}


