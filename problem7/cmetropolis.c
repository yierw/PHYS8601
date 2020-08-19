#include "headers.h"
#include<math.h>
#include<stdlib.h>

double dep(double * p_th, double * p_ph, double th, double ph, int x, int L)
{
    /*calculate energy change if flip */
    double dep;
    double th_old,ph_old;
    double e_old,e_flip;
    int spin;
    /*3-term energy before flip*/
    e_old=cfep3(p_th,p_ph,x, L);
    /* store current spin info*/
    th_old=*(p_th+x);
    ph_old=*(p_ph+x);
    /*flip spin*/
    *(p_th+x)=th;
    *(p_ph+x)=ph;
    /*3-term energy after flip*/
    e_flip=cfep3(p_th,p_ph,x, L);
    /*flip spin back*/
    *(p_th+x)=th_old;
    *(p_ph+x)=ph_old;
    /*calculate energy change*/
    dep=e_flip-e_old;
    return dep;
}

void cmetropolis(double * p_th, double * p_ph, double t, double * array, int L, int flag)
{
    double pi=3.1415926;
    double e0,mx,my,mz;
    double th_old,ph_old,th_new,ph_new;
    double r,dE;
    int x;/* position of current spin */
    int N=L*L;
    /* choose a spin to flip */
    x=rand() % N;
    /* store current info */
    th_old=*(p_th+x);
    ph_old=*(p_ph+x);
    e0=*array;
    mx=*(array+1);
    my=*(array+2);
    mz=*(array+3);
    /* choose the new angular components*/
    if(flag==2)
    {
        th_new=pi*(rand()%1000000)/1000000.0;
        ph_new=2.0*pi*(rand()%1000000)/1000000.0;
        
    }
    else{
        th_new=pi-th_old;
    }
    
    /* calculate the energy change */
    dE=dep(p_th, p_ph, th_new, ph_new, x, L);
    
    r=(rand()%1000000)/1000000.0;
    if(r<=exp(-dE/t))
    {
        /* flip spin */
        *(p_th+x)=th_new;
        *(p_ph+x)=ph_new;
        /* update array */
        *array=e0+dE;/*update energy*/
        *(array+1)=mx-sin(th_old)*cos(ph_old)+sin(th_new)*cos(ph_new);
        *(array+2)=my-sin(th_old)*sin(ph_old)+sin(th_new)*sin(ph_new);
        *(array+3)=mz-cos(th_old)+cos(th_new);
    }
    
}


