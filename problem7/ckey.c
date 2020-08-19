#include "headers.h"
#include<stdlib.h>
#include<math.h>

double ss(double th1, double ph1, double th2, double ph2)
{
    /*Dot product for two spin vectors of unit length*/
    /* th: 0~pi */
    /* ph: 0~2pi*/
    /* spin up ~ th = 0; spin down ~ th = pi*/
    double c1,c2,d1,d2;
    double ss;
    
    c1=sin(th1)*sin(th2);
    c2=cos(th1)*cos(th2);
    d1=sin(ph1)*sin(ph2);
    d2=cos(ph1)*cos(ph2);
    ss=c1*(d1+d2)+c2;    
    return ss;
}

void iarr(double * p_th, double * p_ph, double * array,int L)
{
    int i;
    int N=L*L;
    double ener,mx,my,mz;
    double th,ph;

    ener=0.0;
    mx=0.0;my=0.0;mz=0.0;
    for(i=0;i<N;i++)
    {
        ener+=cfepx(p_th,p_ph,i,L);
        th=*(p_th+i);ph=*(p_ph+i);
        mx+=sin(th)*cos(ph);
        my+=sin(th)*sin(ph);
        mz+=cos(ph);
    }
    *(array)=ener;
    *(array+1)=mx;
    *(array+2)=my;
    *(array+3)=mz;

}

double cfepx (double * p_th, double * p_ph, int x, int L)
{
    /* return the value of   - s_x * s_a - s_x * s_r    */
    int m,a,r;/*position in the lattice*/
    int N=L*L;
    int col=x%L;
    int row=(x-col)/L;
    double ma;/* s_x dot product s_above */
    double mr;/* s_x dot product s_right */
    double epx;
    double thm,phm;/*for middle point*/
    double tha,pha;/*for above point*/
    double thr,phr;/*for right point*/
    
    /* apply periodic boundary */
    /* middle site */
    m=x;
    thm=*(p_th+m);
    phm=*(p_ph+m);
    /* above site */
    if(row==0){a=x-L+N;}
    else{a=x-L;}
    tha=*(p_th+a);
    pha=*(p_ph+a);
    /* right site */
    if(col==(L-1)){r=x+1-L;}
    else{r=x+1;}
    thr=*(p_th+r);
    phr=*(p_ph+r);
    
    ma=ss(thm,phm,tha,pha);
    mr=ss(thm,phm,thr,phr);
    epx=-ma-mr;
        
    return epx;
}

double cfep3(double * p_th, double * p_ph, int x, int L)
{
    /*return the part that is going to be influenced by the change of the current spin*/
    int m,b,l;/*position in the lattice*/
    int N=L*L;
    int col=x%L;
    int row=(x-col)/L;
    double ep3;
    
    /* middle site */
    m=x;
    /* below site */
    if(row==(L-1)){b=x+L-N;}
    else{b=x+L;}
    /*left site */
    if(col==0){l=x-1+L;}
    else{l=x-1;}
    
    ep3=(cfepx(p_th,p_ph,m,L)+cfepx(p_th,p_ph,b,L)+cfepx(p_th,p_ph,l,L));
    return ep3;
}
