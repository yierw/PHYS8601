#include "headers.h"
#include<stdlib.h>

void initial(int * ptl, float p, int N)
{
    int nup, rup;
    int a, i, j;
    
    nup = p*N;
    /* firstly all spins down */
    for(i=0;i<N;i++){*(ptl+i)=-1;}
    rup=0;
    
    while(rup<nup)
    {
        a=rand()%N;
        if(*(ptl+a)==-1)
        {
            *(ptl+a)=1;
            rup++;
        }
    }
}

void wrtl (FILE * fp, int *ptl, int printl, int L)
{
    int i,j;
    
    if(printl==1){
        for(i=0;i<L;i++)
        {
            for(j=0;j<L;j++){fprintf(fp,"  %4d  ",*(ptl+i*L+j));}
            fprintf(fp,"\n");
        }
    }
    else{
        fprintf(fp,"lattice unprinted");
        fprintf(fp," (set printl = 1 to print the lattice)\n");
    }
    
    
}

int fmag(int *ptl, int N)
{
    int i;
    int mag;
    
    mag=0;
    for(i=0;i<N;i++)
    {
        mag+=*(ptl+i);
    }
    return mag;
}

int fener(int *ptl, int L)
{
    int i,ener;
    int N=L*L;
    
    ener=0;
    for(i=0;i<N;i++)
    {
        ener+=fepx(ptl,i,L);
    }
    return ener;
}


int fepx (int *ptl, int x, int L)
{
    int m,a,r;/*position in the lattice*/
    int spm,spa,spr;
    int epx;
    int N=L*L;
    int col=x%L;
    int row=(x-col)/L;
    
    /* middle site */
    m=x;
    /* above site */
    if(row==0){a=x-L+N;}
    else{a=x-L;}
    /* right site */
    if(col==(L-1)){r=x+1-L;}
    else{r=x+1;}
    
    spm=*(ptl+m);
    spa=*(ptl+a);
    spr=*(ptl+r);
    epx=-spm*spa-spm*spr;
    
    return epx;
}

int fep3(int *ptl, int x, int L)
{
    int m,b,l;/*position in the lattice*/
    int ep3;
    int N=L*L;
    int col=x%L;
    int row=(x-col)/L;
    
    /* middle site */
    m=x;
    /* below site */
    if(row==(L-1)){b=x+L-N;}
    else{b=x+L;}
    /*left site */
    if(col==0){l=x-1+L;}
    else{l=x-1;}
    
    ep3=(fepx(ptl,m,L)+fepx(ptl,b,L)+fepx(ptl,l,L));
    
    return ep3;
}

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

/* inserting sorting */
void insert_sort(int * a, int * b, int length)
{
    
    int i,j,m;
    int key,key2;
    
    for(j=1;j<length;j++)
    {
        key=*(a+j);
        key2=*(b+j);
        i=j-1;
        while(i>-1&&*(a+i)>key){
            *(a+i+1)=*(a+i);
            *(b+i+1)=*(b+i);
            i--;
        }
        *(a+i+1)=key;
        *(b+i+1)=key2;
    }
}

/* inserting sorting */
void insert_sort3(int * a, int * b, double * c, int length)
{
    
    int i,j,m;
    int key,key2;
    double key3;
    
    for(j=1;j<length;j++)
    {
        key=*(a+j);
        key2=*(b+j);
        key3=*(c+j);
        i=j-1;
        while(i>-1&&*(a+i)>key){
            *(a+i+1)=*(a+i);
            *(b+i+1)=*(b+i);
            *(c+i+1)=*(c+i);
            i--;
        }
        *(a+i+1)=key;
        *(b+i+1)=key2;
        *(c+i+1)=key3;
    }
}

double min(double * x,  int length)
{
    
    int i,j;
    double a[length];
    double key,min;
    
    /* copy the array */
    for(i=0;i<length;i++)
    {
        a[i]=*(x+i);
    }
    
    /* sort the array */
    for(j=1;j<length;j++)
    {
        key=*(a+j);
        i=j-1;
        while(i>-1&&*(a+i)>key){
            *(a+i+1)=*(a+i);
            i--;
        }
        *(a+i+1)=key;
    }
    
    return a[0];
}






