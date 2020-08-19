#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>

#define flag 0 /*print out lattice =1*/
#define t 3.0
#define L 10
#define MCS 5000
#define burn 2000 /* throw away first burn MCS runs*/

void initial(int *ptl, float p); /* function to randomly occupy sites */
void wrtl (FILE * fp, int *ptl); /*write lattice occupation*/
double fepx(int *ptl, int x);/* calculate the nearest-neighbor interaction energy of site (row,col) periodic boundary conditions applied*/
double fep3(int *ptl, int x);
double fdep(int *ptl, int x);/*return energy change per site*/
double fmag(int *ptl);
double simps3(double *y, double h, int size);

int main()
{
    float p=0.2; /* choose p percent spins up initially */
    int tlen=100;
    int len=MCS-burn;
    int ts;
    double earr[len+tlen],marr[len+tlen];
    double ebar,e2bar,mbar,m2bar,ecor,mcor;
    double psie[tlen+1],psim[tlen+1];
    double inte,intm;
/* do not change */
    int N=L*L;
    int lattice[L][L]; /* L*L lattice; 1 for spin up; -1 for spin down */
    int * ptl = &lattice[0][0];
    double dep,r;
    double ep,mag;/*scaled by lattice size*/
    int x,spin; /* position and value of the current spin */
    int i,j,tmp;
/* no change warning end*/
    FILE * fp;
    FILE * fp1;
    FILE * fp2;

   srand(time(NULL));
    initial(ptl,p);/* initialize the lattice*/
    mag=fmag(ptl)/N;
    ep=0.0;
    for(x=0;x<N;x++){ep+=fepx(ptl,x)/N;}
    
    fp1=fopen("out.dat","w");
    fprintf(fp1,"# L= %d  t= %f \n",L,t);
    fprintf(fp1,"# t(MCS)   psi_e       psi_m \n");

    fp2=fopen("out2.dat","w");
    fprintf(fp2,"#  %f  %d\n ",t,L);
    fprintf(fp2,"#  burn in first %d MCS\n",burn);
    fprintf(fp2,"#  MCS    scaled E  scaled M \n");

/* samplings to burn */
    for(i=0;i<burn;i++)
    {
        for(j=0;j<N;j++)
        {
            x=rand() % N;
            dep=fdep(ptl,x);/* flip spin; calculate energy change*/
            r=(rand()%1000000)/1000000.0;
            if(r<exp(-N*dep/t)){
                /* accept this flip; update physical quantities */
                spin=*(ptl+x);
                mag+=2.0*spin/N;
                ep+=dep;
            }
            else{
                /* reject; flip x back; do not update physical quantities */
                tmp=*(ptl+x);*(ptl+x)=-tmp;
            }
        }
        fprintf(fp2,"%4d %10.6f %10.6f \n",i+1,ep,fabs(mag));
    }
    
/* sampling keep*/
    for(i=0;i<(len+tlen);i++)
    {
        for(j=0;j<N;j++)
        {
            x=rand() % N;
            dep=fdep(ptl,x);/* flip spin; calculate energy change*/
            r=(rand()%1000000)/1000000.0;
            if(r<exp(-N*dep/t)){
                /* accept this flip; update physical quantities */
                spin=*(ptl+x);
                mag+=2.0*spin/N;
                ep+=dep;
            }
            else{
                /* reject; flip x back; do not update physical quantities */
                tmp=*(ptl+x);*(ptl+x)=-tmp;
            }
        }
        fprintf(fp2,"%4d %10.6f %10.6f \n",i+1+burn,ep,fabs(mag));
        earr[i]=ep;
        marr[i]=fabs(mag);
    }

    if(flag==1)
    {
        fp=fopen("lattice.dat","w");
        wrtl(fp,ptl);
    }
    
    /* initial value for each physical quantity*/
    ecor=0.0;mcor=0.0;
    ebar=0.0;e2bar=0.0;
    mbar=0.0;m2bar=0.0;
    
    for(ts=0;ts<=tlen;ts++)
    {
        for(i=0;i<=len;i++)
        {
            ecor+=earr[i]*earr[i+ts];
            ebar+=earr[i];
            e2bar+=earr[i]*earr[i];
            mcor+=marr[i]*marr[i+ts];
            mbar+=marr[i];
            m2bar+=marr[i]*marr[i];
        }
        ecor/=len;ebar/=len;e2bar/=len;
        mcor/=len;mbar/=len;m2bar/=len;
        psie[ts]=(ecor-ebar*ebar)/(e2bar-ebar*ebar);
        psim[ts]=(mcor-mbar*mbar)/(m2bar-mbar*mbar);
        fprintf(fp1," %4d %10.6f %10.6f \n",ts,psie[ts],psim[ts]);
    }
 
/*    printf("integrated correlation time\n");
    inte=simps3(psie,1.0,tlen);
    intm=simps3(psim,1.0,tlen);
    printf("for energy: %f\n",inte);
    printf("for magnetization: %f\n",intm);
    */
    return 0;
}


double simps3(double *y, double h, int size)
{
    int i;
    double f;

    f=0.0;
    for(i=0;i<size-2;i++)
    {
        f+=*(y+i);
        f+=*(y+i+1)*4.0;
        f+=*(y+i+2);
    }
    return f*h/3.0;
}


void initial(int * ptl, float p)
{
    int N, nup, rup;
    int a, i, j;
    
    N=L*L;
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

double fepx (int *ptl, int x)
{ 
    int m,a,r;/*position in the lattice*/
    int spm,spa,spr;
    double epx;
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

double fep3(int *ptl, int x)
{
  int m,b,l;/*position in the lattice*/
  double ep3;
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
  
  ep3=(fepx(ptl,m)+fepx(ptl,b)+fepx(ptl,l))/N;
  return ep3;
}

double fdep (int *ptl, int x)
{
   double epo,epf;
   double tmp,dep;

   /*energy before flipping*/
   epo=fep3(ptl,x);
   /* flip spin */
   tmp=*(ptl+x);
   *(ptl+x)=-tmp;
   /*energy after flipping*/
   epf=fep3(ptl,x);
   /*energy change in this flipping*/
   dep=epf-epo;
   return dep; 
}

void wrtl (FILE * fp, int *ptl)
{
    int i,j;

    for(i=0;i<L;i++)
    {
        for(j=0;j<L;j++){fprintf(fp,"  %4d  ",*(ptl+i*L+j));}
        fprintf(fp,"\n");
    }
}

double fmag(int *ptl)
{
    int i;
    int N=L*L;
    double mag;
    mag=0.0;
    for(i=0;i<N;i++)
    {
        mag+=*(ptl+i);
    }
    return mag;
}

