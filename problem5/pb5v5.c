#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>

void initial(int * ptl, float p, int N); /* function to randomly occupy sites */
void wrtl (FILE * fp, int *ptl, int printl, int L); /*write lattice occupation*/
void fspins(int *ptl, double t, int * array, int L);
void mkhist(FILE * fp3,FILE * fp4, int * pointer,double t, int len, double * k,int klen, int N);

int fmag(int *ptl, int N);
int fepx (int *ptl, int x, int L);/* calculate the nearest-neighbor interaction energy of site (row,col) periodic boundary conditions applied*/
int fep3(int *ptl, int x, int L);
int fdep (int *ptl, int x, int L);/*return energy change per site*/
double fcv1( int * E, double * P, double k, int N, int length);
double fcv2( int * E, double k, int N, int length);
void insert_sort(int * a, int * b, int length);

int main()
{
    int L=8; /* lattice size */
    int MCS=50000; /* total sampling steps */
    int burn=2000; /* throw away first burn MCS runs*/
    
    int clen=0; /* time step until tlen for correlation study; if 0 then not perform */
    double tk[]={1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,
        2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,
        3.0,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,
        4.0,4.1,4.2,4.3,4.4,4.5,4.6,4.7,4.8,4.9,
        5.0};/* reweighting temperatures ; if empty then not perform */
   double t_grid[]={-1};/*do sampling at this temperature */
    
    float p=0.3; /* choose p percent spins up initially */
    int printl=0; /* if want to print the spin configuration, set it to 1 */
    /* !!do not change below!! */
    int klen = (int) sizeof(tk)/ (int) sizeof(double);
    int tlen = (int) sizeof(t_grid)/ (int) sizeof(double);
    double k[klen],t;
    int ts,i,j,h;
    int N=L*L;
    int len=MCS-burn;
    int lattice[L][L]; /* L*L lattice; 1 for spin up; -1 for spin down */
    int * ptl = &lattice[0][0];
    int array[2];/* array[0] unscaled energy;array[1]=unscaled mag */
    int earr[len+clen],marr[len+clen];/* store unscaled energy and unscaled mag at each MCS */
    double ebar,e2bar,mbar,m2bar,ecor,mcor,cv,psie[clen+1],psim[clen+1];
    
    /* output files */
    FILE * fp, * fp1, * fp2, * fp3, * fp4, *fp5;
    fp=fopen("lattice.dat","w");
    fp1=fopen("correlation.dat","w");
    fp2=fopen("output.dat","w");
    fp3=fopen("histogram.dat","w");
    fp4=fopen("reweightingcv.dat","w");
    fp5=fopen("samplingcv.dat","w");
    fprintf(fp5,"#       T    Cv ( do sampling with L=%d)\n",L);
    
    /* initialize the lattice*/
    srand(time(NULL));
    initial(ptl,p,N);
    array[0]=0.0;for(i=0;i<N;i++){array[0]+=fepx(ptl,i,L);}
    array[1]=fmag(ptl,N);
    /* start sampling */
    for (h=0;h<tlen;h++)
    {
        t=t_grid[h];
        fprintf(fp1,"# t(MCS)   psi_e       psi_m ( at T=%f with L=%d )\n",t,L);
        fprintf(fp2,"#  MCS    scaled E  scaled M ( at T=%f with L=%d and burn in first %d MCS \n",t,L,burn);
        fprintf(fp3,"#label  unscaled E  counts  probability( at T=%f with L=%d )\n",t,L);
        fprintf(fp4,"#       T    Cv ( do sampling at T=%f with L=%d)\n",t,L);

        for(i=0;i<burn;i++)  /* burn in */
        {
            for(j=0;j<N;j++){fspins(ptl, t, array, L);}
            fprintf(fp2,"%4d %10.4f %10.4f \n",i+1,(float)array[0]/N,(float)abs(array[1])/N);
        }
      
        for(i=0;i<(len+clen);i++)  /* keep*/
        {
            for(j=0;j<N;j++){fspins(ptl, t, array, L);}
            fprintf(fp2,"%4d %10.4f %10.4f \n",i+1+burn,(float)array[0]/N,(float)abs(array[1])/N);
            earr[i]=array[0];marr[i]=abs(array[1]);
        }
        
        wrtl(fp,ptl,printl,L);
        
        if(t>0)
        {
            fprintf(fp4,"#%10.4f  %lf [importance sampling  ]\n",t,fcv2(earr,1.0/t,N,len));
            fprintf(fp5,"%10.4f  %lf \n",t,fcv2(earr,1.0/t,N,len));
        }
        else
        {
            fprintf(fp4,"#%10.4f  %lf [importance sampling  ]\n",t,fcv2(earr,0.0,N,len));
            fprintf(fp5,"%10.4f  %lf \n",t,fcv2(earr,0.0,N,len));
        }
        
        if(clen>0)/* correlation function calculating */
        {/* initial value for each physical quantity*/
            ecor=0.0;mcor=0.0;
            ebar=0.0;e2bar=0.0;
            mbar=0.0;m2bar=0.0;
            for(ts=0;ts<=clen;ts++)
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
        }

        
        if(klen>0)/* histogram reweighting */
        {
            for(i=0;i<klen;i++){k[i]=1.0/tk[i];}
            mkhist(fp3,fp4,earr,t,len,k,klen,N);
        }

        
    }
    
 
    return 0;
}

void mkhist(FILE * fp3,FILE * fp4,int * pointer,double t0, int len,double * k,int klen, int N)
{
    int i,j,flag,count;
    int size=100;
    int x[size],num[size];
    double px[size];
    double k0,dk,sum;

    for(i=0;i<size;i++){x[i]=0.0;px[i]=0.0;num[i]=0;}
    count=0;
    for(i=0;i<len;i++)
    {
        flag=1;/*preset to record the energy*/
        for(j=0;j<=count;j++)
        {
            if(*(pointer+i)==x[j])
            {
                flag=0;/*find repeated energy, so not record*/
                num[j]++;
            }
        }
        if(flag==1)/*only record energy with flag*/
        {
            count++;
            if(count>=size){printf("increase size in mkhist\n");exit(1);}
            x[count]=*(pointer+i);
            num[count]++;
        }
    }
    insert_sort(x,num,count+1);/*sort the histogram*/
    /* calculate probability distribution*/
    if(t0<0){k0=0.0;}else{k0=1.0/t0;}
    for(j=0;j<klen;j++)
    {

        if(k[j]<0){k[j]=0.0;}
        dk=k0-k[j];
        sum=0.0;
        for(i=0;i<=count;i++)
        {
            px[i]=(double) num[i]*exp(dk*x[i]);
            sum+=px[i];
        }
        for(i=0;i<=count;i++){px[i]/=sum;}
        fprintf(fp4," %10.4f  %lf \n",1.0/k[j],fcv1(x,px,k[j],N,size));
        fprintf(fp3,"# histogram reweighting at %f\n",1.0/k[j]);
        for(i=0;i<=count;i++){fprintf(fp3," %4d  %10d  %10d  %10.6f\n",i,x[i],num[i],px[i]);}
    }
}

double fcv1( int * E, double * P, double k, int N, int length)
{
    double cv,e2bar,ebar;
    int i;
    
    for(i=0;i<=length;i++)
    {
        ebar+=(double) E[i] * P[i];
        e2bar+=(double) E[i] * (double) E[i] * P[i];
    }
    
    cv=(e2bar-ebar*ebar)*k*k/(double) N;
    return cv;
}

double fcv2( int * E, double k, int N, int length)
{
    double cv,e2bar,ebar,tmp;
    int i;
    
    for(i=0;i<length;i++)
    {
        ebar+=(double) E[i];
        e2bar+=(double) E[i]*(double) E[i];
    }
    ebar/= (double) length;
    e2bar/= (double) length;

    cv=(e2bar-ebar*ebar)*k*k/(double) N;
    return cv;
}

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



void fspins(int *ptl, double t, int * array, int L)
{
    double r;
    int dep,tmp,ep,mag;
    int x,spin; /* position and value of the current spin */
    int N=L*L;
    
    ep=*array;
    mag=*(array+1);
    
    x=rand() % N;
    dep=fdep(ptl,x,L);/* flip spin; calculate energy change*/
    
    if(t<0)/*infinite temperature; randomly generating*/
    {
        /* accept this flip; update physical quantities */
        spin=*(ptl+x);
        mag+=2.0*spin;
        ep+=dep;
    }
    else
    {
        r=(rand()%1000000)/1000000.0;
        if(r<exp(-dep/t)){
            /* accept this flip; update physical quantities */
            spin=*(ptl+x);
            mag+=2.0*spin;
            ep+=dep;
        }
        else{
            /* reject; flip x back; do not update physical quantities */
            tmp=*(ptl+x);*(ptl+x)=-tmp;
        }
        
    }
    
    *array=ep;
    *(array+1)=mag;
    
}

/* inserting sorting */
void insert_sort(int * a, int * b,int length)
{
    
    int i,j,m;
    int key,key2;

    for(j=1;j<length;j++)
    {
        key=*(a+j);key2=*(b+j);
        i=j-1;
        while(i>-1&&*(a+i)>key){
            *(a+i+1)=*(a+i);*(b+i+1)=*(b+i);
            i--;
        }
        *(a+i+1)=key;*(b+i+1)=key2;
    }
}

