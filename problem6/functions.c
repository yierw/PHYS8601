#include "headers.h"
#include<math.h>
#include<stdlib.h>

double fcv1( int * E, double * P, double k, int N, int length)
{
    double cv,e2bar,ebar;
    int i;
    
    ebar=0.0;e2bar=0.0;
    
    for(i=0;i<=length;i++)
    {
        ebar+=(double) E[i] * P[i];
        e2bar+=(double) E[i] * (double) E[i] * P[i];
    }
    
    cv=(e2bar-ebar*ebar)*k*k/(double) N;
    
    return cv;
}

void reweight(FILE * fp3, FILE * fp4, int * pointer,double t0, int len,double * k,int klen, int N)
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
        fprintf(fp3,"# reweighting at T = %f\n",1.0/k[j]);
        for(i=0;i<=count;i++){fprintf(fp3," %4d  %10d  %10d  %10.6f\n",i,x[i],num[i],px[i]);}
        fprintf(fp4," %10.4f  %lf \n",1.0/k[j],fcv1(x,px,k[j],N,size));
    }
}

void correlation(FILE *fp1, int * earr, int * marr, int clen, int len)
{
    double ebar,e2bar,mbar,m2bar,ecor,mcor,psie[clen+1],psim[clen+1];
    int ts,i;
    ecor=0.0;mcor=0.0;
    ebar=0.0;e2bar=0.0;
    mbar=0.0;m2bar=0.0;
    for(ts=0;ts<=clen;ts++)
    {
        for(i=0;i<len;i++)
        {
            ecor+=(*(earr+i)) * (*(earr+i+ts));
            ebar+=*(earr+i);
            e2bar+=(*(earr+i))*(*(earr+i));
            mcor+=(*(marr+i)) * (*(marr+i+ts));
            mbar+=*(marr+i);
            m2bar+=(*(marr+i)) *(*(marr+i)) ;
        }
        ecor/=len;ebar/=len;e2bar/=len;
        mcor/=len;mbar/=len;m2bar/=len;
        psie[ts]=(ecor-ebar*ebar)/(e2bar-ebar*ebar);
        psim[ts]=(mcor-mbar*mbar)/(m2bar-mbar*mbar);
        fprintf(fp1," %4d %10.6f %10.6f \n",ts,psie[ts],psim[ts]);
    }
    
}

void thermal1(FILE *fp5, int * E, int * M,double t, int N, int length)
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
        ebar+=(double) E[i] ;
        e2bar+=(double) E[i] * (double) E[i];
        mbar+=(double) M[i];
        m2bar+=(double) M[i] * (double) M[i];
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

void thermal2(FILE *fp5, double * slng, int * E, double t, int N, int length)
{
    double cv,e2bar,ebar,k;
    double sum,tmp,mean,P[length];
    int i;
    
    if(t<0){k=0.0;/* infinite temperature*/}
    else if (t==0.0){printf("cannot calculate zero temperature\n");exit(1);}
    else{k=1.0/t;}
    /*find mean*/
    sum=0.0;
    for(i=0;i<length;i++)
    {
        tmp=*(slng+i) ;
        sum+=tmp;
    }
    mean=sum/length;
    /* find weighted probability*/
    sum=0.0;
    for(i=0;i<length;i++)
    {
        tmp=exp( *(slng+i) - (double) (*(E+i)) * k -mean);/*subtract mean to avoid too big values after take exp() */
        P[i]=tmp;
        sum+=tmp;
    }
    /* normalize */
    for(i=0;i<length;i++){tmp=P[i];P[i]=tmp/sum;}
    
    ebar=0.0;e2bar=0.0;
    
    for(i=0;i<length;i++)
    {/* unscaled by lattice size */
        fprintf(fp5," %d  %d %10.6f\n", i+1, E[i],P[i]);
        ebar+=(double) E[i] * P[i];
        e2bar+=(double) E[i] * (double) E[i] * P[i];
    }
    cv= (e2bar-ebar*ebar)*k*k;
    /* scaled by lattice size */
    ebar/=(double) N;
    cv/=(double) N;
    
    fprintf(fp5," %10.6f %10.6f %10.6f\n", t, ebar, cv);
}



