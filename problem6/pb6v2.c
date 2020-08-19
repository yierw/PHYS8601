#include "headers.h"
#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>

int main()
{
    int type=2; /* =1 metropolis; =2 wang-landau */
    int L=16; /* lattice size */
    int MCS=500000; /* total sampling steps */
    int burn=5000; /* burn in samplings */
    int clen=0; /* for correlation function calculation; if 0 then not perform */
    double t_grid[]={2.3,5.0};/* do sampling at this temperature (-1 means infinity) */
    double tk[]={};/* reweighting temperatures ; if empty then not perform */
    float p=0.5; /* choose p percent spins up initially */
    int printl=0; /* if want to print the spin configuration, set it to 1 */
    /* histogram and density function*/
    int size=300;
    int count,el,dep,s,l;
    int x[size],num[size];
    double slng[size],lng[size],lng1,lng2,lnf,r,tmp;
    float dif;
    /* !!do not change below!! */
    int klen = (int) sizeof(tk)/ (int) sizeof(double);
    int tlen = (int) sizeof(t_grid)/ (int) sizeof(double);
    int len=MCS-burn;
    double k[klen],t;
    int i,j,h;
    int N=L*L;
    int lattice[L][L]; /* L*L lattice; 1 for spin up; -1 for spin down */
    int * ptl = &lattice[0][0];
    int array[2];/* array[0] unscaled energy;array[1]=unscaled mag */
    int earr[len+clen],marr[len+clen];/* store unscaled energy and unscaled mag at each MCS */
    FILE * fp, * fp1, * fp2, * fp3, * fp4, *fp5;    /* output files */

    /* srand(time(NULL)); */
    
    if(type==1)/* metropolis */
    {
        printf("Metropolis Sampling\n");
        fp1=fopen("correlation.dat","w");
        fp2=fopen("output.dat","w");
        fp3=fopen("histogram.dat","w");
        fp4=fopen("reweightingcv.dat","w");
        fp5=fopen("thermal.dat","w");
        fprintf(fp5,"#      T         E         Cv         M           Chi\n");
        for (h=0;h<tlen;h++)/* START: the loop of different temperatures */
        {
            t=t_grid[h];printf("Do Sampling at T = %f\n",t);
            fprintf(fp2,"#  MCS    scaled E  scaled M ( T = %f , L = %d )\n",t,L);
            fprintf(fp3,"#  Draw simulations at       ( T = %f , L = %d )\n",t,L);

            initial(ptl,p,N);
            array[0]=fener(ptl,L);array[1]=fmag(ptl,N);
            
            for(i=0;i<(MCS+clen);i++)/* START: metropolis sampling */
            {
                for(j=0;j<N;j++){metropolis(ptl, t, array, L);}
                fprintf(fp2,"%4d %10.4f %10.4f \n",i,(float)array[0]/N,(float)abs(array[1])/N);
                if(i>=burn){earr[i-burn]=array[0];marr[i-burn]=abs(array[1]);}
            }
            
            if(clen>0)/* calculate correlation function */
            {
                fprintf(fp1,"#  MCS    psi_e        psi_m ( T = %f , L = %d )\n",t,L);
                correlation(fp1,earr,marr,clen,len);
            }
            if(klen>0)/* histogram reweighting */
            {
                fprintf(fp4,"#       T    Cv              ( T = %f , L = %d )\n",t,L);
                for(i=0;i<klen;i++){k[i]=1.0/tk[i];printf("Do reweighting at T = %f\n",tk[i]);}
                reweight(fp3,fp4,earr,t,len,k,klen,N);
            }
            thermal1(fp5,earr,marr,t,N,len); /* START: thermal properties */
        }
        
    }/* END: type 1 */
    else if(type==2)/* wang-landau  */
    {
        printf("Wang-Landau Sampling\n");
        fp3=fopen("histogram.dat","w");
        fp5=fopen("thermal.dat","w");
        fprintf(fp5,"#      T         E         Cv\n");
        
        initial(ptl,p,N);
        array[0]=fener(ptl,L);array[1]=fmag(ptl,N);
        
        for(l=0;l<size;l++){x[l]=0;num[l]=0;lng[l]=0.0;}
        lnf=1.0;lng1=0.0;
        count=0;
     
        for(i=0;i<MCS;i++)/* START: wang-landau sampling */
        {
            for(j=0;j<N;j++)
            {
                s=rand()% N;
                dep=flip(ptl,array,s,L);  /* flip spin return energy change */
                count=hist(array[0], x, num, lng, count);
                if(count>=size){ printf("increase size\n");exit(1); }/*count is the length of x */
                for(l=0;l<count;l++)
                {
                    if(x[l]==array[0]){lng2=lng[l];}
                }
                
                r=(rand()%1000000)/1000000.0;
                if(r < exp(lng1-lng2) ){/* accept this flip*/}
                else{dep=flip(ptl,array,s,L); /*reject, so flip it back */}
                
                for(l=0;l<count;l++)
                {/* re-search to make sure we use the current E=array[0]*/
                    if(x[l]==array[0])
                    {
                        lng1=lng[l]+lnf;/* update g1 */
                        lng[l]=lng1;/* update density anyway */
                    }
                }
            }/* END:one MCS  */

            if((i+1)%1000==0)
            {
                dif=diff(num,count);
                if(dif>=0.7)
                {
                    fprintf(fp3,"histogram of MCS %d: (current E=%d, M=%d)\n", i+1,array[0],array[1]);
                    for(l=0;l<count;l++){
                        slng[l]=lng[l]-min(lng,count);
                        fprintf(fp3," %4d  %f   %10d   %12.6f\n",l+1,(double)x[l]/(double)N,num[l],slng[l]);
                    }
                    lnf*=0.5;
                    for(l=0;l<count;l++){num[l]=0;}
                    printf("%f percent; reset lnf to %lf and zero H(E) at MCS :%d\n",dif,lnf,i+1);
                }
            }
        }/* END: wang-lamdau sampling*/
        for (h=0;h<tlen;h++)
        {
            t=t_grid[h];
            thermal2(fp5,slng,x,t,N,count);
        }
    }/* END: type 2*/
    
        fp=fopen("lattice.dat","w");
        wrtl(fp,ptl,printl,L);

    return 0;
}


