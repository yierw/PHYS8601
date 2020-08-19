#include "headers.h"
#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>

int main()
{
    int L=12; /* lattice size */
    int N=L*L;
    int flag=2;/* = 1 for quantum spin; = 2 for classical spin; = 3 for quantum like classical spin */
    int MCS=22000; /* total sampling steps */
    int burn=2000; /* burn in samplings */
    int len=MCS-burn;
    double earr[len],marr[len];/* store unscaled energy and unscaled mag at each MCS */
    double t;
    double t_grid[]={
        0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,
        1.21,1.22,1.23,1.24,1.25,1.26,1.27,1.28,1.29,
        1.3,1.4,1.5,1.6,1.7,1.8,1.9,
        2.0,2.5,3.0,3.5,4.0,4.5,5.0
    };/* do sampling at this temperature (-1 means infinity) */
    int tlen = (int) sizeof(t_grid)/ (int) sizeof(double);
    float p=1.0; /* choose p percent spins up initially */
    int printl=0; /* if want to print the spin configuration, set it to 1 */
    int i,j,h;
    /* for classical spins*/
    double th[L][L],ph[L][L];
    double * p_th = &th[0][0];
    double * p_ph = &ph[0][0];
    double carr[4],nm;
    double ENER[len],MAGN[len];
    /* for quantum spin */
    int lattice[L][L]; /* L*L lattice; 1 for spin up; -1 for spin down */
    int * ptl = &lattice[0][0];
    int array[2];/* array[0] unscaled energy;array[1]=unscaled mag */
    
    FILE * fp, * fp1, * fp2, * fp3, * fp4, *fp5;    /* output files */
    srand(time(NULL));
    if(flag == 1)
    {
        printf("quantum spin\n");
        fp2=fopen("output.dat","w");
        fp5=fopen("thermal.dat","w");
        fprintf(fp5,"#      T         E         Cv         M           Chi\n");
        
        for (h=0;h<tlen;h++)/* START: the loop of different temperatures */
        {
            t=t_grid[h];printf("Do Sampling at T = %f\n",t);
            fprintf(fp2,"#  MCS    scaled E  scaled M ( T = %f , L = %d )\n",t,L);
            
            initial(ptl,p,N);
            array[0]=fener(ptl,L);array[1]=fmag(ptl,N);
            
            for(i=0;i<MCS;i++)/* START: metropolis sampling */
            {
                for(j=0;j<N;j++){metropolis(ptl, t, array, L);}
                fprintf(fp2,"%4d %10.4f %10.4f \n",i,(float)array[0]/N,(float)abs(array[1])/N);
                if(i>=burn){earr[i-burn]=(double)array[0];marr[i-burn]=(double)abs(array[1]);}
            }
            
            thermal1(fp5,earr,marr,t,N,len); /* START: thermal properties */
        }
        
        fp=fopen("lattice.dat","w");
        wrtl(fp,ptl,printl,L);
    }
    else
    {
        /*fp2=fopen("output.dat","w");*/
        fp5=fopen("thermal.dat","w");
        fprintf(fp5,"#      T         E         Cv         M           Chi\n");

        for (h=0;h<tlen;h++)/* START: the loop of different temperatures */
        {
            t=t_grid[h];printf("Do Sampling at T = %f\n",t);
           /* fprintf(fp2,"#  MCS    scaled E  scaled M ( T = %f , L = %d )\n",t,L);*/
            
            for(i=0;i<L;i++)
            {
                for(j=0;j<L;j++)
                {
                    th[i][j]=0.0;ph[i][j]=0.0;        /* set all spins up */
                    /*printf("%f %f\n",th[i][j], ph[i][j]);*/
                }
            }
            iarr(p_th,p_ph,carr,L);
            nm=sqrt(carr[1]*carr[1]+carr[2]*carr[2]+carr[3]*carr[3]);

            for(i=0;i<MCS;i++)/* START: metropolis sampling */
            {
                for(j=0;j<N;j++){cmetropolis(p_th,p_ph, t, carr, L,flag);}
                nm=sqrt(carr[1]*carr[1]+carr[2]*carr[2]+carr[3]*carr[3]);
                /*fprintf(fp2,"%4d %10.4f %10.4f \n",i,carr[0]/N,nm/N);*/
                if(i>=burn){earr[i-burn]=carr[0];marr[i-burn]=nm;}
            }
            thermal1(fp5,earr,marr,t,N,len); /* START: thermal properties */
        }
    }
    return 0;
}
