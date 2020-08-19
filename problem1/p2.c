/* This code is used to evaluate the value of Pi 
by apply a simpling Monte Carlo method. 

This code is used for performing a simgle Monte Carlo simulation
with different random number streams for each run.

The number of trials included in one run is set up 
by choosing defferent value of length*/

#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>

#define length 10000
#define n 2 /* print out first n points*/

int main()
{
 double x,y;
 double r2,r,est,cons;
 int i,t1,t2;
 
 cons=4.0;
 srand(time(NULL));/* initialize the random number stream*/
 printf("\nSeed of this run:%u\n",time(0));
 printf("The first %d random generated point:\n",n);
 t1=0;
 t2=0;
 for(i=0;i<n;i++){
  x=(rand()%1000000)/1000000.0; /*generate numbers*/
  y=(rand()%1000000)/1000000.0;
  printf("(%10.6f, %10.6f)\n",x,y);
  r2=pow(x,2)+pow(y,2);
  if(r2<=1){t1++;}
  else{t2++;}
  }

 for(i=n;i<length;i++){
  x=(rand()%1000000)/1000000.0;
  y=(rand()%1000000)/1000000.0;
  r2=pow(x,2)+pow(y,2);
  if(r2<=1){t1++;}
  else{t2++;}
  }
 
  printf("%d simulations done.\n",length);
  est=cons*t1/(t1+t2);
  printf("Estimated pi:%10.6f\n",est);
  return 0;
}

