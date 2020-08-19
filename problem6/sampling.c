#include "headers.h"
#include<math.h>
#include<stdlib.h>

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
        }
    }
    
}


int hist(int E, int * x, int * num, double * g, int count)
{
    int i,j,flag,tmp;
    /* accumulate the new energy in histogram */
    
    if(count==0){
        *x=E;
        *num=1;
        count++;
        return count;
    }

    flag=1;/* pre-set -- to record this energy*/
    
    for(j=0;j<count;j++)
    {
        if(E==*(x+j))
        {
            flag=0;/* find repeated energy, so not add new record*/
            tmp=*(num+j);
            *(num+j)=tmp+1;
        }
    }
    
    if(flag==1)/* only record energy with flag = 1 */
    {
        *(x+count)=E;
        tmp=*(num+count);
        *(num+count)=tmp+1;
        count++;
    }
    
    insert_sort3(x,num,g,count);/* sort the histogram along with the density function */
    return count;
}


double diff(int * x,  int length)
{
    
    int i,j,m,a[length];
    int key,mean,tmp;
    
    /* copy the array */
    mean=0.0;
    for(i=0;i<length;i++)
    {
        a[i]=*(x+i);
        tmp=mean;
        mean=tmp+a[i];
    }
    mean/=length;
    
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
    
    for(i=0;i<length;i++)
    {
        if(a[i]>=mean)
        {
            return 1.0-(double) i/length;
        }
    }
    
    return -1.0;
}


