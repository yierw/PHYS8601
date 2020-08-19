#include<stdio.h>
#include<math.h>

int main()
{
    FILE *fp;
    FILE *fo;
    int len=21;
    float x[len],y[len],f;
    int i;
    
    fp=fopen("L20.dat","r");
    fscanf(fp, "%*[^\n]");
    for(i=0;i<len;i++){
        fscanf(fp,"%f %f %*[^\n]",&x[i],&y[i]);
    }
    fclose(fp);

    fo=fopen("ndcv.dat","w");
    for(i=1;i<len-1;i++){
       f=(y[i+1]-y[i-1])/(x[i+1]-x[i-1]);
       fprintf(fo,"%f %f\n",x[i],f);
    }
    fclose(fo);
    
    return 0;
}
