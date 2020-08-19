#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
/* two ways array[i][j]  or *(ptr+i*L+j) */
#define L 50   /* length of the lattice */
#define N L*L  /* number of sites in the lattice */
#define p 0.3  /* concentration of the sites */
#define sim 1  /* number of simulations */
#define flag 1 /* print out detailed results (=1) or not (=0)*/

void rocp(int *ptl, FILE * fp ); /* function to randomly occupy sites */
int lab(int * ptl, int * pts, FILE * fp );/* function to label sites */
void comb(int index, int olab, int rlab, int * ptl);/* function to combine clusters; called by lab */


int findvert(int ncluster, int * ptl,FILE * fp);
int findhori(int ncluster, int * ptl,FILE * fp);
int compare(int * arr, int n, int x);



int main()
{
    int label[L][L]; /* occupation/labels of sites in the lattice*/
    int size[N/2];
    int * ptl = &label[0][0]; /* pointer to the first element of label */
    int * pts = &size[0];
    int ncluster;
    int infv,infh,inf;

    
    int i,j,k;
    float tot=sim;
    FILE * fp;
   
    srand(time(NULL));/* use the current time as seed for random number generator */

    fp=fopen("out.dat","w");/* name of the output file*/

    inf=0;
    for(k=0;k<sim;k++)
    {
        /* initialize*/
        for(i=0;i<N;i++){*(ptl+i)=0;}
        for(i=0;i<N/2;i++){*(pts+i)=0;}
        
        rocp(ptl,fp);  /* randomly occupy the sites in the lattice with free edge boundary conditions*/
        ncluster=lab(ptl, pts,fp); /* label the sites */
        
        infv=findvert(ncluster,ptl,fp);/*search vertical infinite cluster*/
        if(infv==1){inf++;continue;}/*if vertical infinite cluster found, then skip horizontal search*/
        else
        {
            infh=findhori(ncluster,ptl,fp);/*search horizontal infinite cluster*/
            inf+=infh;
        }
        
        printf("%d th sim done\n",k);
    }

    printf("probability of having an infinite cluster is %10.4f\n",(float)inf/tot);
    fclose(fp);

    return 0;
}

void rocp(int * ptl, FILE * fp)
{
    int i,j,k,a,tot;
    float u,prel,noc;

    tot=N;
    noc = p*N;
    /*printf("%d %f\n",tot,noc); */

    if(flag)
    {
        fprintf(fp,"initial values of sites\n");
        for(i=0;i<L;i++){
            for(j=0;j<L;j++){fprintf(fp,"  %4d  ",*(ptl+i*L+j));}
            fprintf(fp,"\n");
        }
    }
    /* generate pseudo-random number between 0 and 1 for each site;
     compare the random number with the concentration probability p;
     if smaller, occupy this site; if bigger, not cooupy*/
    for(i=0;i<N;i++){
        u=(rand()%1000000)/1000000.0;
        if(u<=p){*(ptl+i)=1;}
    }
   
    k=0;
    for(i=0;i<N;i++){k+=*(ptl+i);}
 
    prel=(float)k/(float)tot;
    if(flag){fprintf(fp,"real occupation probability is %10.4f\n",prel);}

    while(k>noc)/* need to unoccupy*/
    {
        a=rand()%tot; /* generate random number from zero to tot=N=L*L */
        if(*(ptl+a)==1){*(ptl+a)=0;k--;}
        prel=(float)k/(float)tot;
        if(flag){fprintf(fp,"real occupation probability is %10.4f\n",prel);}
    }
    while(k<noc)/* need to occupy */
    {
        a=rand()%tot; /* generate random number from zero to tot=N=L*L */
        if(*(ptl+a)==0){*(ptl+a)=1;k++;}
        prel=(float)k/(float)tot;
        if(flag){fprintf(fp,"real occupation probability is %10.4f\n",prel);}
    }

    if(flag)
    {
        fprintf(fp,"occupied lattice sites\n");
        for(i=0;i<L;i++)
        {
            for(j=0;j<L;j++){fprintf(fp,"  %4d  ",*(ptl+i*L+j));}
            fprintf(fp,"\n");
        }
    }
    
}

int lab(int * ptl, int * pts, FILE * fp)
{
    int max,left,above;
    int i,j,tmp;
    int index;
    int count;
    
    max=0;
    for(i=0;i<L;i++)
    {
        for(j=0;j<L;j++)
        {
            if(i==0){above=0;}
            else{above=*(ptl+(i-1)*L+j);}
            if(j==0){left=0;}
            else{left=*(ptl+i*L+j-1);}
            if(*(ptl+i*L+j)==1)
            {
                if(left==0&&above==0)
                {
                    max++;
                    *(ptl+i*L+j)=max;
                    tmp=*(pts+max);*(pts+max)=tmp+1;
                }
                else if(above!=0&&left==0)
                {
                    *(ptl+i*L+j)=above;
                    tmp=*(pts+above);*(pts+above)=tmp+1;
                    
                }
                else if(above==0&&left!=0)
                {
                    *(ptl+i*L+j)=left;
                    tmp=*(pts+left);*(pts+left)=tmp+1;
                    
                }
                else if(above==left)
                {
                    *(ptl+i*L+j)=left;
                    tmp=*(pts+left);*(pts+left)=tmp+1;
                }
                else if(above>left)
                {
                    *(ptl+i*L+j)=left;/*label the current site according to left neighbor*/

                    tmp=*(pts+left);
                    *(pts+left)=*(pts+above)+tmp+1;
                    *(pts+above)=0;
                    
                    index=i*L+j;
                    comb(index,above,left,ptl);
                }
                else if(left>above)
                {
                    *(ptl+i*L+j)=above;/*label the current site according to above neighbor*/
                    
                    tmp=*(pts+above);
                    *(pts+above)=*(pts+left)+tmp+1;
                    *(pts+left)=0;
                    
                    index=i*L+j;
                    comb(index,left,above,ptl);
                }
            }
            else{*(ptl+i*L+j)=0;}
        }
    }
    
    if(flag)
    {
        fprintf(fp,"label of occupied sites\n");
        for(i=0;i<L;i++)
        {
            for(j=0;j<L;j++){fprintf(fp,"  %4d  ",*(ptl+i*L+j));}
            fprintf(fp,"\n");
        }
        
        fprintf(fp,"label size\n");
        count=0;
        for(i=0;i<N/2;i++)
        {
            if(*(pts+i)>0)
            {
                fprintf(fp,"%4d  %4d  \n",i,*(pts+i));
                count++;
            }
        }
        fprintf(fp,"%d clusters\n",count);
    }
    else
    {
        count=0;
        for(i=0;i<N/2;i++)
        {
            if(*(pts+i)>0){count++;}
        }
    }
    
    return count;
    
}

void comb(int index, int olab, int rlab, int * ptl)
{
    int k;

    for(k=0;k<index;k++)
    {
        if(*(ptl+k)==olab){*(ptl+k)=rlab;}
    }
    
}

int findvert(int ncluster, int * ptl,FILE * fp)
{
    int arrv1[ncluster];
    int arrv2[ncluster];
    int n1,n2;
    int i,j,ind,x;
    int exist;
    
    n1=0;n2=0;
    for(i=0;i<ncluster;i++){arrv1[i]=-1;}
    for(i=0;i<ncluster;i++){arrv2[i]=-1;}


    for(i=0;i<L;i++)
    {
        x=*(ptl+i);
        if(x==0){continue;}
        ind=compare(arrv1,n1,x);
        if(ind==-1)
        {
            arrv1[n1]=x;
            n1++;
        }
    }
    
    
    for(i=0;i<L;i++)
    {
        x=*(ptl+L*L-L+i);
        if(x==0){continue;}
        ind=compare(arrv2,n2,x);
        if(ind==-1)
        {
            arrv2[n2]=x;
            n2++;
        }
    }
    
    exist=0;
     for(i=0;i<n1;i++)
     {
         ind=compare(arrv2,n2,arrv1[i]);
         if(ind!=-1)
         {
             if(flag){fprintf(fp,"vertical infinite cluster label:%d\n",arrv1[i]);}
             exist++;
         }
         
     }
    
    if(exist>0){return 1;}
    else{return 0;}
}

int findhori(int ncluster, int * ptl,FILE * fp)
{
    int arrv1[ncluster];
    int arrv2[ncluster];
    int n1,n2;
    int i,j,ind,x;
    int exist;
    
    n1=0;n2=0;
    for(i=0;i<ncluster;i++){arrv1[i]=-1;}
    for(i=0;i<ncluster;i++){arrv2[i]=-1;}
    
    
    for(i=0;i<L;i++)
    {
        x=*(ptl+i*L);
        if(x==0){continue;}
        ind=compare(arrv1,n1,x);
        if(ind==-1)
        {
            arrv1[n1]=x;
            n1++;
        }
    }
    
    
    for(i=0;i<L;i++)
    {
        x=*(ptl+i*L+L-1);
        if(x==0){continue;}
        ind=compare(arrv2,n2,x);
        if(ind==-1)
        {
            arrv2[n2]=x;
            n2++;
        }
    }
    
    exist=0;
    for(i=0;i<n1;i++)
    {
        ind=compare(arrv2,n2,arrv1[i]);
        if(ind!=-1)
        {
            if(flag){fprintf(fp,"horizontal infinite cluster label:%d\n",arrv1[i]);}
            exist++;
        }
        
    }
    
    if(exist>0){return 1;}
    else{return 0;}
}


int compare(int * arr, int n, int x)
{
    int i;
    
    for(i=0;i<n;i++)
    {
        if(x==*(arr+i)){return i;}
    }
    return -1;
    
}


