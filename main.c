//
//  main.c
//  Bratu-C
//
//  Created by Tim Ioannidis on 12/19/11.
//  Copyright 2011 Chemeng NTUA. All rights reserved.
//
    //Bratu 2D Parallel C implementation
    //NTUA, School of Chemical Engineering
    
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "parameters.h"
#include "gmres_config.h"
#include "extern.h"
#include "structs.h"
#include "methods_decl.h"
#include "cuda_wrapper.h"
#include "abfind.c"
#include "Back_sub.c"
#include "CSR.c"
#include "GMRES_m.c"
#include "Krylov.c"
#include "Least_sq.c"
#include "nodnumb.c"
#include "tsfun.c"
#include "xycoord.c"
#include "xydiscr.c"


//gia metrisi xronou

clock_t start1, end1, start2, end2;
double elapsed_GMRES;
    //global variable N->unknowns number definition
    int N=0,dim_buf=0,dim=0,Nz=0;
    //global
    struct common3 gelim;

int main (int argc,char* argv[])
{
    int indexAA[2]={0,0},*inner;
    int nell=0,i=0,j=0,p=0,q=0,iterat=0,k=0,l=0,init_el=0,fin_el=0;
    int *ncod;
    float *xpt,*ypt;
    int temp_div=0,temp_mod=0;
    double *r1,*u,*d;
    struct common1 mesh1;
    struct common2 elem1;
    struct common4 sparse;
    struct common5 param;
    
    nodnumb(&elem1);
    N=elem1.np;
    //dynamic allocation gia na min exw stack-overflow->heap
    //calloc kanei kai initialization me 0
    xpt=(float*)calloc(N,sizeof(double));
    ypt=(float*)calloc(N,sizeof(double));
    ncod=(int*)calloc(N,sizeof(int));
    r1=(double*)calloc(N,sizeof(double));
    u=(double*)calloc(N,sizeof(double));
    d=(double*)calloc(N,sizeof(double));
    //mesh construction
    xydiscr(&mesh1);
    xycoord(&mesh1,&elem1,xpt,ypt);
    //Parametroi provlimatos
    param.alfa=0;
    param.lamda=2;
    //orismos nzeros-sparsity
    sparse.nzeros=(long int)(sparsity*N*N);
    dim=N+1;
       //allocation of dynamic arrays
    sparse.AA=(double*)malloc(sparse.nzeros*sizeof(double));
    sparse.JA=(int*)malloc(sparse.nzeros*sizeof(int));
    sparse.IA=(int*)malloc(dim*sizeof(int));
        //Prepare for Dirichlet boundary conditions
        //markarisma aristera
    for (i=0;i<elem1.nny;i++){
         ncod[i]=1;
    }
        //markarisma katw
    for (i=0;i<(N-elem1.nny+1);(i+=elem1.nny)){
        ncod[i]=1;
    }
        //markarisma panw
    for (i=(elem1.nny-1);i<N;i+=elem1.nny){
        ncod[i]=1;
    }
       //markarisma dexia
    for (i=(N-elem1.nny);i<N;i++){
        ncod[i]=1;
    }
        //ksekinima NEWTON
        iterat=1;
    while (iterat<=Newton_iter){
        for (i=0;i<sparse.nzeros;i++){
            sparse.AA[i]=0;
            sparse.JA[i]=0;
        }
        for (i=0;i<dim;i++){
            sparse.IA[i]=i;
        }
        for (i=0;i<N;i++){
            r1[i]=0;
        }
        Nz=0;
        //MESH SCANNING
            //kathe process sarwnei ena kommati tou mesh->sugekrimena elements
        for (nell=0;nell<elem1.ne;nell++){
            abfind(nell,r1,u,&elem1,xpt,ypt,&sparse,&param,ncod);
        }
            //Impose essential boundary conditions gia ton pinaka r1
            //gia ton sk einai mesa stin CSR
        for (i=0; i<N; i++) {
            if (ncod[i]==1) {
                r1[i]=0;
                r1[i]=-(u[i]-0);
            }
        }
        if ((linear_solver[0]=='s') || (linear_solver[0]=='S')) {
            start1 = clock();
            GMRES_m(d,r1,&sparse);
            end1 = clock();
            elapsed_GMRES = ((double) (end1 - start1)) / CLOCKS_PER_SEC;
            if (iterat==1) {
                printf("\nXronos gia serial GMRES=%.5lfs\n",elapsed_GMRES);
            }
        }
        else if ((linear_solver[0]=='c') || (linear_solver[0]=='C')){
            cuda_GMRES(d,r1,&sparse);
        }
        else if ((linear_solver[0]=='b') || (linear_solver[0]=='B')){
            cuda_GMRES(d,r1,&sparse);
            start1 = clock();
            GMRES_m(d,r1,&sparse);
            end1 = clock();
            elapsed_GMRES = ((double) (end1 - start1)) / CLOCKS_PER_SEC;
            printf("\nXronos gia serial GMRES=%.5lfs\n",elapsed_GMRES);
        }
        for (i=0; i<=N; i++) {
            u[i]+=d[i];
        }
        iterat++;
    }
        //end while
        printf("##################################################\n");
        printf("##################################################\n");
        printf("\t\t Solving of 2-D non linear Bratu Problem\n");
        printf("\n Solver: GMRES Iterative\tMatrix format:CSR\t\n");
        printf("\n nex=%d ney=%d ne=%d N=%d\n",nex,ney,elem1.ne,N);
        printf("\n m(Krylov)=%d\tGMRES iterations=%d\tNewton Iterations=%d\n",m,GMRES_iter,Newton_iter);
        printf("\n nzeros=%ld, Nz for process zero=%d\n",sparse.nzeros,Nz);
        printf("\n Problem Solved, Newton Converged\n");
        printf("\n#############  SOLUTION  #############");
        printf("\n\n\n u=%.12lf\t lamda=%lf \t\n\n\n",u[(N)/2],param.lamda);
return 0;
}

