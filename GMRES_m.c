//
//  gmres.c
//  Bratu-C
//
//  Created by Tim Ioannidis on 12/19/11.
//  Copyright 2011 Chemeng NTUA. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h> 
#include "structs.h"
#include "extern.h"
#include "gmres_config.h"
#include "parameters.h"
#include "methods_decl.h"


void GMRES_m(double *d,double *r1,struct common4 *sparse)
{
    int i=0,j=0,k=0,l=0,k1=0,k2=0,iter=0,tag1=0,tag2=0,right=0,left=0,ierr;
    double tot_rest=1,vita=0,eps=0;
    double e[m+1],y[m],g[m+1];
    double *x,*x0,*r0,*uj,*temp;
    //allocation wste na apothikeutoun sto heap
    x=(double*)calloc(N,sizeof(double));
    x0=(double*)calloc(N,sizeof(double));
    r0=(double*)malloc(N*sizeof(double));
    uj=(double*)malloc(N*sizeof(double));
    temp=(double*)malloc(N*sizeof(double));    
    iter=1;
    printf("Serial GMRES started computation\n");
    while (iter<=GMRES_iter) {
        for (i=0; i<N; i++) {
            x0[i]=x[i];
            uj[i]=0;
            r0[i]=0;
            temp[i]=0;
        }
            //midenismos pinakwn-anusmatwn Krylov
        for (i=0; i<N; i++) {
            for (j=0; j<m; j++) {
                gelim.u_base[i][j]=0;
            }
        }
        for (i=0; i<(m+1); i++) {
            e[i]=0;
            g[i]=0;
            for (j=0; j<m; j++) {
                gelim.Hm[i][j]=0;
                y[j]=0;
            }
        }
        e[0]=1;
        //upologismos r0=b-A*x0 opou b=r1 me matmul se CSR format parallel
        for (i=0; i<N; i++) {
            k1=sparse->IA[i];
            k2=sparse->IA[i+1]-1;
            for (j=k1; j<=k2; j++) {
                r0[i]+=sparse->AA[j]*x0[sparse->JA[j]];
            }
        }
        for (i=0; i<N; i++) {
            r0[i]=r1[i]-r0[i];
        }
        //end CSR multi
        vita=0;
            //upologismos normas tou r0 sto vita->prosoxi sta sunoriaka
        for (i=0; i<N; i++) {
                vita+=pow(r0[i],2);  
        }
        vita=sqrt(vita);
            //upologismos arxikou dianusmatos uj k apothikeusi ston u+base
        for (i=0;i<N; i++) {
            uj[i]=r0[i]/vita;
        }
            //apothikeusi tou u1 ston pinaka vasewn tou Krylov
        for (i=0; i<N; i++) {
            gelim.u_base[i][0]=uj[i];
        }
            //dianusma g=b*e
        for (i=0; i<(m+1); i++) {
            g[i]=vita*e[i];
        }
            //kalesma synartisewn
        Krylov(sparse,uj);
        Least_sq(g);
        Back_sub(y,g);
        //x=x0+MATMUL(u_base(N,m),y(m,1)))
        for (i=0; i<N; i++) {
            for (k=0; k<m; k++) {
                temp[i]+=gelim.u_base[i][k]*y[k];
            }
        }
        for (i=0; i<N; i++) {
            x[i]=x0[i]+temp[i];
        }
        iter++;
    }
    if ((linear_solver[0]=='s') || (linear_solver[0]=='S')) {
        for (i=0; i<N; i++) {
            d[i]=x[i];
        }
    }
    printf("SERIAL=%.15lf\n",x[N/2]);
    free(x);
    free(x0);
    free(r0);
    free(uj);
    free(temp);
}