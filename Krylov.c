//
//  Krylov.c
//  Bratu-C
//
//  Created by Tim Ioannidis on 12/20/11.
//  Copyright 2011 Chemeng NTUA. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h> 
#include "structs.h"
#include "parameters.h"
#include "extern.h"
#include "methods_decl.h"

void Krylov(struct common4 *sparse,double *uj)
{
    int  i=0,j=0,k=0,k1=0,k2=0,p=0,tag1=0,tag2=0,count=0,l=0;   
    double *w;
    
    w=(double*)malloc(N*sizeof(double));
    for (j=0;j<m;j++){
            //j einai to count g mas
            //vazei to kommati tou uj pou antistoixei sto process
        for(k=0;k<N;k++){
                uj[k]=gelim.u_base[k][j];
        }
            //matmul me CSR format w=MATMUL(A,uj)
        for(i=0;i<N;i++){
            w[i]=0;
            k1=sparse->IA[i];
            k2=sparse->IA[i+1];
            for (k=k1;k<k2;k++){
                w[i]+=(sparse->AA[k])*(uj[sparse->JA[k]]);
            }
        }
        for (i=0;i<=j;i++){
            for(k=0;k<N;k++){
                    uj[k]=gelim.u_base[k][i];
            }
                //DOT_PRODUCT(w,uj)
            gelim.Hm[i][j]=0;
            for (k=0;k<N;k++){
                gelim.Hm[i][j]+=w[k]*uj[k];
            }
                //end DOT product
            //w=w-gelim.Hm*uj
            for(k=0;k<N;k++){
                w[k]-=(gelim.Hm[i][j])*(uj[k]);
            }
        }
        for (k=0;k<N;k++){
                gelim.Hm[j+1][j]+=pow(w[k],2);
        }
        gelim.Hm[j+1][j]=sqrt(gelim.Hm[j+1][j]);
        if(j<(m-1)){
            for (k=0;k<N;k++){
                    gelim.u_base[k][j+1]=w[k]/(gelim.Hm[j+1][j]);
            }
        } 
    }
}   

