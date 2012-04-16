//
//  Least_sq.c
//  Bratu-C
//
//  Created by Tim Ioannidis on 12/20/11.
//  Copyright 2011 Chemeng NTUA. All rights reserved.
//

#include <stdio.h>
#include <math.h> 
#include "structs.h"
#include "extern.h"
#include "parameters.h"
#include "methods_decl.h"

void Least_sq(double *g)
{
    int i=0,j=0,k=0,p=0,q=0;
    double W[m+1][m+1],si=0,ci=0;
    double temp1[m+1][m],temp2[m+1];
    
    for (i=0;i<m;i++){
            //upologismos si kai ci
        si=(gelim.Hm[i+1][i])/sqrt(pow(gelim.Hm[i][i],2)+pow(gelim.Hm[i+1][i],2));
        ci=(gelim.Hm[i][i])/sqrt(pow(gelim.Hm[i][i],2)+pow(gelim.Hm[i+1][i],2));
        for (j=0;j<(m+1);j++){
            for (p=0; p<m; p++) {
                temp1[j][p]=0;
            }
            for (k=0;k<(m+1);k++){
                if(j==k){
                    W[j][k]=1;
                }
                else {
                    W[j][k]=0;
                }
            }
            temp2[j]=0;
        }
            //eisagwgi timwn pinaka stin i,i+1 grammi
        W[i][i]=ci;
        W[i][i+1]=si;
        W[i+1][i]=-si;
        W[i+1][i+1]=ci;
            //Hm=MATMUL(W(m+1,m+1),Hm(m+1,m)
        for (j=0;j<(m+1);j++){
            for(k=0;k<m;k++){
                for(p=0;p<(m+1);p++){
                    temp1[j][k]+=(W[j][p]*(gelim.Hm[p][k]));
                }
            }
        }
        for (j=0; j<(m+1); j++) {
            for (k=0; k<m; k++) {
                gelim.Hm[j][k]=temp1[j][k];
            }
        }
        for (j=0;j<(m+1);j++){
            for(p=0;p<(m+1);p++){
                temp2[j]+=W[j][p]*g[p];
            }
        }
        for (j=0; j<(m+1); j++) {
            g[j]=temp2[j];
        }
    }
}
