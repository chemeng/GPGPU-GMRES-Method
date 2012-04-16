//
//  abfind.c
//  Bratu-C
//
//  Created by Tim Ioannidis on 12/19/11.
//  Copyright 2011 Chemeng NTUA. All rights reserved.
//


#include <stdio.h>
#include <math.h> 
#include "structs.h"
#include "extern.h"
#include "parameters.h"
#include "methods_decl.h"

void abfind(int nell,double *r1,double *u,struct common2 *elem1,float *xpt,float *ypt,struct common4 *sparse,struct common5 *param,int *ncod)
{
    double phi[9],phic[9],phie[9],usol=0,usol1=0,usol2=0,x1=0,x2=0,y1=0,y2=0,dett=0;
    int ngl[9];
    double tphx[9],tphy[9],value=0;
    double estifm[9][9];
    int i=0,j=0,k=0,p=0,q=0,t=0,r=0,ind=0,l=0,c1=0;
    double w[]={0.27777777777778,0.444444444444,0.27777777777778};
    double gp[]={0.1127016654,0.5,0.8872983346};
    double x_gp=0,y_gp=0;
    for (i=0; i<9; i++) {
        for (j=0; j<9; j++) {
            estifm[i][j]=0; 
        }
    }

    for (i=0; i<9; i++) {
        ngl[i]=elem1->nop[nell][i];  
    }
    for (j=0; j<3; j++) {
        for (k=0; k<3; k++) {
            x_gp=gp[j];
            y_gp=gp[k];
            tsfun(x_gp,y_gp,phi,phic,phie);
                //isoparametric transformation
            x1=0;
            x2=0;
            y1=0;
            y2=0;
            usol=0;
            usol1=0;
            usol2=0;
            for (c1=0; c1<9; c1++) {
                x1+=xpt[ngl[c1]]*phic[c1];
                x2+=xpt[ngl[c1]]*phie[c1];
                y1+=ypt[ngl[c1]]*phic[c1];
                y2+=ypt[ngl[c1]]*phie[c1];
            }
            dett=x1*y2-x2*y1;
            for (ind=0; ind<9; ind++) {
                tphx[ind]=(y2*phic[ind]-y1*phie[ind])/dett;
                tphy[ind]=(x1*phie[ind]-x2*phic[ind])/dett;
            }
            for (c1=0; c1<9; c1++) {
                    //athroisma sto oloklirwma
                usol+=u[ngl[c1]]*phi[c1];
                usol1+=u[ngl[c1]]*tphx[c1];
                usol2+=u[ngl[c1]]*tphy[c1];
            }
            
                //residuals  
            for (l=0; l<9; l++) {
                r1[ngl[l]]+=w[j]*w[k]*dett*(usol1*tphx[l]+usol2*tphy[l]);
                r1[ngl[l]]-=w[j]*w[k]*dett*(param->lamda)*phi[l]*exp(usol/(1+(param->alfa)*usol));
                for (c1=0; c1<9; c1++) {
                    estifm[l][c1]-=w[j]*w[k]*dett*((tphx[l]*tphx[c1])+(tphy[l]*tphy[c1]));
                    estifm[l][c1]+=w[j]*w[k]*dett*(param->lamda)*phi[l]*phi[c1]*exp(usol/(1+((param->alfa)*usol)))*(1/(pow((1+((param->alfa)*usol)),2)));
                }
            }
        }
    }
    for (t=0; t<9; t++) {
        for (r=0; r<9; r++) {
            value=estifm[t][r];
            p=ngl[t];
            q=ngl[r];
            CSR(nell,elem1,sparse,value,p,q,ncod);
        }
    }
    
}


