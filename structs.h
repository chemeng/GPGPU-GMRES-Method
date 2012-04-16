//
//  structs.h
//  Bratu-C
//
//  Created by Tim Ioannidis on 12/24/11.
//  Copyright 2011 Chemeng NTUA. All rights reserved.
//

#ifndef Bratu_C_structs_h
#define Bratu_C_structs_h

#include "parameters.h"

    //structs declaration
struct common1
{
    float xorigin,yorigin,xlast,ylast,deltax,deltay;
};
struct common2
{
    int nnx,nny,ne,np,nop[nex*ney][9];
};
struct common3
{
    double Hm[m+1][m];
    //u_base[N][m]
    double u_base[(2*nex+1)*(2*ney+1)][m];
};
struct common4
{
    double *AA;
    int *JA,*IA;
    long int nzeros;
};
struct common5
{
    double alfa,lamda;
};

#endif
