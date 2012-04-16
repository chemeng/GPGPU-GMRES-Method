//
//  nodnumb.c
//  Bratu-C
//
//  Created by Tim Ioannidis on 12/24/11.
//  Copyright 2011 Chemeng NTUA. All rights reserved.
//

#include <stdio.h>
#include <math.h> 
#include "structs.h"
#include "extern.h"
#include "parameters.h"
#include "methods_decl.h"

void nodnumb(struct common2 *elem1)
{ 
	int i=0,j=0,k=0,l=0,nel=0;
    elem1->ne=nex*ney;
    elem1->nnx=2*nex+1;
    elem1->nny=2*ney+1;
    elem1->np=elem1->nnx*(elem1->nny);
    // nodal numbering
    nel=-1;
    for (i=1;i<=nex;i++){
        for (j=1;j<=ney;j++){
            nel++;
            for (k=1;k<=3;k++){
                l=3*k-3;
                elem1->nop[nel][l]=(elem1->nny)*(2*i+k-3)+2*j-2;
                elem1->nop[nel][l+1]=elem1->nop[nel][l]+1;
                elem1->nop[nel][l+2]=elem1->nop[nel][l]+2;
            }
        }
    }
}