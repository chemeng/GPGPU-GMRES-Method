//
//  xycoord.c
//  Bratu-C
//
//  Created by Tim Ioannidis on 12/24/11.
//  Copyright 2011 Chemeng NTUA. All rights reserved.
//

#include <stdio.h>
#include <math.h> 
#include "extern.h"
#include "structs.h"
#include "methods_decl.h"
#include "parameters.h"

void xycoord(struct common1 *mesh1,struct common2 *elem1,float *xpt,float *ypt)
{
	int nnode=0,i=0,j=0;
    // (x,y) coordinates of each node
    xpt[0]=(*mesh1).xorigin;
    ypt[0]=(*mesh1).yorigin;
    for (i=1;i<=(*elem1).nnx;i++){	
        nnode=(i-1)*((*elem1).nny);
        xpt[nnode]=(xpt[0])+(i-1)*(((*mesh1).deltax)/2);
        ypt[nnode]=(ypt[0]);
        for (j=2;j<=((*elem1).nny);j++){
            xpt[nnode+j-1]=xpt[nnode];
            ypt[nnode+j-1]=ypt[nnode]+(j-1)*(((*mesh1).deltay)/2);
        }
    }
}