//
//  mesh.c
//  Bratu-C
//
//  Created by Tim Ioannidis on 12/19/11.
//  Copyright 2011 Chemeng NTUA. All rights reserved.
//

#include <stdio.h>
#include <math.h> 
#include "extern.h"
#include "structs.h"
#include "parameters.h"
#include "methods_decl.h"

void xydiscr(struct common1 *mesh1)
{
	mesh1->xorigin=0.0;
	mesh1->yorigin=0.0;
	mesh1->xlast=1.0;
	mesh1->ylast=1.0;
	mesh1->deltax=(float)((mesh1->xlast-(mesh1->xorigin))/nex);
	mesh1->deltay=(float)((mesh1->ylast-(mesh1->yorigin))/ney);
}






