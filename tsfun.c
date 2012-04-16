//
//  tsfun.c
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

void tsfun(double x,double y,double *phi,double *phic,double *phie)
{
	double l1x,l2x,l3x,l1y,l2y,l3y,dl1x,dl2x,dl3x,dl1y,dl2y,dl3y;
    int i,j;
	l1x=2*pow(x,2)-3*x+1;
	l1y=2*pow(y,2)-3*y+1;
	l2x=-4*pow(x,2)+4*x;
	l2y=-4*pow(y,2)+4*y;	
	l3x=2*pow(x,2)-x;  
	l3y=2*pow(y,2)-y;
	dl1x=4*x-3;
	dl1y=4*y-3;
	dl2x=-8*x+4;
	dl2y=-8*y+4;
	dl3x=4*x-1;
	dl3y=4*y-1;
	phi[0]=l1x*l1y;
	phi[1]=l1x*l2y;
	phi[2]=l1x*l3y;
	phi[3]=l2x*l1y;
	phi[4]=l2x*l2y;
	phi[5]=l2x*l3y;
	phi[6]=l3x*l1y;
	phi[7]=l3x*l2y;
	phi[8]=l3x*l3y;
	phic[0]=l1y*dl1x;
	phic[1]=l2y*dl1x;
	phic[2]=l3y*dl1x;
	phic[3]=l1y*dl2x;
	phic[4]=l2y*dl2x;
	phic[5]=l3y*dl2x;
	phic[6]=l1y*dl3x;
	phic[7]=l2y*dl3x;
	phic[8]=l3y*dl3x;
	phie[0]=l1x*dl1y;
	phie[1]=l1x*dl2y;
	phie[2]=l1x*dl3y;
	phie[3]=l2x*dl1y;
	phie[4]=l2x*dl2y;
	phie[5]=l2x*dl3y;
	phie[6]=l3x*dl1y;
	phie[7]=l3x*dl2y;
	phie[8]=l3x*dl3y;
}  

