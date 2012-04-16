//
//  cuda_vec_replace.cu
//  Cuda GMRES
//
//  Created by Tim Ioannidis on 2/18/12.
//  Copyright 2012 Chemeng NTUA. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cuda_config.h"
#include "cuda_methods.h"


//vector replace y1=a*x1 opou y,x vectors kai a arithmos kai y2=x2/a
__global__ void cuda_vec_replace_kernel(int n1,double *y1 , double *a, double *x1, int n2 ,double *y2, double *x2)
{
    int global_tid=0;
    //orismos indexing
    global_tid = threadIdx.x + blockIdx.x * blockDim.x; 
    while (global_tid < n1) {
        y1[global_tid] = (1/(*a))*x1[global_tid];
        global_tid += blockDim.x * gridDim.x;
    }
    global_tid = threadIdx.x + blockIdx.x * blockDim.x; 
    while (global_tid < n2) {
        y2[global_tid] = (*a)*x2[global_tid];
        global_tid += blockDim.x * gridDim.x;
    }

}
 
 
