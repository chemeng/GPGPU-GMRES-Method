//
//  cuda_vec_update.cu
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


//vector update y=y+ax opou y,x vectors kai a arithmos
__global__ void cuda_vec_update_kernel(int n,double *y,double *a,double *x)
{
    int global_tid=0;
    //orismos indexing
    global_tid = threadIdx.x + blockIdx.x * blockDim.x; 
    while (global_tid < n) {
        y[global_tid] -= ((*a)*x[global_tid]);
        global_tid += blockDim.x * gridDim.x;
    }
}
 
