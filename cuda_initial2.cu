//
//  cuda_initial2.cu
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

//kernel initialization-> cuda_initial
//sets y=value

__global__ void cuda_initial2_kernel(int n1,double *y,int n2, double *x )
{
    int global_tid=0;
    //orismos indexing
    global_tid = threadIdx.x + blockIdx.x * blockDim.x; 
    while (global_tid < n1) {
        y[global_tid] = 0 ;
        if (global_tid==0) {
            y[global_tid]=1;
        }
        global_tid += blockDim.x * gridDim.x;
    }
    
    global_tid = threadIdx.x + blockIdx.x * blockDim.x; 
    while (global_tid < n2) {
        x[global_tid] = 0 ;
        global_tid += blockDim.x * gridDim.x;
    }
}
 



