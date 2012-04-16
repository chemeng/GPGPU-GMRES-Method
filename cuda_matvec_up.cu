//
//  cuda_matvec_up.cu
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


//Matvec+vector update-> x=A*y+x opou A non-sparse matrix
__global__ void cuda_matvec_up_kernel(int Ndim,int mdim, double *x,
                                      double *A, double *y )
{
    int i=0,k=0;
    int thread_id=threadIdx.x+blockIdx.x*blockDim.x; //global thread 
    while (thread_id < Ndim) {
        for(i = 0  ; i < Ndim*mdim ; i+=Ndim){
            x[thread_id]+= A[thread_id+i]*y[k] ; 
            k++;
        }
        thread_id += blockDim.x * gridDim.x ;
    }
} 
 
