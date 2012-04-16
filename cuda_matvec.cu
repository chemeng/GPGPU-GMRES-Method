//
//  cuda_matvec.cu
//  Cuda GMRES
//
//  Created by Tim Ioannidis on 2/18/12.
//  Copyright 2012 Chemeng NTUA. All rights reserved.
//
//CSR y=A*x mutliplication using CUDA
//ptr->IA       indices->JA     data->AA

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cuda_config.h"
#include "cuda_methods.h"

//PREPEI NA UPARXOUN STIN GLOBAL TIS DEVICE AA,JA,IA,x,y

//y = A*x upologizei

///////////////////////////////////////////////////////
//////       MONO GIA SPARSE MATRICES           ///////
//////       ME <32 NON-ZEROS PER ROW           ///////
///////////////////////////////////////////////////////
__global__ void cuda_matvec_kernel(int dev_dim, double *y, double *AA, int *JA, int *IA, 
                                   double *x )
{
    int thread_id=threadIdx.x+blockIdx.x*blockDim.x; //global thread index
    int i = thread_id , jj=0;
    if (i<dev_dim) {
        while (i < dev_dim) {
            y[i] = 0;
            for( jj = IA[i] ; jj < IA[i+1]; jj ++ ){
               y[i] += AA[jj] * x[JA[jj]];
            }
            i += blockDim.x * gridDim.x;
        }
    }
}
