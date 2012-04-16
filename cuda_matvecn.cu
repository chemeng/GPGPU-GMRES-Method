//
//  cuda_matvecn.cu
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

//y = b - A*x upologizei

///////////////////////////////////////////////////////
//////       MONO GIA SPARSE MATRICES           ///////
//////       ME <32 NON-ZEROS PER ROW           ///////
///////////////////////////////////////////////////////
__global__ void cuda_matvecn_kernel(int dev_dim, double *y,double *AA, int *JA,
                                    int *IA, double *x, double *b,double *dot_res)

{
    __shared__ double cache[threadsPerBlock]; //thread shared memory
    int thread_id=threadIdx.x+blockIdx.x*blockDim.x; //global thread index
    int cacheIndex=threadIdx.x;
    double temp1 = 0;
    double temp2 = 0;
    int i = thread_id , jj=0;
    if (i<(dev_dim)) {
        while (i < (dev_dim)) {
            temp1 = 0;
            temp2 = 0;
            for( jj = IA[i] ; jj < IA[i+1]; jj ++ ){
                temp1 += AA[jj] * x[JA[jj]];
            }
            y[i] = b[i]-temp1;
            temp2 += y[i] * y[i];
            i += blockDim.x * gridDim.x;
        }
    }
    // set the cache values
    cache[cacheIndex] = temp2;
    // synchronize threads in this block
    __syncthreads();
    if (blockDim.x >= 512  && threadIdx.x < 256) {
        cache[threadIdx.x] += cache[threadIdx.x + 256];
        __syncthreads();
    }
    if (blockDim.x >= 256  && threadIdx.x < 128) {
        cache[threadIdx.x] += cache[threadIdx.x + 128];
        __syncthreads();
    }
    if (blockDim.x >= 128  && threadIdx.x < 64) {
        cache[threadIdx.x] += cache[threadIdx.x + 64];
        __syncthreads();
    }
    //unroll last warp no sync needed
    if (threadIdx.x <32 ) {
        if (blockDim.x >= 64) cache[threadIdx.x] += cache[threadIdx.x +32];
        if (blockDim.x >= 32) cache[threadIdx.x] += cache[threadIdx.x +16];
        if (blockDim.x >= 16) cache[threadIdx.x] += cache[threadIdx.x +8];
        if (blockDim.x >= 8) cache[threadIdx.x] += cache[threadIdx.x +4];
        if (blockDim.x >= 4) cache[threadIdx.x] += cache[threadIdx.x +2];
        if (blockDim.x >= 2) cache[threadIdx.x] += cache[threadIdx.x +1];
    }   
    if (cacheIndex==0) {
        dot_res[blockIdx.x]=cache[0];
    }
}

