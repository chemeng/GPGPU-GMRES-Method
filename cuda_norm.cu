//
//  cuda_norm.cu
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


//euclidian norm of vector vec
__global__ void cuda_norm_kernel(int n,double *a, double *dot_res)
{
    
    __shared__ double cache[threadsPerBlock]; //thread shared memory
    int thread_id=0,cacheIndex=0;
    double temp1 = 0; 
    //orismos indexing
    thread_id = threadIdx.x + blockIdx.x * blockDim.x; 
    cacheIndex = threadIdx.x;
    while (thread_id < n) {
        temp1 += a[thread_id] * a[thread_id];
        thread_id += blockDim.x * gridDim.x;
    }
    // set the cache values
    cache[cacheIndex] = temp1;
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
