//
//  cuda_reduction.cu
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

__global__ void cuda_reduction_kernel(int blocks,double *dev_res, double *red_res,int choice)
{
    __shared__ double cache[threadsPerBlock]; //thread shared memory
    int j=threadIdx.x,l=0;
    int p=blocks - blockDim.x ;
    if ( p>0 ) {
        cache[j]=dev_res[j];
        l=1;
        while ((j+l*blockDim.x) < blocks) {
            cache[j]+=dev_res[j+l*blockDim.x];
            l++;
        }
    }
    else {
        if (j<blocks) {
            cache[j]=dev_res[j];
        }
        else{
            cache[j]=0;
        }
    }
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
    if (j==0) {
        if (choice==1) {
            *red_res=sqrt(cache[0]);
        }
        else
        {
            *red_res=cache[0];        
        }
    }
}
