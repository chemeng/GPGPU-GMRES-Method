//
//  cuda_leastsq_kernel.cu
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

__global__ void cuda_leastsq_kernel(int m1,double *Hm,double *g,double *y,double *W)
{
    //Serial Least Squares problem
    //Hm kai W kata grammi flattened!!!!
    //shared variables
    __shared__ double si;
    __shared__ double ci;
    double temp[threadsPerBlock],help;
    int thread_id=threadIdx.x + blockIdx.x*blockDim.x;
    int j=thread_id,p=0,q=0,i=0;
    for (i=0; i<(m1); i++) {
        j=thread_id;
        //upologismos si kai ci
        if (j==0) {
            si=Hm[(i+1)*(m1) + i]/sqrt(pow(Hm[i*(m1) + i],2) + pow(Hm[(i+1)*(m1) + i],2));
            ci=Hm[i*(m1) + i]/sqrt(pow(Hm[i*(m1) + i],2) + pow(Hm[(i+1)*(m1) + i],2));
        }
        while (j < ((m1)+1)*((m1)+1)) {
            if ((j%((m1)+2))==0) {
                W[j]=1; //stoixeia tis kurias diagwniou
            }
            else
            {
                W[j]=0;
            }
            j+=blockDim.x*gridDim.x;
        }
        __syncthreads();
        j=thread_id;
        //eisagwgi timwn pinaka stin i,i+1 grammi
        if (j==0) {
            W[i*((m1)+1)+i]=ci;
            W[i*((m1)+1)+i+1]=si;
            W[(i+1)*((m1)+1)+i]=-si;
            W[(i+1)*((m1)+1)+i+1]=ci;
        }
        __syncthreads();
        //Hm=MATMUL(W(m+1,m+1),Hm(m+1,m)
        if (j<((m1)+1)) {
            for (q=0; q < (m1); q++) {
                temp[q]=0;
                for (p=0; p<((m1)+1); p++) {
                    temp[q]+=W[j*((m1)+1)+p]*Hm[q+p*(m1)];
                }
            }
            //g=MATVEC(W,g)
            help=0;
            for(p=0;p<((m1)+1);p++){
                help+=W[j*((m1)+1)+p]*g[p];
            }
        }
        __syncthreads();
        if (j<((m1)+1)) {
            for (q=0; q < (m1); q++) {
                Hm[j*(m1)+q]=temp[q];
            }
            g[j]=help;
        }
    }
    __syncthreads();
    //Serial Back_Sub
    if (j==0) {
        for (i=(m-1);i>=0;i--){
            help=g[i];
            if (i<(m-1)){
                for (p=(i+1);p<m;p++){
                    help-=(Hm[i*m + p])*y[p];
                }
            }
            y[i]=help/(Hm[i*m + i]);
        }
        
    }
} 
    
    
    
    
    
    

