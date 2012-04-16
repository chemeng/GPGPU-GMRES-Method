//
//  cuda_comm.cu
//  Cuda GMRES
//
//  Created by Tim Ioannidis on 3/06/12.
//  Copyright 2012 Chemeng NTUA. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cuda_config.h"
#include "cuda_methods.h"


//dot product dot_res=a<dot>b me diastasi dim
__global__ void cuda_comm_kernel(double *dest,double *source,int choice)
{
    if (choice==0) {
        (*dest)=sqrt((*dest)+(*source));
    }
    else
    {    
        *dest += *source;
    }
}

 
