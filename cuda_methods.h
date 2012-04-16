//
//  cudamethods.h
//  Cuda GMRES
//
//  Created by Tim Ioannidis on 2/18/12.
//  Copyright 2012 Chemeng NTUA. All rights reserved.
//

#ifndef Cuda_GMRES_cudamethods_h
#define Cuda_GMRES_cudamethods_h

//declaration of cuda kernels used
__global__ void cuda_dot_kernel();
__global__ void cuda_matvecn_kernel();
__global__ void cuda_matvec_kernel();
__global__ void cuda_vec_update_kernel();
__global__ void cuda_vec_replace_kernel();
__global__ void cuda_vec_replace2_kernel();
__global__ void cuda_initial2_kernel();
__global__ void cuda_matvec_up_kernel();
__global__ void cuda_main_kernel();
__global__ void cuda_leastsq_kernel();
__global__ void cuda_reduction_kernel();
__global__ void cuda_comm_kernel();


#endif
