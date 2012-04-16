//
//  cuda_config.h
//  Cuda GMRES
//
//  Created by Tim Ioannidis on 2/18/12.
//  Copyright 2012 Chemeng NTUA. All rights reserved.
//

#ifndef Cuda_GMRES_cuda_config_h
#define Cuda_GMRES_cuda_config_h

//orismos stoixeiwn gia CUDA gmres
//threadsPerBlock POLU simantiko sunithws 128 i 256, na to checkarw
extern const int threadsPerBlock=512;

//ta blocks pou xrisimopoiountai kathorizontai apo to calling function->cuda_GMRES
//ws blocksPerGrid=(N+threadsPerBlock-1)/threadsPerBlock

#endif
