//
//  cudaGMRES.cu
//  Cuda GMRES 
//
//  Created by Tim Ioannidis on 2/18/12.
//  Copyright 2012 Chemeng NTUA. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <cuda.h>
#include "structs.h"
#include "parameters.h"
#include "extern.h"
#include "cuda_config.h"
#include "cuda_methods.h"
#include "cuda_dot.cu"
#include "cuda_initial2.cu"
#include "cuda_matvec.cu"
#include "cuda_matvecn.cu"
#include "cuda_vec_update.cu"
#include "cuda_reduction.cu"
#include "cuda_norm.cu"
#include "cuda_vec_replace.cu"
#include "cuda_vec_replace2.cu"
#include "cuda_matvec_up.cu"
#include "cuda_leastsq.cu"


///////////////////////////////////////////////////////
//////          SOS!!!!! O DISDIASTATOS         ///////
//////        PINAKAS u_base EINAI KATA STILI   ///////
//////          E.G. A[i][j]=A[j*ROWS+i]        ///////
///////////////////////////////////////////////////////

extern "C"
{
    //lunei to provlima A*x=r1 kai gemizei to d=x
    void cuda_GMRES( double *d, double *r1, struct common4 *sparse)
    {   
        clock_t start1, end1 ;
        float cuda_GMRES_time=0;
        int blocksPerGrid;
        blocksPerGrid=((N+threadsPerBlock-1)/threadsPerBlock);
        if (blocksPerGrid > 65530) {
            printf("WARNING,block number exceeded hardware limit");
            blocksPerGrid=65530;
        }
        //streams and devices
        cudaStream_t stream0;
        cudaSetDevice(0);
        cudaStreamCreate( &stream0 );
        //de metraw ta apo panw sto xrono giati einai to initialization tis kartas
        printf("ThreadsPerBlock=%d\n",threadsPerBlock);
        printf("\nCuda GMRES started computation\n");
        start1 = clock(); 
        //variables declaration
        int iter=0,i=0,j=0;
        double *dev0_AA,*dev0_r1,*dev0_help;
        int *dev0_JA,*dev0_IA;
        double *dev0_x,*dev0_r0,*dev0_w,*dev0_res,*dev0_vita,*dev0_Wm;
        double *dev0_Hm,*dev0_u_base,*dev0_e,*dev0_y,*dev0_g;
    //allocation sto device arrays me dedomena
        cudaMalloc((void**)&dev0_r1,(N)*sizeof(double));
        cudaMalloc((void**)&dev0_AA,(Nz)*sizeof(double));
        cudaMalloc((void**)&dev0_JA,(Nz)*sizeof(int));
        cudaMalloc((void**)&dev0_IA,(N+1)*sizeof(int));
    //allocation sto device voithikwn arrays
        cudaMalloc((void**)&dev0_x,(N)*sizeof(double));
        cudaMalloc((void**)&dev0_r0,(N)*sizeof(double));
        cudaMalloc((void**)&dev0_Hm,((m+1)*m)*sizeof(double));
        cudaMalloc((void**)&dev0_u_base,(N*m)*sizeof(double));
        cudaMalloc((void**)&dev0_e,(m+1)*sizeof(double));
        cudaMalloc((void**)&dev0_y,(m)*sizeof(double));
        cudaMalloc((void**)&dev0_g,(m+1)*sizeof(double));
        cudaMalloc((void**)&dev0_w,(N)*sizeof(double));
        cudaMalloc((void**)&dev0_res,(blocksPerGrid)*sizeof(double));
        cudaMalloc((void**)&dev0_vita,sizeof(double));
        cudaMalloc((void**)&dev0_Wm,((m+1)*(m+1))*sizeof(double));
        cudaMalloc((void**)&dev0_help,sizeof(double));
    //perasma dedomenwn stin global device memory
        cudaMemcpy(dev0_AA, sparse->AA, Nz*sizeof(double), cudaMemcpyDefault );
        cudaMemcpy(dev0_JA, sparse->JA, Nz*sizeof(int), cudaMemcpyDefault );
        cudaMemcpy(dev0_IA, sparse->IA, (N+1)*sizeof(int), cudaMemcpyDefault );
        cudaMemcpy(dev0_r1, r1, N*sizeof(double), cudaMemcpyDefault );

/////////////////////////////////////////////////////////////////////////////////                
        //ksekinima epanaliptikis
        iter=1;
        cuda_initial2_kernel<<<blocksPerGrid,threadsPerBlock,0,stream0>>>((m+1),dev0_e,N,
                                                                          dev0_x);
        while (iter<=GMRES_iter) {
            //upologismos r0=b-A*x opou b=r1 me MATMUL se CSR format kai NORM r0          
            cuda_matvecn_kernel<<<blocksPerGrid,threadsPerBlock,0,stream0>>>(N, dev0_r0, dev0_AA, dev0_JA,
                                                                       dev0_IA, dev0_x, dev0_r1,dev0_res);
            cuda_reduction_kernel<<<1,threadsPerBlock,0,stream0>>>(blocksPerGrid, dev0_res,dev0_vita,1);
            //upologismos uj[]=r0[]/vita kai apothikeusi ston u_base KATA STILI orismenos
            //tautoxrona g=vita*e
            cuda_vec_replace_kernel<<<blocksPerGrid,threadsPerBlock,0,stream0>>>(  N,
                        dev0_u_base, dev0_vita, dev0_r0, m+1, dev0_g, dev0_e);    
            //KATASKEUI UPOXWROU Krylov
            for (j=0; j<m; j++) {   //j einai to count g mas
                if (j >= 1) {
                    //u_base[][j+1]=w[]/Hm[j+1][j]
                    cuda_vec_replace2_kernel<<<blocksPerGrid,threadsPerBlock,0,stream0>>>(N,
                                            &dev0_u_base[(j)*N], &dev0_Hm[(j)*m + j-1], dev0_w);
                }
                //uj[i]=u_base[i][j] 
                //matmul me CSR w=matvec(A,uj)
                cuda_matvec_kernel<<<blocksPerGrid,threadsPerBlock,0,stream0>>>( N,
                                    dev0_w, dev0_AA, dev0_JA, dev0_IA,&dev0_u_base[j*N]);  
                for (i=0; i<=j; i++) {
                    //uj[k]=u_base[k][i]
                    //DOT PRODUCT w*uj-> kai eisagwgi sto Hm[i][j]
                    cuda_dot_kernel<<<blocksPerGrid,threadsPerBlock,0,stream0>>>( N,dev0_w,
                                                                       &dev0_u_base[i*N],dev0_res);
                    cuda_reduction_kernel<<<1,threadsPerBlock,0,stream0>>>(blocksPerGrid, dev0_res,&dev0_Hm[i*m + j],2);   
                    //w=w-Hm(i,j)*uj
                    cuda_vec_update_kernel<<<blocksPerGrid,threadsPerBlock,0,stream0>>>( N,
                                                        dev0_w, &dev0_Hm[i*m + j], &dev0_u_base[i*N] );       
                }
                cuda_norm_kernel<<<blocksPerGrid,threadsPerBlock,0,stream0>>>(N,
                                                    dev0_w,dev0_res);
                cuda_reduction_kernel<<<1,threadsPerBlock,0,stream0>>>(blocksPerGrid, dev0_res,&dev0_Hm[(j+1)*m + j],1);
               /* if (j<(m-1)) {
                //u_base[][j+1]=w[]/Hm[j+1][j]
                    cuda_vec_replace2_kernel<<<blocksPerGrid,threadsPerBlock,0,stream0>>>(N,
                                        &dev0_u_base[(j+1)*N], &dev0_Hm[(j+1)*m + j], dev0_w);
                }*/
            }
            //Least Squares problem
            if (threadsPerBlock<m) {
                printf("ERROR, threadsPerBlock should be greater than m");
            }
            cuda_leastsq_kernel<<<1,threadsPerBlock,0,stream0>>>(m,dev0_Hm,dev0_g,dev0_y,dev0_Wm);
            //TELOS Least Squares problem
            //upologismos x = x0 + matvec(u_base(N,m),y(m))
            cuda_matvec_up_kernel<<<blocksPerGrid,threadsPerBlock,0,stream0>>>(N, m, dev0_x, dev0_u_base, dev0_y);
            iter++;
        }
        //Copy result back to CPU
        cudaMemcpy(d, dev0_x, N*sizeof(double), cudaMemcpyDefault );
        //Free memory
        printf("CUDA=%.15lf\n",d[N/2]);
        cudaFree(dev0_AA);
        cudaFree(dev0_JA);
        cudaFree(dev0_IA);
        cudaFree(dev0_r1);
        cudaFree(dev0_x);
        cudaFree(dev0_r0);
        cudaFree(dev0_u_base);
        cudaFree(dev0_Hm);
        cudaFree(dev0_y);
        cudaFree(dev0_g);
        cudaFree(dev0_e);
        cudaFree(dev0_w);
        cudaFree(dev0_Wm);
        cudaFree(dev0_vita);
        cudaFree(dev0_res);
        end1 = clock();      
        cuda_GMRES_time = ((double) (end1 - start1)) / CLOCKS_PER_SEC;
        printf("\nXronos gia Cuda_GMRES=%.5lfs\n\n",cuda_GMRES_time);
    }
}


