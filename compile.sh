#!bin/sh

icc -c main.c
nvcc -c cuda_GMRES.cu -arch=compute_20

icc main.o cuda_GMRES.o -L/share/apps/cuda/lib64 -lcudart
