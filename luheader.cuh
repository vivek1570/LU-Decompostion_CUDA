#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <time.h>

#define BLOCK_SIZE 512

void printarray(double *mat,int width,int height){
	int i,j;
	for (i=0;i<height;i++){
		for (j=0;j<width;j++){
			printf("%f ",mat[width*i+j]);
			}
		printf("\n");
	}
}

int hfact(int size,int maxsize){ 
  int i=1;int factor=1;
  for (i=2;i<=maxsize;i++){
	if(size%i==0){
		factor=i;
	}
  }
  return factor;
}

__global__ void initialize(double *matL,double *matU,int size){
	int i,j;
	
	i= blockIdx.y * blockDim.y + threadIdx.y;
	j= blockIdx.x * blockDim.x + threadIdx.x;

	if(i<size && j<size){
		matU[i*size+j]=0;
		matL[i*size+j]=0;
		if (i==j){
			matL[i*size+j]=1;
		}
	}


}


//do LU decomposition in parallel nth row nth col using 1 dimensional thread array	
__global__ void luop(double *matL, double *matU, double *matA, int size, int curr, float *upperTime, float *lowerTime) {
	int i,j,k,n;
	double u,l;
	n = blockIdx.x * blockDim.x + threadIdx.x;

	if (n >= curr && n < size) {
			i=curr; 
			j=n;
			u =matA[i*size+j];
			if (j >= i) {
					float startU =clock();  
					for (k= 0; k < i; k++) {
							u= u-matU[k*size+j] * matL[i*size+k]; 
					}
					matU[i * size + j] = u;
					float endU = clock();
					atomicAdd(upperTime,(endU-startU)/CLOCKS_PER_SEC); 
			}

			i=n;
			j=curr;
			l= matA[i*size+j];
			if (j <= i) {
					float startL =clock();
					for (k = 0; k < j; k++) {
							l=l-matU[k*size+j] * matL[i*size+k];
					}
					matL[i*size+j] = l/(double)matU[j*size+j];
					float endL=clock();
					atomicAdd(lowerTime,(endL-startL)/CLOCKS_PER_SEC);
			}
	}
}



void ludecompose(double *matL,double *matU,double *matA,int size,FILE *fp1){
	double *dev_L,*dev_U,*dev_A;

  cudaMalloc((void**)&dev_L, size*size*sizeof(double));
  cudaMalloc((void**)&dev_U, size*size*sizeof(double));
  cudaMalloc((void**)&dev_A, size*size*sizeof(double));
	
  cudaMemcpy(dev_A,matA, size*size*sizeof(double),cudaMemcpyHostToDevice);
	
	//for initialization
	dim3 grid(ceil(size/(double)BLOCK_SIZE),ceil(size/(double)BLOCK_SIZE));
	dim3 block(BLOCK_SIZE,BLOCK_SIZE);

	//for calculation
	int block1=hfact(size,128);
	int grid1=ceil(size/block1);

	//Time
	clock_t x=clock();

	initialize<<<grid,block>>>(dev_L,dev_U, size);
	int i=0;



// Allocate memory for timing variables on the device
float *d_ut, *d_lt;
float h_ut = 0.0f, h_lt = 0.0f;

cudaMalloc(&d_ut, sizeof(float));
cudaMalloc(&d_lt, sizeof(float));
cudaMemcpy(d_ut, &h_ut, sizeof(float), cudaMemcpyHostToDevice);
cudaMemcpy(d_lt, &h_lt, sizeof(float), cudaMemcpyHostToDevice);

for (i=0;i<size; i++) {
    luop<<<grid1,block1>>>(dev_L, dev_U, dev_A, size, i, d_ut, d_lt);
    cudaDeviceSynchronize();
}

cudaMemcpy(&h_ut, d_ut, sizeof(float), cudaMemcpyDeviceToHost);
cudaMemcpy(&h_lt, d_lt, sizeof(float), cudaMemcpyDeviceToHost);

fprintf(fp1,"\nHere both h_ut and h_lt have getting some overheads due to the synchronization of atmoicadd\n it is clear that lower taking more time than upper\n");
fprintf(fp1,"Total time for upper triangular matrix calculations: %.10f s\n", h_ut);
fprintf(fp1,"Total time for lower triangular matrix calculations: %.10f s\n\n", h_lt);



cudaFree(d_ut);
cudaFree(d_lt);

	clock_t y=clock();
	double ct=(double)((y-x)/(double)CLOCKS_PER_SEC);
  fprintf(fp1,"Time taken for LU decomposition in total : %1.10f s\n",ct);
	
  cudaMemcpy(matL,dev_L,size*size*sizeof(double),cudaMemcpyDeviceToHost);
  cudaMemcpy(matU,dev_U,size*size*sizeof(double),cudaMemcpyDeviceToHost);
	cudaFree(dev_A);
	cudaFree(dev_L);
	cudaFree(dev_U);

}
