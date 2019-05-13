/**********************************************************************/
// A now optimized Multigrid Solver for the Heat Equation             //
// Course Material for HPCSE-II, Spring 2019, ETH Zurich              //
// Authors: Sergio Martin, Georgios Arampatzis                        //
// License: Use if you like, but give us credit.                      //
/**********************************************************************/

#include <stdio.h>
#include <math.h>
#include <limits>
#include "heat2d_gpu.hpp"
#include "string.h"
#include <chrono>
#include <iostream>

pointsInfo __p;

void copyToDevice(gridLevel * g, size_t s,size_t gridCount){

  
  for (size_t l = s; l < gridCount; l++){
    cudaMemcpy(g[l].gU, g[l].U, sizeof(double)*g[l].N*g[l].N, cudaMemcpyHostToDevice); checkCUDAError("Error copying U");
    cudaMemcpy(g[l].gRes, g[l].Res, sizeof(double)*g[l].N*g[l].N, cudaMemcpyHostToDevice); checkCUDAError("Error copying Res");
    cudaMemcpy(g[l].gUn, g[l].Un, sizeof(double)*g[l].N*g[l].N, cudaMemcpyHostToDevice); checkCUDAError("Error copying Un");
    cudaMemcpy(g[l].gf, g[l].f, sizeof(double)*g[l].N*g[l].N, cudaMemcpyHostToDevice); checkCUDAError("Error copying f"); 
    cudaDeviceSynchronize();

  }
}


void copyToHost(gridLevel * g,size_t s , size_t gridCount){
  
  for (size_t l = s; l < gridCount; l++){
    cudaMemcpy(g[l].U, g[l].gU, sizeof(double)*g[l].N*g[l].N, cudaMemcpyDeviceToHost); checkCUDAError("Error copying U back");
    cudaMemcpy(g[l].Res, g[l].gRes, sizeof(double)*g[l].N*g[l].N, cudaMemcpyDeviceToHost); checkCUDAError("Error copying Res back");
    cudaMemcpy(g[l].Un, g[l].gUn, sizeof(double)*g[l].N*g[l].N, cudaMemcpyDeviceToHost); checkCUDAError("Error copying Un back");
    cudaMemcpy(g[l].f, g[l].gf, sizeof(double)*g[l].N*g[l].N, cudaMemcpyDeviceToHost); checkCUDAError("Error copying f back"); 
    cudaDeviceSynchronize();

  }
}

int main(int argc, char* argv[])
{
  double tolerance = 1e-0; // L2 Difference Tolerance before reaching convergence.
  size_t N0 = 7; // 2^N0 + 1 elements per side

  // Multigrid Parameters -- Find the best configuration!
  size_t gridCount       = N0-1;     // Number of Multigrid levels to use
  size_t downRelaxations = 5; // Number of Relaxations before restriction
  size_t upRelaxations   = 0;   // Number of Relaxations after prolongation

  gridLevel* g = generateInitialConditions(N0, gridCount);

  
  cudaDeviceProp deviceProp;
  cudaGetDeviceProperties(&deviceProp, 0);
  std::cout << "Using up to " << deviceProp.maxThreadsPerBlock << " threads per block." << std::endl;
  
  auto startTime = std::chrono::system_clock::now();
  
  
  //copyToDevice(g, gridCount);
  
  while (g[0].L2NormDiff > tolerance)  // Multigrid solver start
    {
  
      copyToDevice(g, 0,1);
      applyJacobi(g, 0, downRelaxations); // Relaxing the finest grid first
      calculateResidual(g, 0); // Calculating Initial Residual
    copyToHost(g, 0,1);
  

      for (size_t grid = 1; grid < gridCount; grid++) // Going down the V-Cycle
	{
	  //copyToHost(g, gridCount);
	  applyRestriction(g, grid); // Restricting the residual to the coarser grid's solution vector (f)
	  copyToDevice(g,grid, grid+1);
	  applyJacobi(g, grid, downRelaxations); // Smoothing coarser level
	  calculateResidual(g, grid); // Calculating Coarse Grid Residual
	  copyToHost(g, grid,grid+1);
  	  //  copyToDevice(g, gridCount);
	}

      for (size_t grid = gridCount-1; grid > 0; grid--) // Going up the V-Cycle
	{
	  //copyToHost(g, gridCount);
	  applyProlongation(g, grid); // Prolonging solution for coarser level up to finer level
	  copyToDevice(g, grid, grid+1);
	  applyJacobi(g, grid, upRelaxations); // Smoothing finer level
	  copyToHost(g, grid, grid+1);
	  
	}
      //copyToHost(g, gridCount);
      copyToDevice(g, 0, 1);
	  
      calculateL2Norm(g, 0); // Calculating Residual L2 Norm
    }  // Multigrid solver end

  //copyToHost(g, gridCount);
	
  auto endTime = std::chrono::system_clock::now();
  totalTime = std::chrono::duration<double>(endTime-startTime).count();
  printTimings(gridCount);
  printf("L2Norm: %.4f\n",  g[0].L2Norm);
  freeGrids(g, gridCount);
  return 0;
}

bool failed = 0;

void checkCUDAError(const char *msg)
{
  cudaError_t err = cudaGetLastError();
  if( cudaSuccess != err && !failed)
    {
      fprintf(stderr, "CUDA Error: %s: %s.\n", msg, cudaGetErrorString(err) );
      failed = 1;
      //exit(EXIT_FAILURE);
    }
}

__global__
void jacobi(double h2, double* U, double* Un, double* f, size_t N)
{
  size_t j = blockIdx.x*blockDim.x+threadIdx.x;
  if (j % N == 0 || j % N == N - 1 || j / N == 0 || j / N >= N - 1) return;
  U[j] = (Un[j-N] + Un[j+N] + Un[j-1] + Un[j+1] + f[j]*h2)*0.25;
  }

void applyJacobi(gridLevel* g, size_t l, size_t relaxations)
{
  auto t0 = std::chrono::system_clock::now();

  for (size_t r = 0; r < relaxations; r++)
    {
      double* tmp = g[l].Un; g[l].Un = g[l].U; g[l].U = tmp;
      tmp = g[l].gUn; g[l].gUn = g[l].gU; g[l].gU = tmp;
      jacobi<<<g[l].blocksPerGrid, g[l].threadsPerBlock>>>(g[l].h*g[l].h, g[l].gU, g[l].gUn, g[l].gf, g[l].N);
      checkCUDAError("Error running Jacobi Kernel");
      cudaDeviceSynchronize();
      
    }
  auto t1 = std::chrono::system_clock::now();
  smoothingTime[l] += std::chrono::duration<double>(t1-t0).count();
}

///*
__global__
void residual(double h2, double* U, double* Res, double* f, size_t N)
{
  size_t j = blockIdx.x*blockDim.x+threadIdx.x;
  if (j % N == 0 || j % N == N - 1 || j / N == 0 || j / N >= N - 1) return;
  Res[j] = f[j] + (U[j-N] + U[j+N] - 4*U[j] + U[j-1] + U[j+1]) * h2;
}


void calculateResidual(gridLevel* g, size_t l)
{
  auto t0 = std::chrono::system_clock::now();

  double h2 = 1.0 / pow(g[l].h,2);
  
  
  residual<<<g[l].blocksPerGrid, g[l].threadsPerBlock>>>(h2, g[l].gU, g[l].gRes, g[l].gf, g[l].N);
  checkCUDAError("Error running Residual Kernel");
  cudaDeviceSynchronize();

  auto t1 = std::chrono::system_clock::now();
  residualTime[l] += std::chrono::duration<double>(t1-t0).count();
}
//*/
/*
void calculateResidual(gridLevel* g, size_t l)
{
  auto t0 = std::chrono::system_clock::now();

  double h2 = 1.0 / pow(g[l].h,2);

  size_t N = g[l].N;
  double s = 0;
  for (size_t i = 1; i < N-1; i++)
    for (size_t j = i*N+1; j < i*N+N-1; j++){
      g[l].Res[j] = g[l].f[j] + (g[l].U[j-N] + g[l].U[j+N] - 4*g[l].U[j] + g[l].U[j-1] + g[l].U[j+1]) * h2;
      //s+= g[l].Res[j];

}
  //std::cerr<<"sum " << s << std::endl;
  auto t1 = std::chrono::system_clock::now();
  residualTime[l] += std::chrono::duration<double>(t1-t0).count();
}

//*/
/*
//parallel reduction adapted from https://developer.download.nvidia.com/assets/cuda/files/reduction.pdf by Mark Harris
template <unsigned int blockSize>
__device__ void warpReduce(volatile double *sdata, unsigned int tid) {
  if (blockSize >= 64) sdata[tid] += sdata[tid + 32];
  if (blockSize >= 32) sdata[tid] += sdata[tid + 16];
  if (blockSize >= 16) sdata[tid] += sdata[tid + 8];
  if (blockSize >= 8) sdata[tid] += sdata[tid + 4];
  if (blockSize >= 4) sdata[tid] += sdata[tid + 2];
  if (blockSize >= 2) sdata[tid] += sdata[tid + 1];
}
template <unsigned int blockSize>
__global__ void reduce5(double *g_idata, double *g_odata) {
  extern __shared__ double sdata[];
  unsigned int tid = threadIdx.x;
  unsigned int i = blockIdx.x*(blockSize*2) + tid;
  sdata[tid] = 0;
  sdata[tid] += g_idata[i] + g_idata[i+blockSize];
  __syncthreads();
  if (blockSize >= 2048) { if (tid < 1024) { sdata[tid] += sdata[tid + 1024]; } __syncthreads(); }
  if (blockSize >= 1024) { if (tid < 512) { sdata[tid] += sdata[tid + 512]; } __syncthreads(); }
  if (blockSize >= 512) { if (tid < 256) { sdata[tid] += sdata[tid + 256]; } __syncthreads(); }
  if (blockSize >= 256) { if (tid < 128) { sdata[tid] += sdata[tid + 128]; } __syncthreads(); }
  if (blockSize >= 128) { if (tid < 64) { sdata[tid] += sdata[tid + 64]; } __syncthreads(); }
  if (tid < 32) warpReduce<blockSize>(sdata, tid);
  if (tid == 0) g_odata[blockIdx.x] = sdata[0];
}

void runReduction(int dimGrid, int dimBlock, int smemSize, double * d_idata, double * d_odata)
{
  switch (dimBlock)
    { 
    case 2048:
      reduce5<2048><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata); break;
    case 1024:
      reduce5<1024><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata); break;
    case 512:
      reduce5<512><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata); break;
    case 256:
      reduce5<256><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata); break;
    case 128:
      reduce5<128><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata); break;
    case 64:
      reduce5< 64><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata); break;
    case 32:
      reduce5< 32><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata); break;
    case 16:
      reduce5< 16><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata); break;
    case 8:
      reduce5< 8><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata); break;
    case 4:
      reduce5< 4><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata); break;
    case 2:
      reduce5< 2><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata); break;
    case 1:
      reduce5< 1><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata); break;
    }
}
*/
__global__
void square( double* Res, int N)
{
  size_t j = blockIdx.x*blockDim.x+threadIdx.x;
  if (j >= N) return;
  //printf("res was %f\n", Res[j]);
  Res[j] = Res[j] * Res[j];
  // printf("res is %f\n", Res[j]);
}


void calculateL2Norm(gridLevel* g, size_t l)
{
  auto t0 = std::chrono::system_clock::now();

  
  
  double tmp = 0.0;
  square<<<g[l].blocksPerGrid, g[l].threadsPerBlock>>>(g[l].gRes, g[l].N * g[l].N);
  checkCUDAError("Error running Square Kernel");
  cudaDeviceSynchronize();
  //g[l].threadsPerBlock><<<g[l].blocksPerGrid, g[l].threadsPerBlock, g[l].N * g[l].N
  cudaMemcpy(g[l].Res, g[l].gRes, sizeof(double)*g[l].N * g[l].N, cudaMemcpyDeviceToHost); checkCUDAError("Error copying res squares back");
  cudaDeviceSynchronize();
  
  /*runReduction(g[l].blocksPerGrid, g[l].threadsPerBlock, g[l].N * g[l].N, g[l].gRes, g[l].gResSum);
    checkCUDAError("Error running Reduction Kernel");
  
    cudaMemcpy(g[l].ResSum, g[l].gResSum, sizeof(double)*g[l].blocksPerGrid, cudaMemcpyDeviceToHost); checkCUDAError("Error copying sum results back");
    cudaDeviceSynchronize();
    for (size_t i = 0; i < g[l].blocksPerGrid; i ++){
    tmp+= g[l].ResSum[i];
    }*/
  /*
  for (size_t i = 0; i < g[l].N*g[l].N; i++){
    g[l].Res[i] = g[l].Res[i]*g[l].Res[i];
    // std::cerr << g[l].Res[i] << std::endl;
  }
  //*/
  for (size_t i = 0; i < g[l].N*g[l].N; i++)
    tmp += g[l].Res[i];
  //*/
  g[l].L2Norm = sqrt(tmp);
  g[l].L2NormDiff = fabs(g[l].L2NormPrev - g[l].L2Norm);
  g[l].L2NormPrev = g[l].L2Norm;
  // printf("L2Norm: %.4f\n",  g[0].L2Norm);

  auto t1 = std::chrono::system_clock::now();
  L2NormTime[l] += std::chrono::duration<double>(t1-t0).count();
}

void applyRestriction(gridLevel* g, size_t l)
{
  auto t0 = std::chrono::system_clock::now();
  size_t N = g[l-1].N;
  for (size_t i = 1; i < g[l].N-1; i++)
    for (size_t j = 1; j < g[l].N-1; j++)
      g[l].f[i*g[l].N+j] = ( 1.0*( g[l-1].Res[(2*i-1)*N+2*j-1] + g[l-1].Res[(2*i-1)*N+2*j+1] + g[l-1].Res[(2*i+1)*N+2*j-1]   + g[l-1].Res[(2*i+1)*N+2*j+1] )   +
			     2.0*( g[l-1].Res[(2*i-1)*N+2*j]   + g[l-1].Res[(2*i)*N+2*j-1]   + g[l-1].Res[(2*i+1)*N+2*j]     + g[l-1].Res[(2*i)*N+2*j+1] ) +
			     4.0*( g[l-1].Res[(2*i)*N+2*j] ) ) * 0.0625;

  
  for (size_t i = 0; i < g[l].N*g[l].N; i++)// Resetting U vector for the coarser level before smoothing -- Find out if this is really necessary.
    g[l].U[i] = 0;

  auto t1 = std::chrono::system_clock::now();
  restrictionTime[l] += std::chrono::duration<double>(t1-t0).count();
}

void applyProlongation(gridLevel* g, size_t l)
{
  auto t0 = std::chrono::system_clock::now();

  for (size_t i = 1; i < g[l].N-1; i++)
    for (size_t j = 1; j < g[l].N-1; j++)
      g[l-1].U[2*i*g[l-1].N+2*j] += g[l].U[i*g[l].N+j];

  for (size_t i = 1; i < g[l].N; i++)
    for (size_t j = 1; j < g[l].N-1; j++)
      g[l-1].U[(2*i-1)*g[l-1].N+2*j] += ( g[l].U[(i-1)*g[l].N+j] + g[l].U[i*g[l].N+j] ) *0.5;

  for (size_t i = 1; i < g[l].N-1; i++)
    for (size_t j = 1; j < g[l].N; j++)
      g[l-1].U[(2*i)*g[l-1].N+2*j-1] += ( g[l].U[i*g[l].N+j-1] + g[l].U[i*g[l].N+j] ) *0.5;

  for (size_t i = 1; i < g[l].N; i++)
    for (size_t j = 1; j < g[l].N; j++)
      g[l-1].U[(2*i-1)*g[l-1].N+2*j-1] += ( g[l].U[(i-1)*g[l].N+j-1] + g[l].U[(i-1)*g[l].N+j] + g[l].U[i*g[l].N+j-1] + g[l].U[i*g[l].N+j] ) *0.25;

  auto t1 = std::chrono::system_clock::now();
  prolongTime[l] += std::chrono::duration<double>(t1-t0).count();
}

gridLevel* generateInitialConditions(size_t N0, size_t gridCount)
{
  // Default values:
  __p.nCandles = 4;
  std::vector<double> pars;
  pars.push_back(0.228162);
  pars.push_back(0.226769);
  pars.push_back(0.437278);
  pars.push_back(0.0492324);
  pars.push_back(0.65915);
  pars.push_back(0.499616);
  pars.push_back(0.59006);
  pars.push_back(0.0566329);
  pars.push_back(0.0186672);
  pars.push_back(0.894063);
  pars.push_back(0.424229);
  pars.push_back(0.047725);
  pars.push_back(0.256743);
  pars.push_back(0.754483);
  pars.push_back(0.490461);
  pars.push_back(0.0485152);

  // Allocating Timers
  smoothingTime = (double*) calloc (gridCount, sizeof(double));
  residualTime = (double*) calloc (gridCount, sizeof(double));
  restrictionTime = (double*) calloc (gridCount, sizeof(double));
  prolongTime = (double*) calloc (gridCount, sizeof(double));
  L2NormTime = (double*) calloc (gridCount, sizeof(double));

  // Allocating Grids
  gridLevel* g = (gridLevel*) malloc(sizeof(gridLevel) * gridCount);
  for (size_t i = 0; i < gridCount; i++)
    {
      g[i].N = pow(2, N0-i) + 1;
      g[i].h = 1.0/(g[i].N-1);
	  
      cudaDeviceProp deviceProp;
      cudaGetDeviceProperties(&deviceProp, 0);
      int blockSize = std::min((size_t)deviceProp.maxThreadsPerBlock,g[i].N*g[i].N);
      g[i].threadsPerBlock = (blockSize);
      g[i].blocksPerGrid = ((g[i].N*g[i].N+blockSize - 1)/blockSize);

      g[i].U   = (double*) malloc(sizeof(double) * g[i].N * g[i].N);
      g[i].Un  = (double*) malloc(sizeof(double) * g[i].N * g[i].N);
      g[i].Res = (double*) malloc(sizeof(double) * g[i].N * g[i].N);
      g[i].f   = (double*) malloc(sizeof(double) * g[i].N * g[i].N);
      g[i].ResSum = (double*) malloc(sizeof(double) * g[i].blocksPerGrid);

      cudaMalloc(&g[i].gU, sizeof(double) * g[i].N * g[i].N); checkCUDAError("Error allocating gU");
      cudaMalloc(&g[i].gUn, sizeof(double) * g[i].N * g[i].N); checkCUDAError("Error allocating gUn");
      cudaMalloc(&g[i].gRes, sizeof(double) * g[i].N * g[i].N); checkCUDAError("Error allocating gRes");
      cudaMalloc(&g[i].gf, sizeof(double) * g[i].N * g[i].N); checkCUDAError("Error allocating gf");
      cudaMalloc(&g[i].gResSum, sizeof(double) * g[i].blocksPerGrid); checkCUDAError("Error allocating gResSum");
      
      g[i].L2Norm = 0.0;
      g[i].L2NormPrev = std::numeric_limits<double>::max();
      g[i].L2NormDiff = std::numeric_limits<double>::max();
	  
	  
    }

  // Initial Guess
  for (size_t j = 0; j < g[0].N*g[0].N; j++) g[0].U[j] = 1.0;

  // Boundary Conditions
  for (size_t i = 0; i < g[0].N; i++) g[0].U[i]        = 0.0;
  for (size_t i = 0; i < g[0].N; i++) g[0].U[(g[0].N-1)*g[0].N + i] = 0.0;
  for (size_t i = 0; i < g[0].N; i++) g[0].U[i*g[0].N]        = 0.0;
  for (size_t i = 0; i < g[0].N; i++) g[0].U[i*g[0].N+g[0].N-1] = 0.0;

  // F
  for (size_t i = 0; i < g[0].N; i++)
    for (size_t j = 0; j < g[0].N; j++)
      {
	double h = 1.0/(g[0].N-1);
	double x = i*h;
	double y = j*h;

	g[0].f[i*g[0].N+j] = 0.0;

	for (size_t c = 0; c < __p.nCandles; c++)
	  {
	    double c3 = pars[c*4  + 0]; // x0
	    double c4 = pars[c*4  + 1]; // y0
	    double c1 = pars[c*4  + 2]; c1 *= 100000;// intensity
	    double c2 = pars[c*4  + 3]; c2 *= 0.01;// Width
	    g[0].f[i*g[0].N+j] += c1*exp(-(pow(c4 - y, 2) + pow(c3 - x, 2)) / c2);
	  }
      }

  return g;
}

void freeGrids(gridLevel* g, size_t gridCount)
{
  for (size_t i = 0; i < gridCount; i++)
    {
      free(g[i].U);
      free(g[i].Un);
      free(g[i].f);
      free(g[i].Res);
      cudaFree(g[i].gU);
      cudaFree(g[i].gUn);
      cudaFree(g[i].gf);
      cudaFree(g[i].gRes);
      cudaFree(g[i].gResSum);
    }
  free(g);
}

void printTimings(size_t gridCount)
{
  double* timePerGrid = (double*) calloc (sizeof(double), gridCount);
  double totalSmoothingTime = 0.0;
  double totalResidualTime = 0.0;
  double totalRestrictionTime = 0.0;
  double totalProlongTime = 0.0;
  double totalL2NormTime = 0.0;

  for (size_t i = 0; i < gridCount; i++) timePerGrid[i] = smoothingTime[i] + residualTime[i] + restrictionTime[i] + prolongTime[i] + L2NormTime[i];
  for (size_t i = 0; i < gridCount; i++) totalSmoothingTime += smoothingTime[i];
  for (size_t i = 0; i < gridCount; i++) totalResidualTime += residualTime[i];
  for (size_t i = 0; i < gridCount; i++) totalRestrictionTime += restrictionTime[i];
  for (size_t i = 0; i < gridCount; i++) totalProlongTime += prolongTime[i];
  for (size_t i = 0; i < gridCount; i++) totalL2NormTime += L2NormTime[i];

  double totalMeasured = totalSmoothingTime + totalResidualTime + totalRestrictionTime + totalProlongTime + totalL2NormTime;

  printf("   Time (s)    "); for (size_t i = 0; i < gridCount; i++) printf("Grid%lu   ", i);                    printf("   Total  \n");
  printf("-------------|-"); for (size_t i = 0; i < gridCount; i++) printf("--------"); printf("|---------\n");
  printf("Smoothing    | "); for (size_t i = 0; i < gridCount; i++) printf("%2.3f   ", smoothingTime[i]);    printf("|  %2.3f  \n", totalSmoothingTime);
  printf("Residual     | "); for (size_t i = 0; i < gridCount; i++) printf("%2.3f   ", residualTime[i]);     printf("|  %2.3f  \n", totalResidualTime);
  printf("Restriction  | "); for (size_t i = 0; i < gridCount; i++) printf("%2.3f   ", restrictionTime[i]);  printf("|  %2.3f  \n", totalRestrictionTime);
  printf("Prolongation | "); for (size_t i = 0; i < gridCount; i++) printf("%2.3f   ", prolongTime[i]);      printf("|  %2.3f  \n", totalProlongTime);
  printf("L2Norm       | "); for (size_t i = 0; i < gridCount; i++) printf("%2.3f   ", L2NormTime[i]);       printf("|  %2.3f  \n", totalL2NormTime);
  printf("-------------|-"); for (size_t i = 0; i < gridCount; i++) printf("--------"); printf("|---------\n");
  printf("Total        | "); for (size_t i = 0; i < gridCount; i++) printf("%2.3f   ", timePerGrid[i]); printf("|  %2.3f  \n", totalMeasured);
  printf("-------------|-"); for (size_t i = 0; i < gridCount; i++) printf("--------"); printf("|---------\n");
  printf("\n");
  printf("Running Time      : %.3fs\n", totalTime);
}

