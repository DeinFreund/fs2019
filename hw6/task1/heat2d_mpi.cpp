/**********************************************************************/
// A now optimized Multigrid Solver for the Heat Equation             //
// Course Material for HPCSE-II, Spring 2019, ETH Zurich              //
// Authors: Sergio Martin, Georgios Arampatzis                        //
// License: Use if you like, but give us credit.                      //
/**********************************************************************/

#include <stdio.h>
#include <math.h>
#include <limits>
#include "heat2d_mpi.hpp"
#include "string.h"
#include <chrono>
#include <mpi.h>
#include <iostream>

pointsInfo __p;

int size, rank, rank_2d;
int coord_2d[2];
int N, p, n;
MPI_Comm comm2d;
int left, right, down, up;

MPI_Datatype rowType, colType;


int main(int argc, char* argv[])
{
MPI_Init(&argc, &argv);
MPI_Comm_size(MPI_COMM_WORLD, &size);
MPI_Comm_rank(MPI_COMM_WORLD, &rank);

int ndim = 2; 
int periodic[2] = {0,0};
int dimensions[2] = {0,0};
 MPI_Dims_create(size, 2, dimensions);
MPI_Cart_create(MPI_COMM_WORLD,ndim,dimensions,periodic,1,&comm2d);
MPI_Cart_coords(comm2d,rank,ndim,coord_2d);
MPI_Cart_rank(comm2d,coord_2d,&rank_2d);
MPI_Cart_shift(comm2d, 0, 1, &left, &right);
 MPI_Cart_shift(comm2d, 1, 1, &up, &down);

 double tolerance = 1e-0; // L2 Difference Tolerance before reaching convergence.
 size_t N0 = 10; // 2^N0 + 1 elements per side

 // Multigrid parameters -- Find the best configuration!
 size_t gridCount       = 1;     // Number of Multigrid levels to use
 size_t downRelaxations = 3; // Number of Relaxations before restriction
 size_t upRelaxations   = 3;   // Number of Relaxations after prolongation

 gridLevel* g = generateInitialConditions(N0, gridCount);


 N = g[0].N;
 p = sqrt(size);
 n = N / p;

 if (rank == 0)std::cout << "Splitting the " << N << "x" << N << " matrix into " << size << " " << n << "x" << n << " matrices." << std::endl;

 MPI_Type_contiguous(n, MPI_DOUBLE, &rowType);
 MPI_Type_vector (n, 1, n, MPI_DOUBLE, &colType);
 MPI_Type_commit(&rowType);
 MPI_Type_commit(&colType);

 MPI_Datatype subMatrixType;
 MPI_Type_contiguous(n*n, MPI_DOUBLE, &subMatrixType);
 MPI_Type_commit(&subMatrixType);
 
 for (int y = 0; y < n; y++){
   for (int x = 0 ; x < n; x++){
     int u = x + coord_2d[0] * n;
     int v = y + coord_2d[1] * n;
     g[0].Ul[y*n+x] = g[0].U[v * N + u];
     g[0].Unl[y*n+x] = g[0].Un[v * N + u];
     g[0].Resl[y*n+x] = g[0].Res[v * N + u];
     g[0].fl[y*n+x] = g[0].f[v * N + u];
   }
 }
//Ul, fl
/*
  if (rank == 0)  std::cout << "Target is " << std::endl;

   for (int i = 0; i < size; i++){
     MPI_Barrier(comm2d);
     if (i == rank){
       std::cout << "Submatrix " << rank << std::endl;
       for (int y = 0; y < n; y++){
	 for (int x = 0; x < n; x++){
	   std::cout << g[0].fl[y*n+x] << " ";
	 }
	 std::cout << std::endl;
       }
     }
   }MPI_Barrier(comm2d);
*/

 int b = 0;
 auto startTime = std::chrono::system_clock::now();
  while (g[0].L2NormDiff > tolerance)  // Multigrid solver start
 //for (int b = 0; b < 3; b++) 
  {
    b++;
    /*
  if (rank == 0)  std::cout << "Step " << b << std::endl;

   for (int i = 0; i < size; i++){
     MPI_Barrier(comm2d);
     if (i == rank){
       std::cout << "Submatrix " << rank << std::endl;
       for (int y = 0; y < n; y++){
	 for (int x = 0; x < n; x++){
	   std::cout << g[0].Ul[y*n+x] << " ";
	 }
	 std::cout << std::endl;
       }
     }
   }MPI_Barrier(comm2d);
   // */
  

   applyJacobi(g, 0, downRelaxations); // Relaxing the finest grid first
  calculateResidual(g, 0); // Calculating Initial Residual

  calculateL2Norm(g, 0); // Calculating Residual L2 Norm

  MPI_Bcast(&g[0].L2NormDiff, 1, MPI_DOUBLE, 0, comm2d);
  
 }  // Multigrid solver en

 MPI_Barrier(comm2d);
     
 
 auto endTime = std::chrono::system_clock::now();
 totalTime = std::chrono::duration<double>(endTime-startTime).count();
 if (rank == 0){ 
   printTimings(gridCount);
   printf("L2Norm: %.4f\n",  g[0].L2Norm);
 }
 freeGrids(g, gridCount);
 return 0;
}

void applyJacobi(gridLevel* g, size_t l, size_t relaxations)
{
 auto t0 = std::chrono::system_clock::now();

 double h1 = 0.25;
 double h2 = g[l].h*g[l].h;
 for (size_t r = 0; r < relaxations; r++)
 {
  double* tmp = g[l].Unl; g[l].Unl = g[l].Ul; g[l].Ul = tmp;

  bool hasUp, hasDown, hasLeft, hasRight;
  hasUp = coord_2d[1] > 0;
  hasDown = coord_2d[1] < p -1;
  hasLeft = coord_2d[0] > 0;
  hasRight = coord_2d[0] < p - 1;

  //if (r == 0) std::cout << "Rank " << rank << " has left, right, up, down: " << hasLeft << " " << hasRight << " " << hasUp << " " << hasDown << std::endl;
  
  MPI_Request reqs[8];
  int rcnt = 0;
  if (hasLeft) MPI_Irecv(g[l].left, n, MPI_DOUBLE, left, 1, comm2d, &(reqs[rcnt++]));
  if (hasRight) MPI_Irecv(g[l].right, n, MPI_DOUBLE, right, 2, comm2d, &(reqs[rcnt++]));
  if (hasUp) MPI_Irecv(g[l].up, n, MPI_DOUBLE, up , 3, comm2d, &(reqs[rcnt++]));
  if (hasDown)   MPI_Irecv(g[l].down, n, MPI_DOUBLE, down, 4, comm2d, &(reqs[rcnt++]));

  if (hasLeft) MPI_Isend(g[l].Unl, 1, colType, left, 2,comm2d, &(reqs[rcnt++]));
  if (hasRight) MPI_Isend(g[l].Unl + (n-1), 1, colType, right, 1, comm2d, &(reqs[rcnt++]));
  if (hasUp) MPI_Isend(g[l].Unl, 1, rowType, up, 4, comm2d, &(reqs[rcnt++]));
  if (hasDown) MPI_Isend(g[l].Unl + n * (n-1), 1, rowType, down, 3, comm2d, &(reqs[rcnt++]));

  for (size_t i = 1; i < n-1; i++)
   for (size_t k = 1; k < n-1; k++){
     int j = i * n + k;
     g[l].Ul[j] = (g[l].Unl[j - n] + g[l].Unl[j + n] + g[l].Unl[j - 1] + g[l].Unl[j + 1] + g[l].fl[j] * h2) * 0.25;
   }

  MPI_Status stats[8];
  MPI_Waitall(rcnt, reqs, stats);

   //left
   if (hasLeft)
   for (int j = n; j < n*n-n; j+= n){
     g[l].Ul[j] = (g[l].Unl[j - n] + g[l].Unl[j + n] + g[l].left[j/n] + g[l].Unl[j + 1] + g[l].fl[j] * h2) * 0.25;
   }
   //right
   if (hasRight)
   for (int j = n-1+n; j < n*n-n; j+= n){
     g[l].Ul[j] = (g[l].Unl[j - n] + g[l].Unl[j + n] + g[l].Unl[j - 1] + g[l].right[j/n] + g[l].fl[j] * h2) * 0.25;
   }
   //up
   if (hasUp)
   for (int j = 1; j < n-1; j++){
     g[l].Ul[j] = (g[l].up[j] + g[l].Unl[j + n] + g[l].Unl[j - 1] + g[l].Unl[j + 1] + g[l].fl[j] * h2) * 0.25;
   }
   //down
   if (hasDown)
     for (int j = n*n-n+1; j < n*n-1; j++){
     g[l].Ul[j] = (g[l].Unl[j - n] + g[l].down[j - n*n + n] + g[l].Unl[j - 1] + g[l].Unl[j + 1] + g[l].fl[j] * h2) * 0.25;
   }
   if (hasUp && hasRight){
     int j = n - 1;
     g[l].Ul[j] = (g[l].up[n-1] + g[l].Unl[j + n] + g[l].Unl[j - 1] + g[l].right[0] + g[l].fl[j] * h2) * 0.25;
   }
   if (hasRight && hasDown){
     int j = n * n - 1;
     g[l].Ul[j] = (g[l].Unl[j - n] + g[l].down[j - n*n + n] + g[l].Unl[j - 1] + g[l].right[n - 1] + g[l].fl[j] * h2) * 0.25;
   }
   if (hasDown && hasLeft){
     int j = n * n - n;
     g[l].Ul[j] = (g[l].Unl[j - n] + g[l].down[j - n*n + n] + g[l].left[n-1] + g[l].Unl[j + 1] + g[l].fl[j] * h2) * 0.25;
   }
   if (hasLeft && hasUp){
     int j = 0;
     g[l].Ul[j] = (g[l].up[0] + g[l].Unl[j + n] + g[l].left[j/n] + g[l].Unl[j + 1] + g[l].fl[j] * h2) * 0.25;
   }

 }
 auto t1 = std::chrono::system_clock::now();
 smoothingTime[l] += std::chrono::duration<double>(t1-t0).count();
}

void calculateResidual(gridLevel* g, size_t l)
{
 auto t0 = std::chrono::system_clock::now();

 double h2 = 1.0 / pow(g[l].h,2);



  bool hasUp, hasDown, hasLeft, hasRight;
  hasUp = coord_2d[1] > 0;
  hasDown = coord_2d[1] < p -1;
  hasLeft = coord_2d[0] > 0;
  hasRight = coord_2d[0] < p - 1;

  //if (r == 0) std::cout << "Rank " << rank << " has left, right, up, down: " << hasLeft << " " << hasRight << " " << hasUp << " " << hasDown << std::endl;
  
  MPI_Request reqs[8];
  int rcnt = 0;

  if (hasLeft) MPI_Irecv(g[l].left, n, MPI_DOUBLE, left, 1, comm2d, &(reqs[rcnt++]));
  if (hasRight) MPI_Irecv(g[l].right, n, MPI_DOUBLE, right, 2, comm2d, &(reqs[rcnt++]));
  if (hasUp) MPI_Irecv(g[l].up, n, MPI_DOUBLE, up , 3, comm2d, &(reqs[rcnt++]));
  if (hasDown)   MPI_Irecv(g[l].down, n, MPI_DOUBLE, down, 4, comm2d, &(reqs[rcnt++]));

  if (hasLeft) MPI_Isend(g[l].Ul, 1, colType, left, 2,comm2d, &(reqs[rcnt++]));
  if (hasRight) MPI_Isend(g[l].Ul + (n-1), 1, colType, right, 1, comm2d, &(reqs[rcnt++]));
  if (hasUp) MPI_Isend(g[l].Ul, 1, rowType, up, 4, comm2d, &(reqs[rcnt++]));
  if (hasDown) MPI_Isend(g[l].Ul + n * (n-1), 1, rowType, down, 3, comm2d, &(reqs[rcnt++]));

  for (size_t i = 1; i < n-1; i++)
   for (size_t k = 1; k < n-1; k++){
     int j = i * n + k;
     g[l].Resl[j] = g[l].fl[j] + (g[l].Ul[j - n] + g[l].Ul[j + n] + g[l].Ul[j - 1] + g[l].Ul[j + 1] - 4 * g[l].Ul[j]) * h2;
   }

  MPI_Status stats[8];
  MPI_Waitall(rcnt, reqs, stats);

   //left
   if (hasLeft)
   for (int j = n; j < n*n-n; j+= n){
     g[l].Resl[j] = g[l].fl[j] + (g[l].Ul[j - n] + g[l].Ul[j + n] + g[l].left[j/n] + g[l].Ul[j + 1] - 4 * g[l].Ul[j]) * h2;
   }
   //right
   if (hasRight)
   for (int j = n-1+n; j < n*n-n; j+= n){
     g[l].Resl[j] = g[l].fl[j] + (g[l].Ul[j - n] + g[l].Ul[j + n] + g[l].Ul[j - 1] + g[l].right[j/n] - 4 * g[l].Ul[j]) * h2;
   }
   //up
   if (hasUp)
   for (int j = 1; j < n-1; j++){
     g[l].Resl[j] = g[l].fl[j] + (g[l].up[j] + g[l].Ul[j + n] + g[l].Ul[j - 1] + g[l].Ul[j + 1] - 4 * g[l].Ul[j]) * h2;
   }
   //down
   if (hasDown)
     for (int j = n*n-n+1; j < n*n-1; j++){
     g[l].Resl[j] = g[l].fl[j] + (g[l].Ul[j - n] + g[l].down[j - n*n + n] + g[l].Ul[j - 1] + g[l].Ul[j + 1] - 4 * g[l].Ul[j]) * h2;
   }
   if (hasUp && hasRight){
     int j = n - 1;
     g[l].Resl[j] = g[l].fl[j] + (g[l].up[n-1] + g[l].Ul[j + n] + g[l].Ul[j - 1] + g[l].right[0] - 4 * g[l].Ul[j]) * h2;
   }
   if (hasRight && hasDown){
     int j = n * n - 1;
     g[l].Resl[j] = g[l].fl[j] + (g[l].Ul[j - n] + g[l].down[j - n*n + n] + g[l].Ul[j - 1] + g[l].right[n - 1] - 4 * g[l].Ul[j]) * h2;
   }
   if (hasDown && hasLeft){
     int j = n * n - n;
     g[l].Resl[j] = g[l].fl[j] + (g[l].Ul[j - n] + g[l].down[j - n*n + n] + g[l].left[n-1] + g[l].Ul[j + 1] - 4 * g[l].Ul[j]) * h2;
   }
   if (hasLeft && hasUp){
     int j = 0;
     g[l].Resl[j] = g[l].fl[j] + (g[l].up[0] + g[l].Ul[j + n] + g[l].left[j/n] + g[l].Ul[j + 1] - 4 * g[l].Ul[j]) * h2;
   }


 auto t1 = std::chrono::system_clock::now();
 residualTime[l] += std::chrono::duration<double>(t1-t0).count();
}

void calculateL2Norm(gridLevel* g, size_t l)
{
 auto t0 = std::chrono::system_clock::now();

 double tmp = 0.0;

 for (size_t i = 0; i <n* n; i++)
   g[l].Resl[i] = g[l].Resl[i]*g[l].Resl[i];

 for (size_t i = 0; i < n*n; i++)
   tmp += g[l].Resl[i];
 double total = 0;
 MPI_Reduce(&tmp, &total, 1, MPI_DOUBLE, MPI_SUM,0, comm2d);
 
 g[l].L2Norm = sqrt(total);
 g[l].L2NormDiff = fabs(g[l].L2NormPrev - g[l].L2Norm);
 g[l].L2NormPrev = g[l].L2Norm;
// printf("L2Norm: %.4f\n",  g[0].L2Norm);

 auto t1 = std::chrono::system_clock::now();
 L2NormTime[l] += std::chrono::duration<double>(t1-t0).count();
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

    g[i].U = (double * ) malloc(sizeof(double) * g[i].N * g[i].N);
    g[i].Un = (double * ) malloc(sizeof(double) * g[i].N * g[i].N);
    g[i].Res = (double * ) malloc(sizeof(double) * g[i].N * g[i].N);
    g[i].f = (double * ) malloc(sizeof(double) * g[i].N * g[i].N);
    
    g[i].Ul = (double * ) malloc(sizeof(double) * g[i].N * g[i].N / size);
    g[i].Unl = (double * ) malloc(sizeof(double) * g[i].N * g[i].N / size);
    g[i].Resl = (double * ) malloc(sizeof(double) * g[i].N * g[i].N / size);
    g[i].fl = (double * ) malloc(sizeof(double) * g[i].N * g[i].N / size);
    
    g[i].up = (double * ) malloc(sizeof(double) * g[i].N);
    g[i].down = (double * ) malloc(sizeof(double) * g[i].N);
    g[i].left = (double * ) malloc(sizeof(double) * g[i].N);
    g[i].right = (double * ) malloc(sizeof(double) * g[i].N);

  g[i].L2Norm = 0.0;
  g[i].L2NormPrev = std::numeric_limits<double>::max();
  g[i].L2NormDiff = std::numeric_limits<double>::max();
 }

  // Initial Guess
  for (size_t j = 0; j<g[0].N * g[0].N; j++) g[0].U[j] = 1.0;

  // Boundary Conditions
  for (size_t i = 0; i<g[0].N; i++) g[0].U[i] = 0.0;
  for (size_t i = 0; i<g[0].N; i++) g[0].U[(g[0].N - 1) * g[0].N + i] = 0.0;
  for (size_t i = 0; i<g[0].N; i++) g[0].U[i * g[0].N] = 0.0;
  for (size_t i = 0; i<g[0].N; i++) g[0].U[i * g[0].N + g[0].N - 1] = 0.0;

  // F
  for (size_t i = 0; i<g[0].N; i++)
    for (size_t j = 0; j<g[0].N; j++) {
      double h = 1.0 / (g[0].N - 1);
      double x = i * h;
      double y = j * h;

      g[0].f[i * g[0].N + j] = 0.0;

      for (size_t c = 0; c<__p.nCandles; c++) {
        double c3 = pars[c * 4 + 0]; // x0
        double c4 = pars[c * 4 + 1]; // y0
        double c1 = pars[c * 4 + 2];
        c1 *= 100000; // intensity
        double c2 = pars[c * 4 + 3];
        c2 *= 0.01; // Width
        g[0].f[i * g[0].N + j] += c1 * exp(-(pow(c4 - y, 2) + pow(c3 - x, 2)) / c2);
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
free(g[i].Ul);
  free(g[i].Unl);
  free(g[i].Resl);
  free(g[i].left);
  free(g[i].right);
  free(g[i].up);
  free(g[i].down);
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

