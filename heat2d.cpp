#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "auxiliar/auxiliar.hpp"

#include <iostream> //debugging
using namespace std;

class GridLevel
{
public:
  size_t N; // Number of points per dimension in the grid level
  double h; // DeltaX = DeltaY, the distance between points in the discretized [0,1]x[0,1] domain
  double** f; // Right hand side (external heat sources)
  double** U; // Main grid for Jacobi
  double** Un; // Previous' step grid
  double** Res; // Residual Grid
};

int blockWidth = 16;
int blockHeight = 16; 

//helper to make allocations nicer
template<typename T>
T** alloc2D(size_t dim1, size_t dim2){
  while (dim1%blockWidth != 0) dim1++;
  while (dim2%blockHeight != 0) dim2++;
  T** ret = (T**)malloc(sizeof(T*) * dim1 + sizeof(T) * dim1 * dim2);
  T* start = (T*)(ret + dim1);
  for (int i = 0; i < dim1; i++) ret[i] = start + i * dim2;
  for (int i = 0; i < dim1; i++) for (int j = 0; j < dim2; j++) ret[i][j] = 0;
  return ret;
}

void heat2DSolver(Heat2DSetup& s)
{
  // Multigrid parameters -- Find the best configuration!
  s.setGridCount(6);     // Number of Multigrid levels to use
  s.downRelaxations = 3; // Number of Relaxations before restriction
  s.upRelaxations   = 1;   // Number of Relaxations after prolongation

  // Allocating Grids -- Is there a better way to allocate these grids?
  GridLevel* g = (GridLevel*) calloc(sizeof(GridLevel), s.gridCount);
  for (int i = 0; i < s.gridCount; i++)
    {
      g[i].N = pow(2, s.N0-i) + 1;
      g[i].h = 1.0/(g[i].N-1);

      if (true)
	{
	  //make arrays contiguous in memory
	  g[i].U   = alloc2D<double>(g[i].N, g[i].N);
	  g[i].Un  = alloc2D<double>(g[i].N, g[i].N);
	  g[i].Res = alloc2D<double>(g[i].N, g[i].N);
	  g[i].f   = alloc2D<double>(g[i].N, g[i].N);
	}else
	{
	  g[i].U   = (double**) calloc (sizeof(double*), g[i].N);
	  g[i].Un  = (double**) calloc (sizeof(double*), g[i].N);
	  g[i].Res = (double**) calloc (sizeof(double*), g[i].N);
	  g[i].f   = (double**) calloc (sizeof(double*), g[i].N);

	  for (int j = 0; j < g[i].N ; j++)
	    {
	      g[i].U[j]   = (double*) calloc (sizeof(double), g[i].N);
	      g[i].Un[j]  = (double*) calloc (sizeof(double), g[i].N);
	      g[i].Res[j] = (double*) calloc (sizeof(double), g[i].N);
	      g[i].f[j]   = (double*) calloc (sizeof(double), g[i].N);
	    }}
    }

  // Setting up problem.
  for (int i = 0; i < s.N; i++) for (int j = 0; j < s.N; j++) g[0].U[i][j] = s.getInitial(i,j);
  for (int i = 0; i < s.N; i++) for (int j = 0; j < s.N; j++) g[0].f[i][j] = s.getRHS(i,j);

  while (s.L2NormDiff > s.tolerance)  // Multigrid solver start
    {
      s.applyJacobi_(g, 0, s.downRelaxations); // Relaxing the finest grid first
      s.calculateResidual_(g, 0); // Calculating Initial Residual

      for (int grid = 1; grid < s.gridCount; grid++) // Going down the V-Cycle
	{
	  s.applyRestriction_(g, grid); // Restricting the residual to the coarser grid's solution vector (f)
	  s.applyJacobi_(g, grid, s.downRelaxations); // Smoothing coarser level
	  s.calculateResidual_(g, grid); // Calculating Coarse Grid Residual
	}

      for (int grid = s.gridCount-1; grid > 0; grid--) // Going up the V-Cycle
	{
	  s.applyProlongation_(g, grid); // Prolonging solution for coarser level up to finer level
	  s.applyJacobi_(g, grid, s.upRelaxations); // Smoothing finer level
	}

      s.calculateL2Norm_(g, 0); // Calculating Residual L2 Norm
    }  // Multigrid solver end

  // Saving solution before returning
  for (int i = 0; i < g[0].N; i++) for (int j = 0; j < g[0].N; j++) s.saveSolution(i, j, g[0].U[i][j]);
}

void applyJacobi(GridLevel* g, int l, int relaxations)
{
  for (int r = 0; r < relaxations; r++)
    {
      swap(g[l].Un, g[l].U);

      double fac = pow(g[l].h,2)/4;
      int N = g[l].N-1;
      for (int x = 1; x < N; x += blockWidth)
	for (int y = 1; y < N; y += blockHeight)
	  for (int i = x; i < min(N, x + blockWidth); i++)
	    for (int j = y; j < min (N, y + blockHeight); j++) // Perform a Jacobi Iteration
	      g[l].U[i][j] = (g[l].Un[i-1][j] + g[l].Un[i+1][j] + g[l].Un[i][j-1] + g[l].Un[i][j+1])* 0.25 + g[l].f[i][j]*fac;
    }
}

void calculateResidual(GridLevel* g, int l)
{
  int N = g[l].N-1;
  double fac = 1 / (pow(g[l].h,2));
  for (int x = 1; x < N; x += blockWidth)
    for (int y = 1; y < N; y += blockHeight)
      for (int i = x; i < min(N, x + blockWidth); i++)
	for (int j = y; j < min (N, y + blockHeight); j++) 
	  g[l].Res[i][j] = g[l].f[i][j] + (g[l].U[i-1][j] + g[l].U[i+1][j] - 4*g[l].U[i][j] + g[l].U[i][j-1] + g[l].U[i][j+1]) * fac;
}

double calculateL2Norm(GridLevel* g, int l)
{
  int N = g[l].N-1;
  double tmp = 0.0;
  for (int i = 0; i < g[l].N; i++)
    for (int j = 0; j < g[l].N; j++)
      tmp += g[l].Res[i][j]*g[l].Res[i][j];

  return sqrt(tmp);
}

void applyRestriction(GridLevel* g, int l)
{
  int N = g[l].N-1;
  for (int x = 1; x < N; x += blockWidth)
    for (int y = 1; y < N; y += blockHeight)
      for (int i = x; i < min(N, x + blockWidth); i++)
	for (int j = y; j < min (N, y + blockHeight); j++) 
	  g[l].f[i][j] = ( 1.0*( g[l-1].Res[2*i-1][2*j-1] + g[l-1].Res[2*i-1][2*j+1] + g[l-1].Res[2*i+1][2*j-1]   + g[l-1].Res[2*i+1][2*j+1] ) +
			   2.0*( g[l-1].Res[2*i-1][2*j]   + g[l-1].Res[2*i][2*j-1]   + g[l-1].Res[2*i+1][2*j]     + g[l-1].Res[2*i][2*j+1]   ) +
			   4.0*( g[l-1].Res[2*i][2*j] ) ) * 0.0625;

  for (int i = 0; i < g[l].N; i++)
    for (int j = 0; j < g[l].N; j++)
      g[l].U[i][j] = 0;
}

void applyProlongation(GridLevel* g, int l)
{
  for (int j = 1; j < g[l].N-1; j++)
    for (int i = 1; i < g[l].N-1; i++)
      g[l-1].Un[2*i][2*j] = g[l].U[i][j];

  for (int j = 1; j < g[l].N-1; j++)
    for (int i = 1; i < g[l].N; i++)
      g[l-1].Un[2*i-1][2*j] = ( g[l].U[i-1][j] + g[l].U[i][j] ) * 0.5;

  for (int j = 1; j < g[l].N; j++)
    for (int i = 1; i < g[l].N-1; i++)
      g[l-1].Un[2*i][2*j-1] = ( g[l].U[i][j-1] + g[l].U[i][j] )  * 0.5;

  for (int j = 1; j < g[l].N; j++)
    for (int i = 1; i < g[l].N; i++)
      g[l-1].Un[2*i-1][2*j-1] = ( g[l].U[i-1][j-1] + g[l].U[i-1][j] + g[l].U[i][j-1] + g[l].U[i][j] ) * 0.25;

  for (int j = 0; j < g[l-1].N; j++)
    for (int i = 0; i < g[l-1].N; i++)
      g[l-1].U[i][j] += g[l-1].Un[i][j];
}
