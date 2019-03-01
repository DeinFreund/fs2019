#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "auxiliar/auxiliar.hpp"

#include <iostream> //debugging
using namespace std;
using numeric=long long;
using numericF= long long;

class GridLevel
{
public:
  size_t N; // Number of points per dimension in the grid level
  numeric h; // 1/DeltaX = 1/DeltaY, the distance between points in the discretized [0,1]x[0,1] domain
  numericF** f; // Right hand side (external heat sources)
  numeric** U; // Main grid for Jacobi
  numeric** Un; // Previous' step grid
  numericF** Res; // Residual Grid
};

int blockWidth = 32;
int blockHeight = 8;

//the following values represent the decimal places of the fixed point numbers
int* shifts;// = 34;
int* fReduction;// = 6;
int l2NormShift = 7;

numericF intSquare(numericF in, int l){
  numericF t = in >> (shifts[l] - fReduction[l]);
  return t * t;
}

numeric logNum(numeric in){
  int i = 0;
  for (i = 0; in > 0; i++) in >>= 1;
  return i - 1;
}

numeric convertU(double in, int l){
  numeric ret = ((numeric)(in * (1LL << shifts[l])));
  //cerr << in << " -> " << ret << " ";
  return ret;
}

double revertF(numericF in, int l){
  return in / ((double) (1LL << (shifts[l] - fReduction[l])));
}
numericF convertF(double in, int l){
  numericF ret = ((numericF)(in * (1LL << (shifts[l] - fReduction[l]))));
  //cerr << in << " -> " << ret << " ";
  return ret;
}

double revertU(numeric in, int l){
  return in / ((double) (1LL << shifts[l]));
}
void printMinMax(GridLevel * g){
  for (int l = 0; l < 6; l++){
  double minv = 1e8;
  double maxv = -1e8;

 
  for (int j = 0; j < g[l].N; j++)
    for (int i = 0; i < g[l].N; i++){
      minv = min(revertU(g[l].U[i][j], l), minv);
      maxv = max(revertU(g[l].U[i][j], l), maxv);
    }
  cerr << "= Level " << l << " =" << endl;
  cerr << "minu: " << minv << endl;
  cerr << "maxu: " << maxv << endl;
  for (int j = 0; j < g[l].N; j++)
    for (int i = 0; i < g[l].N; i++){
      minv = min(revertF(g[l].f[i][j],l), minv);
      maxv = max(revertF(g[l].f[i][j],l), maxv);
    }
  cerr << "minf: " << minv << endl;
  cerr << "maxf: " << maxv << endl;
  }
}
void adjustShift(GridLevel * g, int grids){
  for (int l = 0; l < grids; l++){
  double minv = 1e8;
  double maxv = -1e8;

 
  for (int j = 0; j < g[l].N; j++)
    for (int i = 0; i < g[l].N; i++){
      minv = min(revertU(g[l].U[i][j],l), minv);
      maxv = max(abs(revertU(g[l].U[i][j],l)), maxv);
    }
  cerr << "= Level " << l << " =" << endl;
  cerr << "minu: " << minv << endl;
  cerr << "maxu: " << maxv << endl;
  int newUShift = 28 - log(maxv)/log(2);
  cerr << "largest U is " << maxv << " use " << newUShift << " shift" << endl;
  for (int j = 0; j < g[l].N; j++)
    for (int i = 0; i < g[l].N; i++){
      minv = min(revertF(g[l].f[i][j],l), minv);
      maxv = max(abs(revertF(g[l].f[i][j],l)), maxv);
    }
  cerr << "minf: " << minv << endl;
  cerr << "maxf: " << maxv << endl;
  maxv /= g[l].h * g[l].h;
  int newFShift = 28 - log(maxv)/log(2);
  cerr << "largest F is " << maxv << " use " << newFShift << " shift" << endl;

  int extraShift = newUShift - shifts[l];
  /*for (int j = 0; j < g[l].N; j++)
    for (int i = 0; i < g[l].N; i++){
      g[l].U[i][j] <<= extraShift;
    }
    }*/
  }
}


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
      g[i].h = (g[i].N-1);

      //make arrays contiguous in memory
      g[i].U   = alloc2D<numeric>(g[i].N, g[i].N);
      g[i].Un  = alloc2D<numeric>(g[i].N, g[i].N);
      g[i].Res = alloc2D<numericF>(g[i].N, g[i].N);
      g[i].f   = alloc2D<numericF>(g[i].N, g[i].N);
    }
  shifts = (int*)calloc(sizeof(int), s.gridCount);
  fReduction = (int*)calloc(sizeof(int), s.gridCount);

  for (int i = 0; i < s.gridCount; i++){
    shifts[i] = 34;
    fReduction[i]= 6;
  }
  
  // Setting up problem.
  double largestU = 0;
  for (int i = 0; i < s.N; i++) for (int j = 0; j < s.N; j++) largestU = max(largestU, s.getRHS(i,j));
  largestU /= g[s.gridCount - 1].h * g[s.gridCount - 1].h;
  cerr << "largest F is " << largestU << " use " << 31 - log(largestU)/log(2) << " shift" << endl;
  largestU = 0;
  for (int i = 0; i < s.N; i++) for (int j = 0; j < s.N; j++) largestU = max(largestU, s.getInitial(i,j));
  cerr << "largest U is " << largestU << " use " << 31 - log(largestU)/log(2) << " shift" << endl;
  for (int i = 0; i < s.N; i++) for (int j = 0; j < s.N; j++) g[0].U[i][j] = convertU(s.getInitial(i,j),0);
  for (int i = 0; i < s.N; i++) for (int j = 0; j < s.N; j++) g[0].f[i][j] = convertF(s.getRHS(i,j),0);

  adjustShift(g, s.gridCount);
  int blub = 0;
  double minDiff = s.L2NormDiff;
  double prevDiff = s.L2NormDiff;
  int divergesCount = 0;
  
  while (s.L2NormDiff > 1e-6)  // Multigrid solver start
    {
      minDiff = min(s.L2NormDiff, minDiff);
      if (s.L2NormDiff > 2*minDiff && s.L2NormDiff > prevDiff){
	divergesCount ++;//diverging, numeric limit reached
      }else if (s.L2NormDiff < 1e-7 + minDiff){
	divergesCount = 0;
      }
      else
	{
	divergesCount = max(0, divergesCount - 1);
      }
      if (divergesCount > 50){
	cerr << "Solver is no longer converging. The best norm difference was " << minDiff << endl;
	break;
      }
      prevDiff = s.L2NormDiff;
      s.applyJacobi_(g, 0, s.downRelaxations); // Relaxing the finest grid first
      s.calculateResidual_(g, 0); // Calculating Initial Residual

      for (int grid = 1; grid < s.gridCount; grid++) // Going down the V-Cycle
	{
	  s.applyRestriction_(g, grid); // Restricting the residual to the coarser grid's solution vector (f)
	  s.applyJacobi_(g, grid, s.downRelaxations); // Smoothing coarser level
	  s.calculateResidual_(g, grid); // Calculating Coarse Grid Residual
	}
      if (blub++ < 50){
	cout << "diff: " << s.L2NormDiff << endl;

	adjustShift(g, s.gridCount);
      }

      for (int grid = s.gridCount-1; grid > 0; grid--) // Going up the V-Cycle
	{
	  s.applyProlongation_(g, grid); // Prolonging solution for coarser level up to finer level
	  s.applyJacobi_(g, grid, s.upRelaxations); // Smoothing finer level
	}

      //       printMinMax(g);
      s.calculateL2Norm_(g, 0); // Calculating Residual L2 Norm
    }  // Multigrid solver end

  // Saving solution before returning
  for (int i = 0; i < g[0].N; i++) for (int j = 0; j < g[0].N; j++) s.saveSolution(i, j, revertU(g[0].U[i][j], 0));
}


void applyJacobi(GridLevel* g, int l, int relaxations)
{
  for (int r = 0; r < relaxations; r++)
    {
      swap(g[l].Un, g[l].U);

      int N = g[l].N-1;
      numeric fac =  (g[l].h * g[l].h);

      numeric logn = logNum(fac) - fReduction[l] + 1;
      for (int x = 1; x < N; x += blockWidth)
	for (int y = 1; y < N; y += blockHeight)
	  for (int i = x; i < min(N, x + blockWidth); i++)
	    for (int j = y; j < min (N, y + blockHeight); j++) // Perform a Jacobi Iteration
	      g[l].U[i][j] = ((((g[l].Un[i-1][j] + g[l].Un[i+1][j]) >> 1) + ((g[l].Un[i][j-1] + g[l].Un[i][j+1]) >> 1)) + (numeric)(g[l].f[i][j] >> logn)) >> 1;
    }
}

void calculateResidual(GridLevel* g, int l)
{
  int N = g[l].N-1;
  numeric fac =  (g[l].h * g[l].h);
  numeric logn = logNum(fac) - fReduction[l];
  for (int x = 1; x < N; x += blockWidth)
    for (int y = 1; y < N; y += blockHeight)
      for (int i = x; i < min(N, x + blockWidth); i++)
	for (int j = y; j < min (N, y + blockHeight); j++) 
	  g[l].Res[i][j] = g[l].f[i][j] + (((numericF)(g[l].U[i-1][j] + g[l].U[i+1][j] - (((numericF)g[l].U[i][j]) << 2) + g[l].U[i][j-1] + g[l].U[i][j+1])) << (logn) );
}

  double calculateL2Norm(GridLevel* g, int l)
  {
    int N = g[l].N-1;
    int s = l2NormShift;
    long long tmp = 0.0;
    for (int i = 0; i < g[l].N; i++)
      for (int j = 0; j < g[l].N; j++){
	tmp += intSquare(g[l].Res[i][j] << s, l);
	//cerr << g[l].Res[i][j] << " | " << intSquare(g[l].Res[i][j]) << " | " << tmp << endl;
      }

    return sqrt((tmp) / (double)(1 << (s << 1)));
  }

    void applyRestriction(GridLevel* g, int l)
    {
      int N = g[l].N-1;
      for (int x = 1; x < N; x += blockWidth)
	for (int y = 1; y < N; y += blockHeight)
	  for (int i = x; i < min(N, x + blockWidth); i++)
	    for (int j = y; j < min (N, y + blockHeight); j++) 
	      g[l].f[i][j] = (( ( ((g[l-1].Res[(i<<1)-1][(j<<1)-1] + g[l-1].Res[(i<<1)-1][(j<<1)+1])>> 1) + ((g[l-1].Res[(i<<1)+1][(j<<1)-1]   + g[l-1].Res[(i<<1)+1][(j<<1)+1]) >> 1) ) >> 3) +
			      (((( g[l-1].Res[(i<<1)-1][(j<<1)]   + g[l-1].Res[(i<<1)][(j<<1)-1]) >> 1)   + ((g[l-1].Res[(i<<1)+1][(j<<1)]     + g[l-1].Res[(i<<1)][(j<<1)+1])>> 1)   ) >> 2) +
			      (( g[l-1].Res[(i<<1)][(j<<1)] )>>2) );

      for (int i = 0; i < g[l].N; i++)
	for (int j = 0; j < g[l].N; j++)
	  g[l].U[i][j] = 0;
    }

void applyProlongation(GridLevel* g, int l)
{
  for (int j = 1; j < g[l].N-1; j++)
    for (int i = 1; i < g[l].N-1; i++)
      g[l-1].Un[(i<<1)][(j<<1)] = g[l].U[i][j];

  for (int j = 1; j < g[l].N-1; j++)
    for (int i = 1; i < g[l].N; i++)
      g[l-1].Un[(i<<1)-1][(j<<1)] = ( g[l].U[i-1][j] + g[l].U[i][j] ) >> 1;

  for (int j = 1; j < g[l].N; j++)
    for (int i = 1; i < g[l].N-1; i++)
      g[l-1].Un[(i<<1)][(j<<1)-1] = ( g[l].U[i][j-1] + g[l].U[i][j] ) >> 1;

  for (int j = 1; j < g[l].N; j++)
    for (int i = 1; i < g[l].N; i++)
      g[l-1].Un[(i<<1)-1][(j<<1)-1] = ((( g[l].U[i-1][j-1] + g[l].U[i-1][j] ) >> 1 )+ ((g[l].U[i][j-1] + g[l].U[i][j] ) >> 1 )) >> 1;

  for (int j = 0; j < g[l-1].N; j++)
    for (int i = 0; i < g[l-1].N; i++)
      g[l-1].U[i][j] += g[l-1].Un[i][j];
}
