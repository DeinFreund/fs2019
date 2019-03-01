#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "auxiliar/auxiliar.hpp"

#include <iostream> //debugging
using namespace std;
using numeric=long ;
using numericF= long;

class GridLevel
{
public:
  size_t N; // Number of points per dimension in the grid level
  numeric h; // 1/DeltaX = 1/DeltaY, the distance between points in the discretized [0,1]x[0,1] domain
  numericF** f; // Right hand side (external heat sources)
  numeric** U; // Main grid for Jacobi
  numeric** Un; // Previous' step grid
  numericF** Res; // Residual Grid

  
  double dh; // 1/DeltaX = 1/DeltaY, the distance between points in the discretized [0,1]x[0,1] domain
  double** df; // Right hand side (external heat sources)
  double** dU; // Main grid for Jacobi
  double** dUn; // Previous' step grid
  double** dRes; // Residual Grid
};

int blockWidth = 32;
int blockHeight = 8;

//the following values represent the decimal places of the fixed point numbers
int* shifts;// = 34;
int* fReduction;// = 6;
int l2NormShift = 7;

long long intSquare(numericF in, int l){
  long long t = in >> (shifts[l] - fReduction[l]);
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

int getInitialShift(GridLevel * g, Heat2DSetup& s){

  int l = 0;
  double maxv = 1;
  for (int j = 0; j < g[l].N; j++)
    for (int i = 0; i < g[l].N; i++){
      maxv = max(abs(s.getRHS(i,j)), maxv);
    }
  int newFShift = 30 - log(maxv)/log(2);
  cerr << "largest F/Res is " << maxv << " use " << newFShift << " F shift" << endl;
  return newFShift;
}

void adjustShift(GridLevel * g, int grids){
  int lowerUShift = 100;
  int lowerFShift = 100;
  int bestShift = 100;
  for (int l = 0; l < grids; l++){
    double minv = 1e8;
    double maxv = 1;

 
    for (int j = 0; j < g[l].N; j++)
      for (int i = 0; i < g[l].N; i++){
	minv = min(revertU(g[l].U[i][j],l), minv);
	maxv = max(abs(revertU(g[l].U[i][j],l)), maxv);
      }
    cerr << "= Level " << l << " =" << endl;
    cerr << "minu: " << minv << endl;
    cerr << "maxu: " << maxv << endl;
    int newUShift = 30 - log(maxv)/log(2);
    cerr << "largest U is " << maxv << " use " << newUShift << " shift" << endl;
    
    minv = 1e8;
    maxv = 1;
    for (int j = 0; j < g[l].N; j++)
      for (int i = 0; i < g[l].N; i++){
	minv = min(revertF(g[l].f[i][j],l), minv);
	maxv = max(abs(revertF(g[l].f[i][j],l)), maxv);
      }
    cerr << "minf: " << minv << endl;
    cerr << "maxf: " << maxv << endl;
    int newUShift2 = 30 - log(maxv/(g[l].h * g[l].h))/log(2);
    cerr << "largest F is " << maxv << " use " << newUShift2 << " shift" << endl;

    bestShift = min(bestShift, min(newUShift2, newUShift));

    double maxf = maxv;
    minv = 1e8;
    maxv = 1;
    for (int j = 0; j < g[l].N; j++)
      for (int i = 0; i < g[l].N; i++){
	minv = min(revertF(g[l].Res[i][j],l), minv);
	maxv = max(abs(revertF(g[l].Res[i][j],l)), maxv);
      }
    cerr << "minres: " << minv << endl;
    cerr << "maxres: " << maxv << endl;
    maxv = max(maxv, maxf);
    int newFShift = 30 - log(maxv)/log(2);
    cerr << "largest F/Res is " << maxv << " use " << newFShift << " F shift" << endl;
  }

  int newUShift = bestShift;
  
  for (int l = 0; l < grids; l++){
    //newUShift = 27;
    lowerUShift = newUShift;
    int extraShift = min(1,max(0, newUShift - shifts[l]));
    cerr << "new U shift is " << shifts[l] + extraShift << " change: " << extraShift << endl;
    shifts[l] = shifts[l] + extraShift;
    fReduction[l] += extraShift;
    for (int j = 0; j < g[l].N; j++)
      {
	for (int i = 0; i < g[l].N; i++){
	  g[l].U[i][j] <<= extraShift;
	}
      }/*
	 newFShift = min(newFShift, lowerFShift);
	 int oldFShift = shifts[l] - fReduction[l];
	 extraShift = newFShift - oldFShift;
	 lowerFShift = newFShift;
	 extraShift = 0;
	 fReduction[l] -= extraShift;

	 cerr << "new F shift is " << newFShift << " ( " << fReduction[l] << ")" << endl;    
	 for (int j = 0; j < g[l].N; j++)
	 {
	 for (int i = 0; i < g[l].N; i++){ //wtf really needed?
	 g[l].f[i][j] <<= extraShift;
	 g[l].Res[i][j] <<= extraShift;
	 }
	 }
       */
  }
  printMinMax(g);
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

bool highPrecMode = 0;

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
      
      g[i].dh = 1.0/(g[i].N-1);

      g[i].h = (g[i].N-1);

      //make arrays contiguous in memory
      g[i].U   = alloc2D<numeric>(g[i].N, g[i].N);
      g[i].Un  = alloc2D<numeric>(g[i].N, g[i].N);
      g[i].Res = alloc2D<numericF>(g[i].N, g[i].N);
      g[i].f   = alloc2D<numericF>(g[i].N, g[i].N);
      g[i].dU   = alloc2D<double>(g[i].N, g[i].N);
      g[i].dUn  = alloc2D<double>(g[i].N, g[i].N);
      g[i].dRes = alloc2D<double>(g[i].N, g[i].N);
      g[i].df   = alloc2D<double>(g[i].N, g[i].N);
    }
  shifts = (int*)calloc(sizeof(int), s.gridCount);
  fReduction = (int*)calloc(sizeof(int), s.gridCount);

  int shift = getInitialShift(g,s);
  for (int i = 0; i < s.gridCount; i++){
    shifts[i] = shift;
    fReduction[i]= 0;
  }
  
  // Setting up problem.
  for (int i = 0; i < s.N; i++) for (int j = 0; j < s.N; j++) g[0].U[i][j] = convertU(s.getInitial(i,j),0);
  for (int i = 0; i < s.N; i++) for (int j = 0; j < s.N; j++) g[0].f[i][j] = convertF(s.getRHS(i,j),0);
  for (int i = 0; i < s.N; i++) for (int j = 0; j < s.N; j++) g[0].df[i][j] = (s.getRHS(i,j),0);

  //adjustShift(g, s.gridCount);
  int blub = 0;
  double minDiff = s.L2NormDiff;
  double prevDiff = s.L2NormDiff;
  int divergesCount = 0;
  
  while (s.L2NormDiff > 1e-2)  // Multigrid solver start
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
      

      for (int grid = s.gridCount-1; grid > 0; grid--) // Going up the V-Cycle
	{
	  s.applyProlongation_(g, grid); // Prolonging solution for coarser level up to finer level
	  s.applyJacobi_(g, grid, s.upRelaxations); // Smoothing finer level
	}
      if (blub++ % 50 == 0){
	cout << "diff: " << s.L2NormDiff << endl;

	adjustShift(g, s.gridCount);
      }
      //       printMinMax(g);
      s.calculateL2Norm_(g, 0); // Calculating Residual L2 Norm
    }  // Multigrid solver end

  cerr << "The last norm difference was " << s.L2NormDiff << endl;
  // Saving solution before returning
  for (int i = 0; i < g[0].N; i++) for (int j = 0; j < g[0].N; j++) s.saveSolution(i, j, revertU(g[0].U[i][j], 0));
}


void applyJacobi(GridLevel* g, int l, int relaxations)
{
  if (highPrecMode){

  }
  else{
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
}

void calculateResidual(GridLevel* g, int l)
{
  
  if (highPrecMode){

  }
  else{
    int N = g[l].N-1;
    numeric fac =  (g[l].h * g[l].h);
    numeric logn = logNum(fac) - fReduction[l];
    for (int x = 1; x < N; x += blockWidth)
      for (int y = 1; y < N; y += blockHeight)
	for (int i = x; i < min(N, x + blockWidth); i++)
	  for (int j = y; j < min (N, y + blockHeight); j++) 
	    g[l].Res[i][j] = g[l].f[i][j] + (((numericF)(g[l].U[i-1][j] + g[l].U[i+1][j] - (((numericF)g[l].U[i][j]) << 2) + g[l].U[i][j-1] + g[l].U[i][j+1])) << (logn) );
  }
}

double calculateL2Norm(GridLevel* g, int l)
{
  if (highPrecMode){

  }
  else{
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
}

void applyRestriction(GridLevel* g, int l)
{
  
  if (highPrecMode){

  }
  else{
    int shift = (shifts[l] - fReduction[l]) - (shifts[l-1] - fReduction[l-1]);
    int N = g[l].N-1;
    for (int x = 1; x < N; x += blockWidth)
      for (int y = 1; y < N; y += blockHeight)
	for (int i = x; i < min(N, x + blockWidth); i++)
	  for (int j = y; j < min (N, y + blockHeight); j++) 
	    g[l].f[i][j] = (( ( ((g[l-1].Res[(i<<1)-1][(j<<1)-1] + g[l-1].Res[(i<<1)-1][(j<<1)+1])>> 1) + ((g[l-1].Res[(i<<1)+1][(j<<1)-1]   + g[l-1].Res[(i<<1)+1][(j<<1)+1]) >> 1) ) >> 3) +
			    (((( g[l-1].Res[(i<<1)-1][(j<<1)]   + g[l-1].Res[(i<<1)][(j<<1)-1]) >> 1)   + ((g[l-1].Res[(i<<1)+1][(j<<1)]     + g[l-1].Res[(i<<1)][(j<<1)+1])>> 1)   ) >> 2) +
			    (( g[l-1].Res[(i<<1)][(j<<1)] )>>2) ) << shift;

    for (int i = 0; i < g[l].N; i++)//wtf
      for (int j = 0; j < g[l].N; j++)
	g[l].U[i][j] = 0;
  }
}

void applyProlongation(GridLevel* g, int l)
{
  
  if (highPrecMode){

  }
  else{
    int shift = shifts[l] - shifts[l-1];
    for (int j = 1; j < g[l].N-1; j++)
      for (int i = 1; i < g[l].N-1; i++)
	g[l-1].Un[(i<<1)][(j<<1)] = g[l].U[i][j] >> shift;

    shift++;
    for (int j = 1; j < g[l].N-1; j++)
      for (int i = 1; i < g[l].N; i++)
	g[l-1].Un[(i<<1)-1][(j<<1)] = ( g[l].U[i-1][j] + g[l].U[i][j] ) >> shift;

    for (int j = 1; j < g[l].N; j++)
      for (int i = 1; i < g[l].N-1; i++)
	g[l-1].Un[(i<<1)][(j<<1)-1] = ( g[l].U[i][j-1] + g[l].U[i][j] ) >> shift;

    for (int j = 1; j < g[l].N; j++)
      for (int i = 1; i < g[l].N; i++)//merge loops
	g[l-1].Un[(i<<1)-1][(j<<1)-1] = ((( g[l].U[i-1][j-1] + g[l].U[i-1][j] ) >> 1 )+ ((g[l].U[i][j-1] + g[l].U[i][j] ) >> 1 )) >> shift;

    for (int j = 0; j < g[l-1].N; j++)//wtf? remove this
      for (int i = 0; i < g[l-1].N; i++)
	g[l-1].U[i][j] += g[l-1].Un[i][j];
  }
}
