#pragma vector aligned
#pragma ivdep

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "auxiliar/auxiliar.hpp"

#include <iostream> //debugging
using namespace std;
using loop_t=unsigned int;
using numeric=long ;
using numericF= long;
#define debug 0
#define SLOW_MODE 1 // set this to one if the code isn't accurate enough

class GridLevel
{
public:
  loop_t N; // Number of points per dimension in the grid level
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

loop_t blockWidth = 32;
loop_t blockHeight = 8;

//the following values represent the decimal places of the fixed point numbers
int* shifts;
int* fReduction;
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
  //if (debug) cerr << in << " -> " << ret << " ";
  return ret;
}

double revertF(numericF in, int l){
  return in / ((double) (1LL << (shifts[l] - fReduction[l])));
}
numericF convertF(double in, int l){
  numericF ret = ((numericF)(in * (1LL << (shifts[l] - fReduction[l]))));
  //if (debug) cerr << in << " -> " << ret << " ";
  return ret;
}

double revertU(numeric in, int l){
  return in / ((double) (1LL << shifts[l]));
}
void printMinMax(GridLevel * g){
  for (loop_t l = 0; l < 6; l++){
    double minv = 1e8;
    double maxv = -1e8;

 
    for (loop_t j = 0; j < g[l].N; j++)
      for (loop_t i = 0; i < g[l].N; i++){
	minv = min(revertU(g[l].U[i][j], l), minv);
	maxv = max(revertU(g[l].U[i][j], l), maxv);
      }
    if (debug) cerr << "= Level " << l << " =" << endl;
    if (debug) cerr << "minu: " << minv << endl;
    if (debug) cerr << "maxu: " << maxv << endl;
    for (loop_t j = 0; j < g[l].N; j++)
      for (loop_t i = 0; i < g[l].N; i++){
	minv = min(revertF(g[l].f[i][j],l), minv);
	maxv = max(revertF(g[l].f[i][j],l), maxv);
      }
    if (debug) cerr << "minf: " << minv << endl;
    if (debug) cerr << "maxf: " << maxv << endl;
  }
}

int getInitialShift(GridLevel * g, Heat2DSetup& s){

  int l = 0;
  double maxv = 1;
  for (loop_t j = 0; j < g[l].N; j++)
    for (loop_t i = 0; i < g[l].N; i++){
      maxv = max(abs(s.getRHS(i,j)), maxv);
    }
  int newFShift = 30 - log(maxv)/log(2);
  if (debug) cerr << "largest F/Res is " << maxv << " use " << newFShift << " F shift" << endl;
  return newFShift;
}

void adjustShift(GridLevel * g, int grids){
  int bestShift = 100;
  for (int l = 0; l < grids; l++){
    double minv = 1e8;
    double maxv = 1;

 
    for (loop_t j = 0; j < g[l].N; j++)
      for (loop_t i = 0; i < g[l].N; i++){
	minv = min(revertU(g[l].U[i][j],l), minv);
	maxv = max(abs(revertU(g[l].U[i][j],l)), maxv);
      }
    if (debug) cerr << "= Level " << l << " =" << endl;
    if (debug) cerr << "minu: " << minv << endl;
    if (debug) cerr << "maxu: " << maxv << endl;
    int newUShift = 30 - log(maxv)/log(2);
    if (debug) cerr << "largest U is " << maxv << " use " << newUShift << " shift" << endl;
    
    minv = 1e8;
    maxv = 1;
    for (loop_t j = 0; j < g[l].N; j++)
      for (loop_t i = 0; i < g[l].N; i++){
	minv = min(revertF(g[l].f[i][j],l), minv);
	maxv = max(abs(revertF(g[l].f[i][j],l)), maxv);
      }
    if (debug) cerr << "minf: " << minv << endl;
    if (debug) cerr << "maxf: " << maxv << endl;
    int newUShift2 = 30 - log(maxv/(g[l].h * g[l].h))/log(2);
    if (debug) cerr << "largest F is " << maxv << " use " << newUShift2 << " shift" << endl;

    bestShift = min(bestShift, min(newUShift2, newUShift));

    double maxf = maxv;
    minv = 1e8;
    maxv = 1;
    for (loop_t j = 0; j < g[l].N; j++)
      for (loop_t i = 0; i < g[l].N; i++){
	minv = min(revertF(g[l].Res[i][j],l), minv);
	maxv = max(abs(revertF(g[l].Res[i][j],l)), maxv);
      }
    if (debug) cerr << "minres: " << minv << endl;
    if (debug) cerr << "maxres: " << maxv << endl;
    maxv = max(maxv, maxf);
    int newFShift = 30 - log(maxv)/log(2);
    if (debug) cerr << "largest F/Res is " << maxv << " use " << newFShift << " F shift" << endl;
  }

  int newUShift = bestShift;
  
  for (int l = 0; l < grids; l++){
    //newUShift = 27;
    int extraShift = min(1,max(0, newUShift - shifts[l]));
    if (debug) cerr << "new U shift is " << shifts[l] + extraShift << " change: " << extraShift << endl;
    shifts[l] = shifts[l] + extraShift;
    fReduction[l] += extraShift;
    for (loop_t j = 0; j < g[l].N; j++)
      {
	for (loop_t i = 0; i < g[l].N; i++){
	  g[l].U[i][j] <<= extraShift;
	}
      }
  }
  printMinMax(g);
}


//helper to make allocations nicer
template<typename T>
T** alloc2D(loop_t dim1, loop_t dim2){
  while (dim1%blockWidth != 0) dim1++;
  while (dim2%blockHeight != 0) dim2++;
  T** ret = (T**)malloc(sizeof(T*) * dim1 + sizeof(T) * dim1 * dim2);
  T* start = (T*)(ret + dim1);
  for (loop_t i = 0; i < dim1; i++) ret[i] = start + i * dim2;
  for (loop_t i = 0; i < dim1; i++) for (loop_t j = 0; j < dim2; j++) ret[i][j] = 0;
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
  for (loop_t i = 0; i < s.N; i++) for (loop_t j = 0; j < s.N; j++) g[0].U[i][j] = convertU(s.getInitial(i,j),0);
  for (loop_t i = 0; i < s.N; i++) for (loop_t j = 0; j < s.N; j++) g[0].f[i][j] = convertF(s.getRHS(i,j),0);
  for (loop_t i = 0; i < s.N; i++) for (loop_t j = 0; j < s.N; j++) g[0].df[i][j] = (s.getRHS(i,j));
  for (loop_t i = 0; i < s.N; i++) for (loop_t j = 0; j < s.N; j++) g[0].dU[i][j] = s.getInitial(i,j);

  highPrecMode = 0;
  //adjustShift(g, s.gridCount);

  for (loop_t mode = 0; mode <(1 + SLOW_MODE) ; mode++){
    
    int blub = 0;
    //s.calculateL2Norm_(g, 0); // Calculating Residual L2 Norm
    highPrecMode = mode;
    double minDiff = s.L2NormDiff;
    double prevDiff = s.L2NormDiff;
    int divergesCount = 0;

    if (mode == 1){
      highPrecMode = 1;
      s.calculateResidual_(g, 0); // Calculating Initial Residual
      s.calculateL2Norm_(g, 0); // Calculating Residual L2 Norm
      cout << "init diff: " << s.L2Norm << endl;
    }
    s.L2NormDiff = 1e100;
    s.L2Norm = 1e100;
    
    while (s.L2NormDiff > (highPrecMode ? s.tolerance : 1e-3))  // Multigrid solver start
      {
	if (s.L2Norm < minDiff){
	  minDiff =s.L2Norm;
	  //cerr << "The new best norm is " << minDiff << endl;
	  if (mode == 0)for (loop_t i = 0; i < s.N; i++) for (loop_t j = 0; j < s.N; j++) g[0].dU[i][j] = revertU(g[0].U[i][j],0);
	  if (mode == 0)for (loop_t i = 0; i < s.N; i++) for (loop_t j = 0; j < s.N; j++) g[0].dRes[i][j] = revertF(g[0].Res[i][j],0);
	} 
	if (s.L2NormDiff > 2*minDiff && s.L2Norm > prevDiff){
	  divergesCount ++;//diverging, numeric limit reached
	}else if (s.L2NormDiff < 1e-7 + minDiff){
	  divergesCount = 0;
	}
	else
	  {
	    divergesCount = max(0, divergesCount - 1);
	  }
	
	prevDiff = s.L2Norm;
	if (divergesCount > 10){
	  cerr << "Solver is no longer converging. The best norm difference was " << minDiff << endl;
	  break;
	}
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
	
	//       printMinMax(g);
	s.calculateResidual_(g, 0); // Calculating Initial Residual
	s.calculateL2Norm_(g, 0); // Calculating Residual L2 Norm
	if (blub++ % 5 == 0 || 0){
	  //cout << "diff: " << s.L2NormDiff << endl;
	  //cout << "norm: " << s.L2Norm << endl;
	  //cout << "mindiff: " << minDiff << endl;

	  if (mode == 0)adjustShift(g, s.gridCount);
	}
      }  // Multigrid solver end
    
    cerr << "The last norm difference was " << s.L2NormDiff << " (norm "  << s.L2Norm<< ") with " << shifts[0] << " shifts" << endl;
  }
  // Saving solution before returning
  for (loop_t i = 0; i < g[0].N; i++) for (loop_t j = 0; j < g[0].N; j++) s.saveSolution(i, j, g[0].dU[i][j]);
}


void applyJacobi(GridLevel* g, int l, int relaxations)
{
  if (highPrecMode){
    double fac = pow(g[l].dh,2)/4;
    loop_t N = g[l].N-1;
    if (l > 0){
      for (loop_t i = 1; i < N; i++)
	for (loop_t j = 1; j < N; j++) // Perform a Jacobi Iteration
	  g[l].dU[i][j] = g[l].df[i][j]*fac;
    }
    for (int r = (l > 0); r < relaxations; r++) {
	swap(g[l].dUn, g[l].dU);

    for (loop_t x = 1; x < N; x += blockWidth)
      for (loop_t y = 1; y < N; y += blockHeight){
	loop_t W = min(N, x + blockWidth);
	loop_t H = min(N, y + blockHeight);
	for (loop_t i = x; i < W; i++){
	  for (loop_t j = y; j < H; j++) {
		g[l].dU[i][j] = (g[l].dUn[i-1][j] + g[l].dUn[i+1][j] + g[l].dUn[i][j-1] + g[l].dUn[i][j+1])* 0.25 + g[l].df[i][j]*fac;
	  }
	}
      }
      }
  }
  else{
    if (l > 0){
      loop_t N = g[l].N-1;
      numeric fac =  (g[l].h * g[l].h);
      numeric logn = logNum(fac) - fReduction[l] + 2;
      for (loop_t i = 1; i < N; i++)
	for (loop_t j = 1; j < N; j++) // Perform a Jacobi Iteration
	  g[l].U[i][j] =  (numeric)(g[l].f[i][j] >> logn);
    }
    for (int r = (l > 0); r < relaxations; r++) {
      
	swap(g[l].Un, g[l].U);

	loop_t N = g[l].N-1;
	numeric fac =  (g[l].h * g[l].h);

	numeric logn = logNum(fac) - fReduction[l] + 1;
    for (loop_t x = 1; x < N; x += blockWidth)
      for (loop_t y = 1; y < N; y += blockHeight){
	loop_t W = min(N, x + blockWidth);
	loop_t H = min(N, y + blockHeight);
	for (loop_t i = x; i < W; i++){
	  for (loop_t j = y; j < H; j++) {
		g[l].U[i][j] = ((((g[l].Un[i-1][j] + g[l].Un[i+1][j]) >> 1) + ((g[l].Un[i][j-1] + g[l].Un[i][j+1]) >> 1)) + (numeric)(g[l].f[i][j] >> logn)) >> 1;
	  }
	}
      }
    }
  }
    }

void calculateResidual(GridLevel* g, int l)
{
  
  if (highPrecMode){
    loop_t N = g[l].N-1;
    double fac = 1 / (pow(g[l].dh,2));
    for (loop_t x = 1; x < N; x += blockWidth)
      for (loop_t y = 1; y < N; y += blockHeight){
	loop_t W = min(N, x + blockWidth);
	loop_t H = min(N, y + blockHeight);
	for (loop_t i = x; i < W; i++){
	  for (loop_t j = y; j < H; j++) {
	    g[l].dRes[i][j] = g[l].df[i][j] + (g[l].dU[i-1][j] + g[l].dU[i+1][j] - 4*g[l].dU[i][j] + g[l].dU[i][j-1] + g[l].dU[i][j+1]) * fac;
	  }
	}
      }
  }
  else{
    loop_t N = g[l].N-1;
    numeric fac =  (g[l].h * g[l].h);
    numeric logn = logNum(fac) - fReduction[l];
    for (loop_t x = 1; x < N; x += blockWidth)
      for (loop_t y = 1; y < N; y += blockHeight){
	loop_t W = min(N, x + blockWidth);
	loop_t H = min(N, y + blockHeight);
	for (loop_t i = x; i < W; i++){
	  for (loop_t j = y; j < H; j++) {
	    g[l].Res[i][j] = g[l].f[i][j] + (((numericF)(g[l].U[i-1][j] + g[l].U[i+1][j] - (((numericF)g[l].U[i][j]) << 2) + g[l].U[i][j-1] + g[l].U[i][j+1])) << (logn) );
	  }
	}
      }
  }
}

double calculateL2Norm(GridLevel* g, int l)
{
  if (highPrecMode){
    double tmp = 0.0;
    for (loop_t i = 0; i < g[l].N; i++)
      for (loop_t j = 0; j < g[l].N; j++)
	tmp += g[l].dRes[i][j]*g[l].dRes[i][j];

    return sqrt(tmp);
  }
  else{
    int s = l2NormShift;
    long long tmp = 0.0;
    for (loop_t i = 0; i < g[l].N; i++)
      for (loop_t j = 0; j < g[l].N; j++){
	tmp += intSquare(g[l].Res[i][j] << s, l);
	//if (debug) cerr << g[l].Res[i][j] << " | " << intSquare(g[l].Res[i][j]) << " | " << tmp << endl;
      }

    return sqrt((tmp) / (double)(1 << (s << 1)));
  }
}

void applyRestriction(GridLevel* g, int l)
{
  
  if (highPrecMode){
    loop_t N = g[l].N-1;
    for (loop_t x = 1; x < N; x += blockWidth)
      for (loop_t y = 1; y < N; y += blockHeight){
	loop_t W = min(N, x + blockWidth);
	loop_t H = min(N, y + blockHeight);
	for (loop_t i = x; i < W; i++){
	  for (loop_t j = y; j < H; j++) {
	    g[l].df[i][j] = ( 1.0*( g[l-1].dRes[2*i-1][2*j-1] + g[l-1].dRes[2*i-1][2*j+1] + g[l-1].dRes[2*i+1][2*j-1]   + g[l-1].dRes[2*i+1][2*j+1] ) +
			      2.0*( g[l-1].dRes[2*i-1][2*j]   + g[l-1].dRes[2*i][2*j-1]   + g[l-1].dRes[2*i+1][2*j]     + g[l-1].dRes[2*i][2*j+1]   ) +
			      4.0*( g[l-1].dRes[2*i][2*j] ) ) * 0.0625;
	  }
	}
      }

  }
  else{
    int shift = (shifts[l] - fReduction[l]) - (shifts[l-1] - fReduction[l-1]);
    loop_t N = g[l].N-1;
    for (loop_t x = 1; x < N; x += blockWidth)
      for (loop_t y = 1; y < N; y += blockHeight){
	loop_t W = min(N, x + blockWidth);
	loop_t H = min(N, y + blockHeight);
	for (loop_t i = x; i < W; i++){
	  for (loop_t j = y; j < H; j++) {
	    g[l].f[i][j] = (( ( ((g[l-1].Res[(i<<1)-1][(j<<1)-1] + g[l-1].Res[(i<<1)-1][(j<<1)+1])>> 1) + ((g[l-1].Res[(i<<1)+1][(j<<1)-1]   + g[l-1].Res[(i<<1)+1][(j<<1)+1]) >> 1) ) >> 3) +
			    (((( g[l-1].Res[(i<<1)-1][(j<<1)]   + g[l-1].Res[(i<<1)][(j<<1)-1]) >> 1)   + ((g[l-1].Res[(i<<1)+1][(j<<1)]     + g[l-1].Res[(i<<1)][(j<<1)+1])>> 1)   ) >> 2) +
			    (( g[l-1].Res[(i<<1)][(j<<1)] )>>2) ) << shift;
	  }
	}
      }

  }
}

void applyProlongation(GridLevel* g, int l)
{
  
  if (highPrecMode){
    
    loop_t i;
    loop_t N = g[l].N-1;

    for (loop_t x = 1; x < N; x += blockWidth)
      for (loop_t y = 1; y < N; y += blockHeight){
	loop_t W = min(N, x + blockWidth);
	loop_t H = min(N, y + blockHeight);
	for (loop_t i = x; i < W; i++){
	  for (loop_t j = y; j < H; j++) {
	g[l-1].dU[2*i][2*j] += g[l].dU[i][j];
	g[l-1].dU[2*i-1][2*j] += ( g[l].dU[i-1][j] + g[l].dU[i][j] ) * 0.5;
	g[l-1].dU[2*i][2*j-1] += ( g[l].dU[i][j-1] + g[l].dU[i][j] )  * 0.5;
	g[l-1].dU[2*i-1][2*j-1] += ( g[l].dU[i-1][j-1] + g[l].dU[i-1][j] + g[l].dU[i][j-1] + g[l].dU[i][j] ) * 0.25;
      }
    }
	  }

    for (loop_t j = 1; j < g[l].N-1; j++){
      i = g[l].N -1;
      g[l-1].dU[2*i-1][2*j] += ( g[l].dU[i-1][j] + g[l].dU[i][j] ) * 0.5;
      g[l-1].dU[2*i-1][2*j-1] += ( g[l].dU[i-1][j-1] + g[l].dU[i-1][j] + g[l].dU[i][j-1] + g[l].dU[i][j] ) * 0.25;
    }
    
    int j =g[l].N -1 ;
    for (i = 1; i < N; i++){
      g[l-1].dU[2*i][2*j-1] += ( g[l].dU[i][j-1] + g[l].dU[i][j] )  * 0.5;
      g[l-1].dU[2*i-1][2*j-1] += ( g[l].dU[i-1][j-1] + g[l].dU[i-1][j] + g[l].dU[i][j-1] + g[l].dU[i][j] ) * 0.25;
      
    }
    i = g[l].N -1;
    g[l-1].dU[2*i-1][2*j-1] += ( g[l].dU[i-1][j-1] + g[l].dU[i-1][j] + g[l].dU[i][j-1] + g[l].dU[i][j] ) * 0.25;

  }
  else{
    int shift = shifts[l] - shifts[l-1];
    int shift2 = shift + 1;
    loop_t i;
    loop_t N = g[l].N-1;

    for (loop_t x = 1; x < N; x += blockWidth)
      for (loop_t y = 1; y < N; y += blockHeight){
	loop_t W = min(N, x + blockWidth);
	loop_t H = min(N, y + blockHeight);
	for (loop_t i = x; i < W; i++){
	  for (loop_t j = y; j < H; j++) {
	g[l-1].U[(i<<1)][(j<<1)] += g[l].U[i][j] >> shift;
	g[l-1].U[(i<<1)-1][(j<<1)] += ( g[l].U[i-1][j] + g[l].U[i][j] ) >> shift2;
	g[l-1].U[(i<<1)-1][(j<<1)-1] += ((( g[l].U[i-1][j-1] + g[l].U[i-1][j] ) >> 1 )+ ((g[l].U[i][j-1] + g[l].U[i][j] ) >> 1 )) >> shift2;
	g[l-1].U[(i<<1)][(j<<1)-1] += ( g[l].U[i][j-1] + g[l].U[i][j] ) >> shift2;
      }
    }
      }
    
    for (loop_t j = 1; j < N; j++) {
      i = g[l].N -1;
      g[l-1].U[(i<<1)-1][(j<<1)] += ( g[l].U[i-1][j] + g[l].U[i][j] ) >> shift2;
      g[l-1].U[(i<<1)-1][(j<<1)-1] += ((( g[l].U[i-1][j-1] + g[l].U[i-1][j] ) >> 1 )+ ((g[l].U[i][j-1] + g[l].U[i][j] ) >> 1 )) >> shift2;
    }
    int j =g[l].N -1 ;
    for (i = 1; i < N; i++){
      g[l-1].U[(i<<1)][(j<<1)-1] += ( g[l].U[i][j-1] + g[l].U[i][j] ) >> shift2;
      g[l-1].U[(i<<1)-1][(j<<1)-1] += ((( g[l].U[i-1][j-1] + g[l].U[i-1][j] ) >> 1 )+ ((g[l].U[i][j-1] + g[l].U[i][j] ) >> 1 )) >> shift2;
    }
    i = g[l].N -1;
    g[l-1].U[(i<<1)-1][(j<<1)-1] += ((( g[l].U[i-1][j-1] + g[l].U[i-1][j] ) >> 1 )+ ((g[l].U[i][j-1] + g[l].U[i][j] ) >> 1 )) >> shift2;

  }
}
