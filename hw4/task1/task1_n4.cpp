#include "model/heat2d.hpp"
#include "korali.h"
#include <upcxx/upcxx.hpp>

int main(int argc, char* argv[])
{

  upcxx::init();
  // Loading temperature measurement data
  FILE* dataFile = fopen("data_n4.in", "r");
  fscanf(dataFile, "%lu", &p.nPoints);

  p.xPos    = (double*) calloc (sizeof(double), p.nPoints);
  p.yPos    = (double*) calloc (sizeof(double), p.nPoints);
  p.refTemp = (double*) calloc (sizeof(double), p.nPoints);

  for (int i = 0; i < p.nPoints; i++)
    {
      fscanf(dataFile, "%le ", &p.xPos[i]);
      fscanf(dataFile, "%le ", &p.yPos[i]);
      fscanf(dataFile, "%le ", &p.refTemp[i]);
    }

  auto problem = Korali::Problem::Posterior(heat2DSolver);
  // 4-Candle Model x4 parameters/candle 
  p.nCandles = 4;
        
  /*
    TODO : fill in all parameter information required 
  */
  Korali::Parameter::Uniform c1x("Candle 1 X", 0.0, 0.5);
  Korali::Parameter::Uniform c1y("Candle 1 Y", 0.0, 1.0);
  Korali::Parameter::Uniform c1w("Candle 1 Width", 0.04, 0.06);
  Korali::Parameter::Uniform c1i("Candle 1 Intensity", 0.4, 0.6);
  Korali::Parameter::Uniform c2x("Candle 2 X", 0.0, 0.5);
  Korali::Parameter::Uniform c2y("Candle 2 Y", 0.0, 1.0);
  Korali::Parameter::Uniform c2w("Candle 2 Width", 0.04, 0.06);
  Korali::Parameter::Uniform c2i("Candle 2 Intensity", 0.4, 0.6);
  Korali::Parameter::Uniform c3x("Candle 3 X", 0.5, 1.0);
  Korali::Parameter::Uniform c3y("Candle 3 Y", 0.0, 1.0);
  Korali::Parameter::Uniform c3w("Candle 3 Width", 0.04, 0.06);
  Korali::Parameter::Uniform c3i("Candle 3 Intensity", 0.4, 0.6);
  Korali::Parameter::Uniform c4x("Candle 4 X", 0.5, 1.0);
  Korali::Parameter::Uniform c4y("Candle 4 Y", 0.0, 1.0);
  Korali::Parameter::Uniform c4w("Candle 4 Width", 0.04, 0.06);
  Korali::Parameter::Uniform c4i("Candle 4 Intensity", 0.4, 0.6);

  problem.addParameter(&c1x);
  problem.addParameter(&c1y);
  problem.addParameter(&c1i);
  problem.addParameter(&c1w);
  problem.addParameter(&c2x);
  problem.addParameter(&c2y);
  problem.addParameter(&c2i);
  problem.addParameter(&c2w);
  problem.addParameter(&c3x);
  problem.addParameter(&c3y);
  problem.addParameter(&c3i);
  problem.addParameter(&c3w);
  problem.addParameter(&c4x);
  problem.addParameter(&c4y);
  problem.addParameter(&c4i);
  problem.addParameter(&c4w);

  problem.setReferenceData(p.nPoints, p.refTemp);


  auto solver = Korali::Solver::CMAES(&problem);

  int Ng = 2000; // max generations for CMAES
	
  solver.setStopMinDeltaX(1e-6);
  solver.setPopulationSize(23); // ~4+3*log(N) //changed for subtask 4
  solver.setMu(4);
  solver.setMaxGenerations(Ng);
  solver.run();

  upcxx::finalize();

  return 0;
}
