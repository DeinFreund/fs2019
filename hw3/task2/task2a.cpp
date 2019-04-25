#include "model/heat2d.hpp"
#include "korali.h"
#include <string>
#include <iostream>

void likelihood(double * pars, double * eval){

  heat2DSolver(pars, eval);

  
}

int main(int argc, char* argv[])
{
  // Loading temperature measurement data

  auto problem = Korali::Problem::Likelihood(likelihood);


  FILE* dataFile = fopen("data.in", "r");
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


  p.nCandles = 3;
  Korali::Parameter::Uniform c1x("Candle 1 X", 0.0, 0.5 + 0.5 * (p.nCandles == 1));
  Korali::Parameter::Uniform c1y("Candle 1 Y", 0.0, 1.0);
  Korali::Parameter::Uniform c2x("Candle 2 X", 0.5, 1.0);
  Korali::Parameter::Uniform c2y("Candle 2 Y", 0.0, 1.0);
  Korali::Parameter::Uniform c3x("Candle 3 X", 0.5, 1.0);
  Korali::Parameter::Uniform c3y("Candle 3 Y", 0.0, 1.0);


  if (p.nCandles > 0){
    problem.addParameter(&c1x);
    problem.addParameter(&c1y);
  }
  if (p.nCandles > 1){
    problem.addParameter(&c2x);
    problem.addParameter(&c2y);
  }
  if (p.nCandles > 2){
    problem.addParameter(&c3x);
    problem.addParameter(&c3y);
  }
  // Very important: dont forget to give the reference data to Korali!
  problem.setReferenceData(p.nPoints, p.refTemp);

  auto sampler = Korali::Solver::TMCMC(&problem);

  // TMCMC-specific configuration.
  // Number of samples to represent the distribution.
  // The more samples, the more precise the representation will be
  // but may take more time to run per generation, and more generations
  // to find a perfectly annealing representation.
  sampler.setPopulationSize(10000); //very high to be accurate

  // Defines the 'sensitivity' of re-sampling. That is, how much the new
  // samples within a chain will scale during evaluation. A higher value
  // is better to explore a larger space, while a lower value will be
  // more precise for small magnitude parameters.
  sampler.setCovarianceScaling(0.2);

  // Run TMCMC to produce tmcmc.txt.
  // Use plotmatrix_hist to see the result of the sampling.
  sampler.run();


  // Use CMAES to find the maximum of -weed(x)
  auto maximizer = Korali::Solver::CMAES(&problem);

  // CMAES-specific configuration.
  // StopMinDeltaX defines how close we want the result to be
  // from the actual global minimum. The smaller this value is, 
  // the more generations it may take to find it.
  maximizer.setStopMinDeltaX(1e-11);

  // Population size defines how many samples per generations we want to run
  // For CMAES, a small number of samples (64-256) will do the trick.
  maximizer.setPopulationSize(128);

  // For this problem, we may need to run more generations
  maximizer.setMaxGenerations(1000);

  // Run CMAES and report the result
  maximizer.run();


  return 0;
}
