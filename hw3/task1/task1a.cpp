#include <iostream>
#include "model/grass.hpp"
#include "korali.h"

#define MM 80.0
#define PH 6.0

// Objective function
double maximize_weed(double* x)
{
  //std::cout << "@ " << x[0] << " | " << x[1] << ": " << getGrassHeight(x[0], x[1], PH, MM) << std::endl;
  return getGrassHeight(x[0], x[1], PH, MM); //why is the grass height negative?
}

int main(int argc, char* argv[])
{
  // We want to maximize the objective function directly.
  auto problem = Korali::Problem::Direct(maximize_weed);

  Korali::Parameter::Uniform x("x(km)", 0.0, 5.0);
  Korali::Parameter::Uniform y("y(km)", 0.0, 5.0);
  problem.addParameter(&x);
  problem.addParameter(&y);

  // Use CMAES to find the maximum of weed
  auto maximizer = Korali::Solver::CMAES(&problem);

  // CMAES-specific configuration.
  // StopMinDeltaX defines how close we want the result to be
  // from the actual global minimum. The smaller this value is, 
  // the more generations it may take to find it.
  maximizer.setStopMinDeltaX(1e-11);

  // Population size defines how many samples per generations we want to run
  // For CMAES, a small number of samples (64-256) will do the trick.
  maximizer.setPopulationSize(128);

  // Run CMAES and report the result
  maximizer.run();

  return 0;
}
