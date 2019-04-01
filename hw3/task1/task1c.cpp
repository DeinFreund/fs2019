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

  // Use TMCMC to sample the weed
  auto sampler = Korali::Solver::TMCMC(&problem);

  // TMCMC-specific configuration.
  // Number of samples to represent the distribution.
  // The more samples, the more precise the representation will be
  // but may take more time to run per generation, and more generations
  // to find a perfectly annealing representation.
  sampler.setPopulationSize(10000);

  // Defines the 'sensitivity' of re-sampling. That is, how much the new
  // samples within a chain will scale during evaluation. A higher value
  // is better to explore a larger space, while a lower value will be
  // more precise for small magnitude parameters.
  sampler.setCovarianceScaling(0.2);

  // Run TMCMC to produce tmcmc.txt.
  // Use plotmatrix_hist to see the result of the sampling.
  sampler.run();
  return 0;
}
