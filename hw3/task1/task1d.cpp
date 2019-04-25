#include "model/grass.hpp"
#include "korali.h"

// Grass Height at different spots, as measured by Herr Kueheli.
size_t  nSpots;
double* xPos;
double* yPos;
double* heights;

void likelihood_weed(double* x, double* fx)
{
	
  double ph = x[0];
  double mm = x[1];

  // Filling fx for Korali to calculate likelihood
  // No need to use the negative because this is not about maximizing or minimizing
  // the model, but the likelihood. Korali will take care of that.
  for (size_t i = 0; i < nSpots; i++)
    fx[i] = getGrassHeight(xPos[i], yPos[i], ph, mm);
}

int main(int argc, char* argv[])
{
  // Loading grass height data

  FILE* dataFile = fopen("grass.in", "r");

  fscanf(dataFile, "%lu", &nSpots);
  xPos     = (double*) calloc (sizeof(double), nSpots);
  yPos     = (double*) calloc (sizeof(double), nSpots);
  heights  = (double*) calloc (sizeof(double), nSpots);

  for (int i = 0; i < nSpots; i++)
    {
      fscanf(dataFile, "%le ", &xPos[i]);
      fscanf(dataFile, "%le ", &yPos[i]);
      fscanf(dataFile, "%le ", &heights[i]);
    }

  // We want to maximize the posterior distribution of the parameters.
  // We do have some prior information now.
  auto problem = Korali::Problem::Posterior(likelihood_weed);


  Korali::Parameter::Uniform ph("pH", 4.0, 9.0); 

  Korali::Parameter::Gaussian mm("mm", 90.0, 20.0);
  mm.setBounds(10, 170); // 4 sigma bound.

  problem.addParameter(&ph);
  problem.addParameter(&mm);

  // Very important: dont forget to give the reference data to Korali!
  problem.setReferenceData(nSpots, heights);

  // Use CMAES to find the maximum of -weed(x)
  auto maximizer = Korali::Solver::CMAES(&problem);

  // CMAES-specific configuration.
  // StopMinDeltaX defines how close we want the result to be
  // from the actual global minimum. The smaller this value is, 
  // the more generations it may take to find it.
  maximizer.setStopMinDeltaX(1e-11);

  // Population size defines how many samples per generations we want to run
  // For CMAES, a small number of samples (64-256) will do the trick.
  maximizer.setPopulationSize(64);

  // For this problem, we may need to run more generations
  maximizer.setMaxGenerations(1000);

  // Run CMAES and report the result
  maximizer.run();

  // Use TMCMC to sample the weed
  auto sampler = Korali::Solver::TMCMC(&problem);

  // TMCMC-specific configuration.
  // Number of samples to represent the distribution.
  // The more samples, the more precise the representation will be
  // but may take more time to run per generation, and more generations
  // to find a perfectly annealing representation.
  sampler.setPopulationSize(100000); //very high to be accurate

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
