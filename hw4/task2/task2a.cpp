#include <stdio.h>
#include <chrono>
#include <upcxx/upcxx.hpp>
#include "sampler/sampler.hpp"

size_t nSamples;
size_t nParameters;

#define NSAMPLES 240
#define NPARAMETERS 2

int main(int argc, char* argv[])
{
  upcxx::init();
  int rankId    = upcxx::rank_me();
  int rankCount = upcxx::rank_n();

  nSamples = NSAMPLES / rankCount;
  if (rankId < NSAMPLES % rankCount) nSamples++;
  nParameters = NPARAMETERS;


  upcxx::global_ptr<double> res_ptr;
  if (rankId == 0) res_ptr = upcxx::new_array<double>(NSAMPLES);
  upcxx::broadcast(&res_ptr, 1, 0).wait();
  
  double* sampleArray = initializeSampler(NSAMPLES, nParameters);
  double* resultsArray = (double*) calloc (nSamples, sizeof(double));

  printf("Rank %1d: Processing %ld Samples each with %ld Parameter(s)...\n", rankId, nSamples, nParameters);


  int offset = (NSAMPLES / rankCount) * rankId + std::min(NSAMPLES % rankCount, rankId);

  auto start = std::chrono::system_clock::now();
  for (size_t i = 0; i < nSamples; i++) resultsArray[i] = evaluateSample(&sampleArray[(i+offset)*nParameters]);
  

  printf("Rank %1d: Storing results at offset %1d\n", rankId, offset);

  upcxx::rput(resultsArray, res_ptr + offset, nSamples).wait();

  upcxx::barrier();

  auto end = std::chrono::system_clock::now();

  if (rankId == 0){
    printf("Checking results...\n");

      double* resultsArray = (double*) calloc (nSamples, sizeof(double));
      for (int i = 0; i < nSamples; i++) resultsArray[i] = *(res_ptr.local()+i);
      for (int i = 0; i < nSamples; i++) printf("Result %1d: %f\n", i, resultsArray[i]);
    checkResults(res_ptr.local());

    double totalTime = std::chrono::duration<double>(end-start).count();
    printf("Total Running Time: %.3fs\n", totalTime);
  }

  upcxx::finalize();
  return 0;
}
