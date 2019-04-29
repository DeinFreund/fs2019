#include <stdio.h>
#include <chrono>
#include <queue>
#include <upcxx/upcxx.hpp>
#include "sampler/sampler.hpp"

#define NSAMPLES 240
#define NPARAMETERS 2

int rankId;
int rankCount;
size_t nSamples;
size_t nParameters;

double* sampleArray;

upcxx::global_ptr<double> res_ptr;

void calcSamples(size_t nSamples, size_t nParameters, int offset){
  double* resultsArray = (double*) calloc (nSamples, sizeof(double));

  printf("Rank %1d: Processing %ld Samples each with %ld Parameter(s)...\n", rankId, nSamples, nParameters);

  for (size_t i = 0; i < nSamples; i++) resultsArray[i] = evaluateSample(&sampleArray[(i+offset)*nParameters]);

  printf("Rank %1d: Storing results at offset %1d\n", rankId, offset);

  upcxx::rput(resultsArray, res_ptr + offset, nSamples).wait();
  free(resultsArray);
}

int processedIndex = 0;
int nodesFinished = 1;
std::pair<int,int> getWork(){
    int batchSize = std::min(NSAMPLES - processedIndex, (NSAMPLES - processedIndex) / 50 + 2);
    printf("Giving out a batch of size %1d, current index is %1d.\n", batchSize, processedIndex);    
    processedIndex += batchSize;
    if (batchSize == 0) nodesFinished ++;
    return {batchSize, processedIndex - batchSize};
}

int main(int argc, char* argv[])
{
  upcxx::init();
  rankId    = upcxx::rank_me();
  rankCount = upcxx::rank_n();
  nSamples = NSAMPLES;
  nParameters = NPARAMETERS;

  if (rankId == 0) printf("Processing %ld Samples each with %ld Parameter(s)...\n", nSamples, nParameters);

  auto t0 = std::chrono::system_clock::now();

  if (rankId == 0) res_ptr = upcxx::new_array<double>(NSAMPLES);
  upcxx::broadcast(&res_ptr, 1, 0).wait();

  sampleArray = initializeSampler(nSamples, nParameters);
  if (rankId > 0){
    //im a consumer
    while(1){
      auto work = upcxx::rpc(0, getWork).wait();
      if (work.first == 0) break;
      calcSamples(work.first, nParameters, work.second);
    }
  }else{
    //im a producer
    while(nodesFinished < rankCount) upcxx::progress();
  }
  
  printf("Rank %1d has finished.\n", rankId);
  auto t1 = std::chrono::system_clock::now();

  upcxx::barrier();
  if (rankId == 0)  t1 = std::chrono::system_clock::now();

  if (rankId == 0)
  {
    printf("Checking results...\n");
    
    checkResults(res_ptr.local()); // Make sure you check results!
    double evalTime = std::chrono::duration<double>(t1-t0).count();
    printf("Total Running Time: %.3fs\n", evalTime);
  }else{
    double evalTime = std::chrono::duration<double>(t1-t0).count();
    printf("Running Time: %.3fs\n", evalTime);

  }
  
  upcxx::finalize();
  return 0;
}
