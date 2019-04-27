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
upcxx::global_ptr<double> sampleArray;

int main(int argc, char* argv[])
{
  upcxx::init();
  rankId    = upcxx::rank_me();
  rankCount = upcxx::rank_n();
  nSamples = NSAMPLES;
  nParameters = NPARAMETERS;

  if (rankId == 0) printf("Processing %ld Samples (24 initially available), each with %ld Parameter(s)...\n", nSamples, nParameters);
  if (rankId == 0) initializeSampler(nSamples, nParameters);

  auto t0 = std::chrono::system_clock::now();
  
  
  upcxx::future<> the_life_the_universe_and_everything = upcxx::make_future();
  bool finished = false;
  if (rankId == 0){
    for (int target = 0; target < nSamples; target++){
      auto dest = upcxx::rpc(target % rankCount, [](){return sampleArray = upcxx::new_array<double>(NPARAMETERS);}).wait();
      double sample[nParameters];
      getSample(target, sample);
      upcxx::rput(sample, dest, nParameters).wait();
      auto fut = upcxx::rpc(target % rankCount, [](){return evaluateSample(sampleArray.local());}).then([target](double res){updateEvaluation(target, res);});
      the_life_the_universe_and_everything = upcxx::when_all(the_life_the_universe_and_everything, fut);
    } 
    the_life_the_universe_and_everything.wait();
    for (int target = 0; target < rankCount; target++){
      upcxx::rpc(target, [&finished](){finished = true;}).wait();
    }
  }

  while (!finished) upcxx::progress();


  auto t1 = std::chrono::system_clock::now();

  if (rankId == 0)
    {
      checkResults(/* This will FAIL */ ); // Make sure you check results!
      double evalTime = std::chrono::duration<double>(t1-t0).count();
      printf("Total Running Time: %.3fs\n", evalTime);
    }

  upcxx::finalize();
  return 0;
}


