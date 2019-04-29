#include "conduits/upcxx.h"
#include "solvers/base.h"
#include <queue>

extern Korali::Solver::Base*  _solver;
extern Korali::Problem::Base* _problem;

int rankId;
int rankCount;
size_t nSamples;
size_t nParameters;
upcxx::global_ptr<double> cache;
upcxx::global_ptr<double> sampleArrayPointer;
int curSample = -1;

upcxx::global_ptr<double> sendDataPre(int sampleId){
  curSample = sampleId;
  return cache;
}

upcxx::future<> sendData(int target, int sampleId,  double* data, int len){
  return upcxx::rpc(target, sendDataPre, sampleId).then([data, len](upcxx::global_ptr<double> t){return upcxx::rput(data, t, len);});
}

double evaluateSample(){
  return _problem->evaluateSample(cache.local());
}

Korali::Conduit::UPCXX::UPCXX(Korali::Solver::Base* solver) : Base::Base(solver) {};

upcxx::future<> combinedFutures;

void Korali::Conduit::UPCXX::processSample(size_t sampleId)
{
  int worker = (sampleId % (rankCount - 1)) + 1; //this could use load balancing if sampling times were more varied
  if (worker == 1){ //ugly, this should be in the sampler code to guarantee samples are split evenly
    combinedFutures.wait();
    combinedFutures = upcxx::make_future();
  }
  if (worker >= rankCount){
    //my turn
    double fitness = _problem->evaluateSample(sampleArrayPointer.local() + (nParameters*sampleId));
    _solver->updateEvaluation(sampleId, fitness);  
  }else{
    auto fut = sendData(worker, sampleId, sampleArrayPointer.local() + (nParameters * sampleId), nParameters).then([worker](){return upcxx::rpc(worker, evaluateSample);}).then([sampleId](double fitness){_solver->updateEvaluation(sampleId, fitness);});
    combinedFutures = upcxx::when_all(combinedFutures, fut);
  }
}

void Korali::Conduit::UPCXX::run()
{
  upcxx::init();
  rankId    = upcxx::rank_me();
  rankCount = upcxx::rank_n();
  nSamples = _solver->_sampleCount;
  nParameters = _solver->N;

  // Creating sample array in global shared memory
  cache = upcxx::new_array<double>(nParameters);
  if (rankId == 0){
    combinedFutures = upcxx::make_future();
    sampleArrayPointer = upcxx::new_array<double>(nParameters * nSamples); //not really necessary for this implementation
    _solver->runSolver();
    combinedFutures.wait();
  }

  upcxx::finalize();
}


