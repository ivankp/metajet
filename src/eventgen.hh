#ifndef METAJET_BENCHMARK_EVENTGEN_HH
#define METAJET_BENCHMARK_EVENTGEN_HH

#include <random>

class eventgen {
  std::mt19937 gen; // mersenne twister random number generator
  std::uniform_real_distribution<double> dist;

public:
  eventgen();
  void operator()(double& px, double& py, double& pz, double& E);
};

#endif
