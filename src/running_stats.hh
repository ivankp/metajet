#ifndef RUNNING_STATS_HH
#define RUNNING_STATS_HH

#include <cmath>

class running_stats {
  int n;
  double oldM, newM, oldS, newS;

public:
  running_stats(): n(0), oldM(0.), newM(0.), oldS(0.), newS(0.) { }

  void clear() noexcept { n=0; }
  void push(double x) noexcept;

  inline int      num() const noexcept { return n; }
  inline double  mean() const noexcept { return (n>0) ? newM : 0.; }
  inline double   var() const noexcept { return ( (n>1) ? newS/(n-1) : 0. ); }
  inline double stdev() const noexcept { return std::sqrt( var() ); }
};

#endif
