#include "running_stats.hh"

void running_stats::push(double x) noexcept {
  n++;
  // See Knuth TAOCP vol 2, 3rd edition, page 232
  if (n == 1) {
    oldM = newM = x;
    oldS = 0.;
  } else {
    newM = oldM + (x - oldM)/n;
    newS = oldS + (x - oldM)*(x - newM);
    // set up for next iteration
    oldM = newM;
    oldS = newS;
  }
}
