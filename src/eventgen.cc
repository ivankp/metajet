#include "eventgen.hh"

#include <chrono>
#include <cmath>

using namespace std;

eventgen::eventgen()
: gen(chrono::system_clock::now().time_since_epoch().count()),
  dist(0.,1.) { }

void eventgen::operator()(double& px, double& py, double& pz, double& E) {
  double m, pt,
         eta = 10.*(acos(1.-2.*dist(gen))/M_PI-0.5),
         phi = 2.*M_PI*dist(gen);
  while (!isfinite(pt = 10.-150.*log(1.-dist(gen)))) { }
  while (!isfinite(m  =     -20.*log(1.-dist(gen)))) { }

  px = pt*cos(phi);
  py = pt*sin(phi);
  pz = pt*sinh(eta);
  E  = sqrt(px*px+py*py+pz*pz+m*m);
}
