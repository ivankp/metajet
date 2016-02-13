#include <iostream>
#include <iomanip>
#include <sstream>
#include <array>
#include <vector>
#include <functional>
#include <chrono>
#include <random>
#include <future>
#include <cstring>
#include <cstdio>
#include <cmath>

#include <fastjet/ClusterSequence.hh>

#include "metajet.hh"

using namespace std;
using namespace fastjet;

template <typename T> inline T sq(T x) noexcept { return x*x; }
template <typename T, typename... TT>
inline T sq(T x, TT... xx) noexcept { return sq(x)+sq(xx...); }

int main(int argc, char **argv)
{
  if (argc!=3 && argc!=4) {
    cout << "Usage: " << argv[0] << " algorithm radius N" << endl;
    return 1;
  }

  fastjet::ClusterSequence::print_banner();

  const double R = atof(argv[2]);
  JetAlgorithm jalg;
  double power;
  function<double(const PseudoJet& a, const PseudoJet& b)> dij;
  function<double(const PseudoJet& a)> diB;

  if (!strcmp(argv[1],"antikt")) {
    jalg = antikt_algorithm;
    power = -1;
    dij = [R](const PseudoJet& a, const PseudoJet& b){
      return min(1./a.kt2(),1./b.kt2())*sq(a.delta_R(b)/R);
    };
    diB = [](const PseudoJet& a){ return 1./a.kt2(); };
  } else if (!strcmp(argv[1],"kt")) {
    jalg = kt_algorithm;
    power = 1;
    dij = [R](const PseudoJet& a, const PseudoJet& b){
      return min(a.kt2(),b.kt2())*sq(a.delta_R(b)/R);
    };
    diB = [](const PseudoJet& a){ return a.kt2(); };
  } else if (!strcmp(argv[1],"cambridge")) {
    jalg = cambridge_algorithm;
    power = 0;
    dij = [R](const PseudoJet& a, const PseudoJet& b){
      return sq(a.delta_R(b)/R);
    };
    diB = [](const PseudoJet& a){ return 1.; };
  } else {
    cerr << "Undefined algorithm " << argv[1] << endl;
    return 1;
  }

  JetDefinition jdef(jalg,R);

  // mersenne twister random number generator
  mt19937 gen(chrono::system_clock::now().time_since_epoch().count());
  uniform_real_distribution<double> dist(0.0,1.0);

  long long event = 0;
  size_t np = (argc>3 ? atoi(argv[3]) : 0);
  double px, py, pz, E;
  vector<PseudoJet> pp;
  while (true) {

    // Generate random physical 4-vectors
    if (np) {
      pp.clear();
      pp.reserve(np);
      for (size_t i=0, n=np; i<n; ++i) {
        double m, pt,
               eta = 10.*(acos(1.-2.*dist(gen))/M_PI-0.5),
               phi = 2.*M_PI*dist(gen);
        while (!isfinite(pt = 10.-150.*log(1.-dist(gen)))) { }
        while (!isfinite(m  =     -20.*log(1.-dist(gen)))) { }

        px = pt*cos(phi);
        py = pt*sin(phi);
        pz = pt*sinh(eta);
        E  = sqrt(sq(px,py,pz,m));

        // cout << px <<' '<< py <<' '<< pz <<' '<< E << endl;

        pp.emplace_back(px,py,pz,E);
      }
    } else {
      while (cin >> px >> py >> pz >> E) pp.emplace_back(px,py,pz,E);
      np = pp.size();
    }

    for (const auto& p : pp) {
      test(p.kt2())
      test(p.rap())
      test(p.phi())
    }
    cout << "R10^2 = " << sq(pp[1].delta_R(pp[0])) << endl;
    cout << "R20^2 = " << sq(pp[2].delta_R(pp[0])) << endl;
    cout << "R21^2 = " << sq(pp[2].delta_R(pp[1])) << endl;
    cout << "d10 = " << dij(pp[1],pp[0]) << endl;
    cout << "d20 = " << dij(pp[2],pp[0]) << endl;
    cout << "d21 = " << dij(pp[2],pp[1]) << endl;
    cout << endl;

    // FastJet **********************************************

    auto fut_fj = async(launch::async, [&jdef,&pp](){
      // cluster jets
      vector<PseudoJet> jets =
        ClusterSequence(pp, jdef, false).inclusive_jets();

      // sort jets by pT
      sort( jets.begin(), jets.end(),
        [](const PseudoJet& i, const PseudoJet& j){ return i.pt() > j.pt(); }
      );

      // return jets' pT
      size_t njets = jets.size();
      vector<double> pts(njets);
      for (size_t i=0; i<njets; ++i) pts[i] = jets[i].pt();
      return pts;
    });

    // metajet **********************************************

    auto fut_mj = async(launch::async, [=]() {
      // cluster jets
      vector<PseudoJet> jets =
        metajet::cluster(pp.begin(),pp.end(),R,power);

      // sort jets by pT
      sort( jets.begin(), jets.end(),
        [](const PseudoJet& i, const PseudoJet& j){ return i.pt() > j.pt(); }
      );

      // return jets' pT
      size_t njets = jets.size();
      vector<double> pts(njets);
      for (size_t i=0; i<njets; ++i) pts[i] = jets[i].pt();
      return pts;
    });

    // Print output *****************************************

    const vector<double> out_fj = fut_fj.get();
    const vector<double> out_mj = fut_mj.get();

    if (out_mj!=out_fj || argc==3) {
      const auto flags = cout.flags();
      cout << right << fixed << scientific << setprecision(8);
      cout << endl << "FJ:";
      for (double pt : out_fj) cout <<' '<< pt;
      cout << endl << "MJ:";
      for (double pt : out_mj) cout <<' '<< pt;
      cout << endl;
      cout.flags(flags);

      for (size_t i=0; i<np; ++i) {
        cout << "p" << i << ": ";
        cout << right << fixed << scientific << setprecision(8);
        cout << pp[i].px() << ' ' << pp[i].py() << ' '
             << pp[i].pz() << ' ' << pp[i].E ();
        cout.flags(flags);
        cout << " diB = " << diB(pp[i]) << endl;
      }
      if (np<21)
        for (size_t i=1; i<np; ++i)
          for (size_t j=0; j<i; ++j)
            cout << "d" << setw(3) << i << setw(3) << j
                 << " = " << dij(pp[i],pp[j]) << endl;
      cout << endl;

      break;

    } else {
      for (int i=0;i<9;++i) cout << '\b';
      cout << setw(9) << event;
      cout.flush();
    }

    ++event;

    if (argc==3) break;
  }

  return 0;
}
