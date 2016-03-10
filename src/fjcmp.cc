#include <iostream>
#include <iomanip>
#include <sstream>
#include <array>
#include <vector>
#include <functional>
#include <chrono>
#include <future>
#include <cstring>
#include <cmath>

#include <fastjet/ClusterSequence.hh>

#include "eventgen.hh"

#ifdef DEBUG
#define BR std::cout << std::endl;
#define test(var) \
  std::cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << std::endl;
#endif

#include "metajet.hh"

using namespace std;
using namespace fastjet;

using metajet::sq;

ostream& operator<<(ostream& s, const PseudoJet& j) {
  s << setw(16) << j.px()
    << setw(16) << j.py()
    << setw(16) << j.pz()
    << setw(16) << j.E ();
  return s;
}

int main(int argc, char **argv)
{
  if (argc!=3 && argc!=4) {
    cout << "Usage: " << argv[0] << " algorithm radius N" << endl;
    return 1;
  }

  ostream coutN(cout.rdbuf());
  coutN << right << fixed << scientific << setprecision(8);

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
  metajet::cluster_sequence<> seq_mj(power,R);

  eventgen genevent;

  long long event = 0;
  size_t np = (argc>3 ? atoi(argv[3]) : 0);
  double px, py, pz, E;
  vector<PseudoJet> pp;
  while (true) {
    for (int i=0; i<9; ++i) cout << '\b';
    cout << setw(9) << event;
    cout.flush();

    // Generate random physical 4-vectors
    if (np) {
      pp.clear();
      pp.reserve(np);
      for (size_t i=0; i<np; ++i) {
        genevent(px, py, pz, E);
        pp.emplace_back(px, py, pz, E);
      }
    } else {
      while (cin >> px >> py >> pz >> E) pp.emplace_back(px,py,pz,E);
      np = pp.size();
    }

    #ifdef DEBUG
    BR
    for (size_t i=1; i<np; ++i) {
      test(pp[i].kt2())
      test(pp[i].rap())
      test(pp[i].phi())
      for (size_t j=0; j<i; ++j) {
        coutN << "R"<<i<<j<<"^2 = " << sq(pp[i].delta_R(pp[j])) << endl;
        coutN << "d"<<i<<j<<" = " << dij(pp[i],pp[j]) << endl;
      }
    }
    BR
    #endif

    // FastJet **********************************************

    // for (auto& p : pp) {
    //   static int i = 0;
    //   p.set_user_index(i++);
    // }

    auto fut_fj = async(launch::async, [&jdef,&pp](){
      // cluster jets
      auto seq = ClusterSequence(pp, jdef, false);
      vector<PseudoJet> jets = seq.inclusive_jets();

      // sort jets by pT
      sort( jets.begin(), jets.end(),
        [](const PseudoJet& i, const PseudoJet& j){ return i.pt() > j.pt(); }
      );

      // #ifdef DEBUG
      // cout <<flush<< "FJ Strategy: " << seq.strategy_string() << endl;
      // #endif

      /*cout << jets.front().cluster_hist_index() << endl;
      PseudoJet partner;
      for (const auto& jet : jets.front().constituents()) {
        cout << jet.cluster_hist_index();
        if (jet.has_partner(partner))
          cout << ' ' << partner.cluster_hist_index();
        cout << endl;
      }*/

      // return jets' pT
      size_t njets = jets.size();
      vector<double> pts(njets);
      for (size_t i=0; i<njets; ++i) pts[i] = jets[i].E();
      return pts;
    });

    // metajet **********************************************

    auto fut_mj = async(launch::async, [&seq_mj,&pp]() {
      // cluster jets
      vector<PseudoJet> jets = seq_mj.cluster(pp.cbegin(),pp.cend());

      // sort jets by pT
      sort( jets.begin(), jets.end(),
        [](const PseudoJet& i, const PseudoJet& j){ return i.pt() > j.pt(); }
      );

      // return jets' pT
      size_t njets = jets.size();
      vector<double> pts(njets);
      for (size_t i=0; i<njets; ++i) pts[i] = jets[i].E();
      return pts;
    });

    // Print output *****************************************

    // fut_fj.wait();
    // fut_mj.wait();

    const vector<double> out_fj = fut_fj.get();
    const vector<double> out_mj = fut_mj.get();

    if (out_mj!=out_fj || argc==3) {
      coutN << endl << "FJ:";
      for (double E : out_fj) coutN << setw(16) << E;
      coutN << endl << "MJ:";
      for (double E : out_mj) coutN << setw(16) << E;
      coutN << endl;

      for (size_t i=0; i<np; ++i) {
        coutN << 'p' << i << ':' << pp[i];
        cout << " diB = " << diB(pp[i]) << endl;
      }
      if (np<21)
        for (size_t i=1; i<np; ++i)
          for (size_t j=0; j<i; ++j)
            coutN << "d" << setw(3) << i << setw(3) << j
                 << " = " << dij(pp[i],pp[j]) << endl;
      cout << endl;

      break;
    }

    if (!np) break;
    ++event;
  }

  return 0;
}
