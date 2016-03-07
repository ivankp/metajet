#include <iostream>
#include <iomanip>
#include <cstring>
#include <chrono>

#include <boost/program_options.hpp>

#include <fastjet/ClusterSequence.hh>
#include "metajet.hh"

#include "eventgen.hh"
#include "running_stats.hh"

using namespace std;
using namespace std::chrono;
using namespace fastjet;
namespace po = boost::program_options;

#define test(var) \
  std::cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << std::endl;

/*struct mom {
  double px, py, pz, E;
  mom(double px, double py, double pz, double E)
  : px(px), py(py), pz(pz), E(E) { }

  mom& operator+=(const mom& other) {
    px += other.px;
    py += other.py;
    pz += other.pz;
    E  += other.E;
  }
};*/

int main(int argc, char **argv)
{
  string cfname;
  double R, time_per_Np;
  int power;
  vector<unsigned> Nps;

  // options ---------------------------------------------------
  try {
    po::options_description desc("Options");
    desc.add_options()
      ("conf,c", po::value(&cfname),
       "configuration file (pos opt)")
      ("radius,r", po::value(&R)->required(),
       "clustering algorithm R parameter")
      ("power,p", po::value(&power)->required(),
       "clustering algorithm p parameter")
      ("time,t", po::value(&time_per_Np)->default_value(5.,"5"),
       "time per particle multiplicity in seconds")
      ("num,n", po::value(&Nps)->multitoken()->default_value(
        {2,4,8,16,32,64,128,256,512,1024},"2^1...2^10"),
       "number of particles per event")
    ;

    po::positional_options_description pos;
    pos.add("conf",1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv)
      .options(desc).positional(pos).run(), vm);
    if (argc == 1) {
      cout << desc << endl;
      return 0;
    }
    if (vm.count("conf")) {
      po::store( po::parse_config_file<char>(
        vm["conf"].as<string>().c_str(), desc), vm);
    }
    po::notify(vm);
  } catch (exception& e) {
    cerr << "Options: " <<  e.what() << endl;
    return 1;
  }
  // end options ---------------------------------------------------

  fastjet::ClusterSequence::print_banner();

  JetDefinition jdef(
    ( power==-1 ? antikt_algorithm
    : ( power==1 ? kt_algorithm
      : ( power==0 ? cambridge_algorithm
        : throw runtime_error("FastJet can only work with power=={1,0,-1}")
    ) ) ), R
  );

  cout << jdef.description() << endl;

  eventgen genevent;
  double px, py, pz, E;

  for (unsigned Np : Nps) {
    running_stats stats_fj, stats_mj;
    vector<PseudoJet> pp;
    pp.reserve(Np);

    long nevents = 0;
    const auto t0 = steady_clock::now();
    while (
      duration_cast<duration<double>>(
        steady_clock::now() - t0
      ).count() < time_per_Np
    ) {
      pp.clear();
      for (unsigned p=0; p<Np; ++p) {
        genevent(px, py, pz, E);
        pp.emplace_back(px, py, pz, E);
      }

      // FastJet
      auto t1 = high_resolution_clock::now();
      vector<PseudoJet> jets_fj = ClusterSequence(pp,jdef,false).inclusive_jets();
      stats_fj.push( duration_cast<nanoseconds>(
        high_resolution_clock::now() - t1
      ).count()/1000.);

      // metajet
      t1 = high_resolution_clock::now();
      vector<PseudoJet> jets_mj = metajet::cluster(pp.begin(),pp.end(),R,power);
      stats_mj.push( duration_cast<nanoseconds>(
        high_resolution_clock::now() - t1
      ).count()/1000.);

      ++nevents;
    }

    cout << setw(5) << Np << ' ' << setw(9) << nevents << ' '
         << setw(9) << stats_fj.mean() << " ± "
         << setw(9) << stats_fj.stdev() << " μs   "
         << setw(9) << stats_mj.mean() << " ± "
         << setw(9) << stats_mj.stdev() << " μs"
         << setw(9) << (stats_mj.mean()/stats_fj.mean()) << endl;
  }

}
