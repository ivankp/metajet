// Written by Ivan Pogrebnyak @ MSU

#ifndef METAJET_HH
#define METAJET_HH

#include <cstdint>
#include <cstdlib>
#include <cmath>
#include <limits>
#include <iterator>
#include <utility>
#include <algorithm>
#include <stdexcept>
#include <string>
#include <vector>

namespace metajet {

constexpr double twopi   = M_PI*2;
constexpr double max_rap = 1e5; // From FastJet
constexpr double max_pt  = 1e100;

template <typename T> inline T sq(T x) noexcept { return x*x; }
template <typename T, typename... TT>
inline T sq(T x, TT... xx) noexcept { return sq(x)+sq(xx...); }

template <typename Iter>
using deref_t = typename std::iterator_traits<Iter>::value_type;

// ###########################################################################

#define MJ_PPROP_PP(name, expr) \
  template <typename T> inline auto \
    name(const T& p) noexcept -> decltype(expr) \
    { return expr; }

MJ_PPROP_PP( px, p[0] )
MJ_PPROP_PP( py, p[1] )
MJ_PPROP_PP( pz, p[2] )
MJ_PPROP_PP( E , p[3] )

MJ_PPROP_PP( pt2, sq(px(p),py(p)) )

template <typename T>
inline double phi(const T& p) noexcept {
  double _phi = std::atan2(py(p),px(p));
  if (__builtin_expect(_phi < 0.,0)) _phi += twopi;
  else if (__builtin_expect(_phi >= twopi,0)) _phi -= twopi;
  return _phi;
}
template <typename T>
inline double phi(const T& p, double pt2) noexcept {
  return ( __builtin_expect(pt2==0.,0) ? 0. : phi(p) );
}
inline double delta_phi(double phi_1, double phi_2) noexcept {
  double dphi = std::abs(phi_1-phi_2);
  return ( __builtin_expect(dphi > M_PI,0) ? twopi-dphi : dphi );
}

template <typename T>
double rap(const T& p, double pt2) noexcept {
  // Replicated FastJet rapidity calculation
  // for compatibility in cases of unphysical 4-momenta

  double _rap;
  const double abs_pz = std::abs(pz(p));

  if (__builtin_expect(E(p) == abs_pz && pt2 == 0.,0)) {
    // Point has infinite rapidity -- convert that into a very large
    // number, but in such a way that different 0-pt momenta will have
    // different rapidities (so as to lift the degeneracy between
    // them) [this can be relevant at parton-level]
    _rap = max_rap + abs_pz;
    if (pz(p) < 0.) _rap = -_rap;
  } else {
    // get the rapidity in a way that's modestly insensitive to roundoff
    // error when things pz,E are large (actually the best we can do without
    // explicit knowledge of mass) and force non tachyonic mass
    double m2_pt2 = (E(p)+pz(p))*(E(p)-pz(p));
    if (m2_pt2 < pt2) m2_pt2 = pt2;
    _rap = 0.5*std::log(m2_pt2/sq(E(p)+abs_pz));
    if (pz(p) > 0.) _rap = -_rap;
  }
  return _rap;
}

// ###########################################################################

struct pseudo_jet {
  double p[4]; // 4-momentum
  double Rij2, diB, dij;

  template <typename P>
  pseudo_jet(const P& p): p{px(p),py(p),pz(p),E(p)} { }

  inline void operator+=(const pseudo_jet& j) noexcept {
    p[0] += j.p[0];
    p[1] += j.p[1];
    p[2] += j.p[2];
    p[3] += j.p[3];
  }
  double calc_Rij2(const pseudo_jet& j) noexcept {
    const double pt2i = pt2(p);
    const double pt2j = pt2(j.p);
    return sq( rap(p,pt2i)-rap(j.p,pt2j),
               delta_phi(phi(p,pt2i),phi(j.p,pt2j)) );
  }
  inline void update_diB(double power) noexcept {
    // p = -1 : anti-kt
    // p =  0 : Cambridge/Aachen
    // p =  1 : kt
    diB = std::pow(pt2(p),power);
  }
  inline void update_dij(const pseudo_jet& j, double R2) noexcept {
    dij = std::min(diB,j.diB) * Rij2 / R2;
  }
};

// ###########################################################################

template <typename PseudoJet=pseudo_jet,
          typename IndexType=std::uint16_t,
          typename FastIndexType=std::uint_fast16_t>
class cluster_sequence {
public:
  using pseudo_jet_t = PseudoJet;
  using index_t = IndexType;
  using fasti_t = FastIndexType;

  static constexpr auto index_size = sizeof(index_t);
  static constexpr auto  pjet_size = sizeof(pseudo_jet_t);

private:
  const double power, R2;

  index_t capacity; // maximum number of particles that
                    // would fit without reallocation

  index_t *dij_heap, *diB_heap, *near;
  pseudo_jet_t *pp; // pseudo-jets

  static index_t validate_capacity(index_t n) {
    if (n <= std::numeric_limits<index_t>::max() ) return n;
    throw std::out_of_range(
      "attempting to allocate cluster_sequence container larger then "
      + std::to_string(std::numeric_limits<index_t>::max())
      + " (" + std::to_string(n) + ')' );
  }

  void alloc() {
    void *mem = malloc(index_size*capacity*3 + pjet_size*capacity);
    if ( !mem ) throw std::runtime_error(
      "cannot allocate memory for cluster_sequence");
    dij_heap = reinterpret_cast<index_t*>(mem);
    diB_heap = dij_heap + capacity;
    near     = diB_heap + capacity;
    pp       = reinterpret_cast<pseudo_jet_t*>(near + capacity);
  }

public:
  cluster_sequence(double power, double R)
  : power(power), R2(R*R), capacity(0), dij_heap(nullptr) { }
  cluster_sequence(double power, double R, index_t n)
  : power(power), R2(R*R), capacity(validate_capacity(n)) { alloc(); }
  ~cluster_sequence() { free(dij_heap); }

private:
  double update_Rij2(fasti_t i, fasti_t j) noexcept(noexcept(
    std::declval<pseudo_jet_t&>().calc_Rij2(
      std::declval<const pseudo_jet_t&>() )))
  {
    const double Rij2 = pp[i].calc_Rij2(pp[j]);
    if (Rij2 < pp[i].Rij2) {
      pp[i].Rij2 = Rij2;
      near[i] = j;
    }
    return Rij2;
  }

  inline void update_diB(fasti_t i) noexcept(noexcept(
    std::declval<pseudo_jet_t&>().update_diB( std::declval<double>() )))
  {
    pp[i].update_diB(power);
  }

  inline void update_dij(fasti_t i) noexcept(noexcept(
    std::declval<pseudo_jet_t&>().update_dij(
      std::declval<const pseudo_jet_t&>(), std::declval<double>() )))
  {
    pp[i].update_dij(pp[near[i]], R2);
  }

public:
  // =========================================================================

  template <typename InputIterator,
            typename OutputType=deref_t<InputIterator>>
  std::vector<OutputType>
  cluster(InputIterator first, InputIterator last) {
    std::vector<OutputType> jj; // complete jets

    const fasti_t nump = validate_capacity(std::distance(first,last));
    if (nump > capacity) {
      // allocate container memory
      if (dij_heap) free(dij_heap);
      capacity = nump;
      alloc();
    }

    fasti_t i = 0;
    for (auto it=first; it!=last; ++it) {
      auto* p = new (pp+i) pseudo_jet_t(*it); // copy particle
      p->update_diB(power); // update diB
      p->Rij2 = std::numeric_limits<double>::max(); // reset Rij
      ++i;
    }

    // update nearest geometric neighbors Rij
    for (fasti_t i=1; i<nump; ++i) {
      for (fasti_t j=0; j<i; ++j) {
        const double Rij2 = update_Rij2(i,j);
        if (Rij2 < pp[j].Rij2) {
          pp[j].Rij2 = Rij2;
          near[j] = i;
        }
      }
    }

    // update nearest scaled neighbors dij
    for (fasti_t i=0; i<nump; ++i) update_dij(i);

    // populate heaps
    for (index_t i=0; i<nump; ++i) dij_heap[i] = i;
    for (index_t i=0; i<nump; ++i) diB_heap[i] = i;

    // loop until pseudo-jets are used up -------------------
    for (fasti_t n=nump; n>1;) {

      // compute order in the heaps
      std::make_heap(dij_heap,dij_heap+n, [this](index_t i, index_t j){
        return (this->pp[i].dij > this->pp[j].dij); });
      std::make_heap(diB_heap,diB_heap+n, [this](index_t i, index_t j){
        return (this->pp[i].diB > this->pp[j].diB); });

      if ( pp[dij_heap[0]].dij < pp[diB_heap[0]].diB ) { // MERGE

        const fasti_t p = dij_heap[0], near_p = near[p];

        // merge p into near_p
        pp[near_p] += pp[p];
        update_diB(near_p); // recompute beam distances diB
        pp[near_p].Rij2 = std::numeric_limits<double>::max(); // reset Rij
        near[near_p] = near_p; // because being near is not
                               // necessarily commutative

        --n;
        // remove p from dij heap
        dij_heap[0] = dij_heap[n];
        // remove p from diB heap
        fasti_t p_in_diB_heap = 0;
        while (diB_heap[p_in_diB_heap]!=p) ++p_in_diB_heap;
        diB_heap[p_in_diB_heap] = diB_heap[n];

        // update nearest geometric neighbors Rij
        if (__builtin_expect(n>1,1)) {
          for (fasti_t i=0; i<n; ++i) {
            const fasti_t pi = dij_heap[i], near_pi = near[pi];
            if (near_pi!=p && near_pi!=near_p) continue;
            pp[pi].Rij2 = std::numeric_limits<double>::max(); // reset Rij
            for (fasti_t j=0; j<n; ++j) {
              if (j!=i) {
                const fasti_t pj = dij_heap[j], near_pj = near[pj];
                const double Rij2 = update_Rij2(pi,pj);
                // NOTE: seems to be necessary only for kt
                if (near_pj!=p && near_pj!=near_p && pp[pj].Rij2>Rij2) {
                  pp[pj].Rij2 = Rij2;
                  near[pj] = pi;
                  update_dij(pj);
                }
              }
            }
            update_dij(pi);
          }
        }

      } else { // JET

        const fasti_t p = diB_heap[0];

        // identify as jet
        const auto* mom = pp[p].p;
        jj.emplace_back(mom[0],mom[1],mom[2],mom[3]);

        --n;
        // remove p from diB heap
        diB_heap[0] = diB_heap[n];
        // remove p from dij heap
        fasti_t p_in_dij_heap = 0;
        while (dij_heap[p_in_dij_heap]!=p) ++p_in_dij_heap;
        dij_heap[p_in_dij_heap] = dij_heap[n];

        // update nearest geometric neighbors Rij
        for (fasti_t i=0; i<n; ++i) {
          const fasti_t pi = dij_heap[i];
          if (near[pi]!=p) continue;
          pp[pi].Rij2 = std::numeric_limits<double>::max(); // reset Rij
          for (fasti_t j=0; j<n; ++j) {
            if (j!=i) {
              const fasti_t pj = dij_heap[j];
              const double Rij2 = update_Rij2(pi,pj);
              // NOTE: seems to be necessary only for kt
              if (near[pj]!=p && pp[pj].Rij2>Rij2) {
                pp[pj].Rij2 = Rij2;
                near[pj] = pi;
                update_dij(pj);
              }
            }
          }
          update_dij(pi);
        }

      }

    } // end loop

    // identify last pseudo-jet as a jet
    if (__builtin_expect(nump>0,1)) {
      const auto* p = pp[dij_heap[0]].p;
      jj.emplace_back(p[0],p[1],p[2],p[3]);
    }

    return std::move(jj);
  }

};

// ###########################################################################

} // end namespace metajet

#endif
