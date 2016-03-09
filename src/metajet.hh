// Written by Ivan Pogrebnyak @ MSU

#ifndef METAJET_HH
#define METAJET_HH

#include <cstdint>
#include <cstdlib>
#include <cmath>
#include <limits>
//#include <iterator>
#include <utility>
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

// template <typename Iter>
// using deref_t = typename std::iterator_traits<Iter>::value_type;

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

template <typename P> struct pseudo_jet {
  using p_type = P;

  p_type p; // particle
  double Rij2, diB, dij;

  pseudo_jet(const P& p): p(p) { }

  inline void operator+=(const pseudo_jet<P>& other) noexcept {
    p += other.p;
  }
  double update_Rij2(const pseudo_jet<P>& j) noexcept {
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
  inline void update_dij(const pseudo_jet<P>& j, double R2) noexcept {
    dij = std::min(diB,j.diB) * Rij2 / R2;
  }
};

// ###########################################################################

template <typename T> struct pseudo_jet_wrap: public T {
  using iter = typename std::list<pseudo_jet_wrap<T>>::iterator;
  iter near;

  pseudo_jet_wrap(const T& p, double power): T(p) {
    T::update_diB(power); // compute beam distances diB
    T::Rij2 = std::numeric_limits<double>::max(); // reset Rij
  }
  inline double update_Rij2(iter j) noexcept {
    const double Rij2 = T::update_Rij2(*j);
    if (Rij2 < this->Rij2) {
      this->Rij2 = Rij2;
      near = j;
    }
    return Rij2;
  }
  inline void update_dij(double R2) noexcept {
    T::update_dij(*near,R2);
  }
  iter merge(double power) noexcept {
    #ifdef DEBUG
    std::cout << std::fixed << std::scientific << std::setprecision(8);
    std::cout << "merged " << E(T::p) << " and " << E(near->p) << std::endl;
    test(T::Rij2)
    #endif

    (*this) += *near;

    #ifdef DEBUG
    test(E(T::p))
    std::cout << std::endl;
    #endif

    T::update_diB(power); // compute beam distances diB
    T::Rij2 = std::numeric_limits<double>::max(); // reset Rij
    return near;
  }
};

// ###########################################################################

template <typename PseudoJet,
          typename IndexTypes=std::pair< std::uint16_t, std::uint_fast16_t >>
class clustering_sequence {
public:
  using pseudo_jet_t = PseudoJet;
  using p_type  = typename pseudo_jet_t::p_type;
  using index_t = typename IndexTypes::first_type;
  using fasti_t = typename IndexTypes::second_type;

  static constexpr auto index_size = sizeof(index_t);
  static constexpr auto  pjet_size = sizeof(pseudo_jet_t);

private:
  const double R2, power;

  index_t capacity; // maximum number of particles that
                    // would fit without reallocation

  void* mem; // memory pool
  index_t *dij_heap, *diB_heap, *near;
  PseudoJet *pp; // pseudo-jets

  static index_t validate_capacity(index_t n) {
    if (n <= std::numeric_limits<index_t>::max() ) return n;
    throw std::out_of_range(
      "attempting to allocate clustering_sequence container larger then "
      + std::to_string(std::numeric_limits<index_t>::max())
      + " (" + std::to_string(n) + ')' );
  }

  void alloc() {
    if ( !(mem = malloc(index_size*capacity*3 + pjet_size*capacity)) )
      throw std::runtime_error("cannot allocate mem for clustering_sequence");
    dij_heap = static_cast<index_t*>(mem);
    diB_heap = dij_heap + index_size*capacity;
    near     = diB_heap + index_size*capacity;
    pp       = static_cast<value_t*>(near + index_size*capacity);
  }

public:
  clustering_sequence_container(double R, double power)
  : R2(R*R), power(power), capacity(0), mem(nullptr) { }
  clustering_sequence_container(double R, double power, index_t n)
  : R2(R*R), power(power), capacity(validate_capacity(n)) { alloc(); }

private:
  inline pseudo_jet_t& min_dij() const noexcept { return *(pp+dij_heap[0]); }
  inline pseudo_jet_t& min_diB() const noexcept { return *(pp+diB_heap[0]); }

  inline double update_Rij2(fasti_t i, fasti_t j)
  noexcept(noexcept(pseudo_jet_t::update_Rij2(
    std::declval<const pseudo_jet_t&>() )))
  {
    return pp[i]->update_Rij2(pp[j]);
  }

  inline void update_diB(fasti_t i)
  noexcept(noexcept(pseudo_jet_t::update_diB( std::declval<double>() )))
  {
    pp[i]->update_diB(power);
  }

  inline void update_dij(fasti_t i, fasti_t j)
  noexcept(noexcept(pseudo_jet_t::update_dij(
    std::declval<const pseudo_jet_t&>(), std::declval<double>() )))
  {
    pp[i]->update_dij(pp[j], R2);
  }

public:
  // =========================================================================

  template <typename InputIterator, typename Cut>
  std::vector<p_type>
  cluster(InputIterator first, InputIterator last, Cut cut) {

    const fasti_t nump = validate_capacity(std::distance(first,last));
    if (needed_capacity > capacity) {
      // allocate container memory
      if (mem) free(mem);
      capacity = needed_capacity;
      alloc();

      fasti_t i = 0;
      for (auto it=first; it!=last; ++it) {
        // copy particle and update diB
        ( new (pp+i) pseudo_jet_t(*it) )->update_diB(power);
        ++i;
      }
    }

    std::vector<p_type> jj; // complete jets

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

    // TODO: compute order in the heaps

    // loop until pseudo-jets are used up -------------------
    for (fasti_t n=nump; n>1; --n) {

      if ( pp[dij_heap[0]].dij < pp[diB_heap[0]].diB ) { // MERGE

        const fasti_t p = dij_heap[0], near_p = near[p];

        // merge particles
        pp[p] += pp[near_p];
        update_diB(p); // recompute beam distances diB
        pp[p].Rij2 = std::numeric_limits<double>::max(); // reset Rij

        // update nearest geometric neighbors Rij
        if (__builtin_expect(n>2,1)) {
          for (fasti_t i=0; i<n; ++i) {
            const fasti_t pi = dij_heap[i], near_pi = near[pi];
            if ((near_pi!=p && near_pi!=near_p) || pi==near_p) continue;
            pp[pi].Rij2 = std::numeric_limits<double>::max(); // reset Rij
            for (fasti_t j=0; j<n; ++j) {
              if (j!=i) {
                const fasti_t pj = dij_heap[j]
                if (pj!=near_p) {
                  const fasti_t near_pj = near[pj];
                  const double Rij2 = update_Rij2(pi,pj);
                  // NOTE: seems to be necessary only for kt
                  if (near_pj!=p && near_pj!=near_p && pp[pj].Rij2>Rij2) {
                    pp[pj].Rij2 = Rij2;
                    near_pj = pi;
                    update_dij(pj);
                  }
                }
              }
            }
            update_dij(pi);
          }
        }

        // TODO: recompute heaps' orders
        // move this above Rij loop

      } else { // JET

        const fasti_t p = diB_heap[0];

        // identify as jet if passes cut
        if ( cut(pp[p].p) ) jj.emplace_back(std::move(pp[p].p));

        // update nearest geometric neighbors Rij
        for (fasti_t i=0; i<n; ++i) {
          const fasti_t pi = diB_heap[i], near_pi = near[pi];
          if (near_pi!=p || pi==p) continue;
          pp[pi].Rij2 = std::numeric_limits<double>::max(); // reset Rij
          for (fasti_t j=0; j<n; ++j) {
            if (j!=i) {
              const fasti_t pj = dij_heap[j]
              if (pj!=p) {
                const fasti_t near_pj = near[pj];
                const double Rij2 = update_Rij2(pi,pj);
                // NOTE: seems to be necessary only for kt
                if (near_pj!=p && pp[pj].Rij2>Rij2) {
                  pp[pj].Rij2 = Rij2;
                  near_pj = pi;
                  update_dij(pj);
                }
              }
            }
          }
          update_dij(pi);
        }

        // TODO: recompute heaps' orders
        // move this above Rij loop

      }

    }

    // identify last pseudo-jet as a jet
    if ( __builtin_expect(nump>0,1) && cut(pp[dij_heap[0]].p) )
      jj.emplace_back(std::move(pp[dij_heap[0]].p));

    // IDEA: check for N==1 at the beginning of the function

    return std::move(jj);
  }

  template <typename InputIterator>
  inline std::vector<p_type>
  cluster(InputIterator first, InputIterator last) {
    return cluster(first, last, [](const p_type& p){ return true; });
  }

};

// ###########################################################################

} // end namespace metajet

#endif
