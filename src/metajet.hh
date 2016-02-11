// Written by Ivan Pogrebnyak @ MSU

#ifndef META_JET_HH
#define META_JET_HH

#include <cmath>
#include <iterator>
#include <limits>
#include <vector>
#include <list>

namespace metajet {

constexpr double twopi   = M_PI*2;
constexpr double max_rap = 1e5; // From FastJet
constexpr double max_pt  = 1e100;

template <typename T> inline T sq(T x) noexcept { return x*x; }
template <typename T, typename... TT>
inline T sq(T x, TT... xx) noexcept { return sq(x)+sq(xx...); }

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
  return ( __builtin_expect(dhi > M_PI,0) ? twopi-dphi : dphi );
}

template <typename T>
double rap(const T& p, double pt2) noexcept {
  // Replicated FastJet rapidity calculation
  // for compatibility in cases of unphysical 4-momenta

  double _rap;
  const double abs_pz = std::abs(pz(p));

  if (E(p) == abs_pz && pt2 == 0.) {
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
  return rap;
}

// --------------------------------------------------------

template <typename Iter>
using deref_t = std::iterator_traits<Iter>::value_type;

/* TODO: use this
template <bool With> struct with_constits { };
template <> struct with_constits<true> {
  std::vector<int> constits;
  with_constits(int id): constits(1,id) { }
};*/

template <typename P> struct pseudo_jet {
  P p; // particle
  double Rij, diB, dij;

  pseudo_jet(const P& p): p(p) { }

  inline void operator+=(const pseudo_jet<P>& other) noexcept {
    p += other.p;
  }
  bool update_Rij(const pseudo_jet<P>& j) noexcept {
    const double pt2i = pt2(p);
    const double pt2j = pt2(j.p);
    const double Rij = sq( rap(p,pt2i)-rap(p,pt2j), delta_phi(phi(p,pt2i),phi(p,pt2j)) );
    if (Rij < this->Rij) {
      this->Rij = Rij;
      return true;
    } else return false;
  }
  inline void update_diB(double power) noexcept {
    // p = -1 : anti-kt
    // p =  0 : Cambridge/Aachen
    // p =  1 : kt
    diB_ = pow(pt2(p),power);
  }
  inline void update_dij(const pseudo_jet<P>& j) noexcept {
    dij_ = std::min(diB,j.diB) * Rij;
  }
};

// --------------------------------------------------------

template <typename PseudoJet=pseudo_jet<deref_t<InputIterator>>,
          typename Cut, typename InputIterator>
std::vector<deref_t<InputIterator>>
cluster(InputIterator first, InputIterator last, double power,
        Cut cut = [](const deref_t<InputIterator>& p){ return true; }
) {

  template <typename T> struct with_nearest: public T {
    using iter = std::list<with_nearest<T>>::iterator;
    iter near;

    with_nearest(const T& p, double power): T(p) {
      update_diB(power); // compute beam distances diB
      Rij = std::numeric_limit<double>::max(); // reset Rij
    }

    inline void update_Rij(iter j) noexcept {
      if ( T::update_Rij(*j) ) near = j;
    }

    inline void update_dij() noexcept { T::update_dij(*near); }

    iter merge(double power) noexcept {
      (*this) += *near;
      update_diB(power); // compute beam distances diB
      Rij = std::numeric_limit<double>::max(); // reset Rij
      return near;
    }
  };

  std::list<with_nearest<PseudoJet>> pp; // particles and pseudo-jets
  std::vector<deref_t<InputIterator>> jj; // complete jets

  for (auto it=first; it!=last; ++it) // fill list of PseudoJets
    pp.emplace_back(*it,power);

  // update nearest geometric neighbors Rij
  for (auto p=++pp.begin(), end=pp.end(); p!=end; ++p) {
    for (auto q=pp.begin(); q!=p; ++q) {
      p->update_Rij(q);
      if (p->Rij < q->Rij) q->near = p;
    }
  }

  // update nearest scaled neighbors dij
  for (auto& p : pp) p.update_dij();

  // TODO: implement grid

  // loop until pseudo-jets are used up -------------------
  while (pp.size()) {

    auto p = pp.begin();
    double dist = p.diB; // minimum distance
    bool merge = false;

    // find smallest distance
    for (auto q=pp.begin(), end=pp.end(); q!=end; ++q) {
      if (q->dij < dist) { p = q; dist = q->dij; merge = true; }
      if (q->diB < dist) { p = q; dist = q->diB; merge = false; }
    }

    // Either merge or identify a jet
    if (merge) {

      // merge particles
      auto x = p->merge(); // x = obsolete pseudo-jet

      // update nearest geometric neighbors Rij
      for (auto p1=++pp.begin(), end=pp.end(); p1!=end; ++p1) {
        bool p1_affected = (p1->near==p || p1->near==x);
        if (p1_affected) // reset Rij
          p1->Rij = std::numeric_limit<double>::max();

        for (auto p2=pp.begin(); p2!=p1; ++p2) {
          bool p2_affected = (p2->near==p || p2->near==x);
          if (p2_affected) {
            p2->Rij = std::numeric_limit<double>::max();

            if (p1_affected) { // only if both affected
              // TODO: <<<<<<<<<<<<<<<<<<<<<<<<<<
            }
          }
          p1->update_Rij(p2);
          if (p1->Rij < p2->Rij) p2->near = p1;
        }
      }

    } else {

      // identify as jet
      if ( cut(p.p) ) jj.emplace_back(std::move(p.p)); // passes cut

      pp.erase(p); // "remove"

      // TODO: <<<<<<<<<<<<<<<<<<<<<<<<<<
      for (auto p1=++pp.begin(), end=pp.end(); p1!=end; ++p1) {
        if (p1.near==p) {
          p1.rm_near();
          for (PseudoJet p2=first; p2!=null; p2=p2.next) {
            if (p1!=p2) p1.update_near(p2,false);
          }
          p1.update_dij();
        }
      }

    }

  }

  return std::move(jj);

}

} // end namespace

#endif
