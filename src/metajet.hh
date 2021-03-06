// Written by Ivan Pogrebnyak @ MSU

#ifndef META_JET_HH
#define META_JET_HH

#include <cmath>
#include <iterator>
#include <limits>
#include <vector>
#include <list>
#include <cmath>

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

// --------------------------------------------------------

template <typename Iter>
using deref_t = typename std::iterator_traits<Iter>::value_type;

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

// --------------------------------------------------------

template <typename InputIterator, typename Cut,
          typename PseudoJet=pseudo_jet<deref_t<InputIterator>>>
std::vector<typename PseudoJet::p_type>
cluster(InputIterator first, InputIterator last,
        double r, double power, Cut cut
) {
  const double R2 = sq(r);

  std::list<pseudo_jet_wrap<PseudoJet>> pp; // particles and pseudo-jets
  std::vector<typename PseudoJet::p_type> jj; // complete jets

  for (auto it=first; it!=last; ++it) // fill list of PseudoJets
    pp.emplace_back(*it,power);

  // update nearest geometric neighbors Rij
  for (auto p=++pp.begin(), end=pp.end(); p!=end; ++p) {
    for (auto q=pp.begin(); q!=p; ++q) {
      const double Rij2 = p->update_Rij2(q);
      if (Rij2 < q->Rij2) {
        q->Rij2 = Rij2;
        q->near = p;
      }
    }
  }

  // update nearest scaled neighbors dij
  for (auto& p : pp) p.update_dij(R2);

  // loop until pseudo-jets are used up -------------------
  while (pp.size()>1) {

    auto p = pp.begin();
    double dist = p->diB; // minimum distance
    bool merge = false;

    // find smallest distance
    for (auto q=pp.begin(), end=pp.end(); q!=end; ++q) {
      if (q->dij < dist) { p = q; dist = q->dij; merge = true; }
      if (q->diB < dist) { p = q; dist = q->diB; merge = false; }
    }

    // Either merge or identify a jet
    if (merge) {

      // merge particles
      auto x = p->merge(power); // x = obsolete pseudo-jet

      // update nearest geometric neighbors Rij
      if (__builtin_expect(pp.size()>2,1)) {
        for (auto p1=pp.begin(), end=pp.end(); p1!=end; ++p1) {
          if ((p1->near!=p && p1->near!=x) || p1==x) continue;
          p1->Rij2 = std::numeric_limits<double>::max(); // reset Rij
          for (auto p2=pp.begin(); p2!=end; ++p2) {
            if (p2!=p1 && p2!=x) {
              const double Rij2 = p1->update_Rij2(p2);
              // NOTE: seems to be necessary only for kt
              if (p2->near!=p && p2->near!=x && p2->Rij2>Rij2) {
                p2->Rij2 = Rij2;
                p2->near = p1;
                p2->update_dij(R2);
              }
            }
          }
          p1->update_dij(R2);
        }
      }

      pp.erase(x);

    } else {

      // identify as jet if passes cut
      if ( cut(p->p) ) jj.emplace_back(std::move(p->p));

      // update nearest geometric neighbors Rij
      for (auto p1=pp.begin(), end=pp.end(); p1!=end; ++p1) {
        if (p1->near!=p || p1==p) continue;
        p1->Rij2 = std::numeric_limits<double>::max(); // reset Rij
        for (auto p2=pp.begin(); p2!=end; ++p2) {
          if (p2!=p1 && p2!=p) {
            const double Rij2 = p1->update_Rij2(p2);
            // NOTE: seems to be necessary only for kt
            if (p2->near!=p && p2->Rij2>Rij2) {
              p2->Rij2 = Rij2;
              p2->near = p1;
              p2->update_dij(R2);
            }
          }
        }
        p1->update_dij(R2);
      }

      pp.erase(p);

    }

  }

  // identify last pseudo-jet as a jet
  if ( __builtin_expect(pp.size(),1) && cut(pp.front().p) )
    jj.emplace_back(std::move(pp.front().p));

  // IDEA: check for N==1 at the beginning of the function

  return std::move(jj);
}

template <typename InputIterator,
          typename PseudoJet=pseudo_jet<deref_t<InputIterator>>>
inline std::vector<typename PseudoJet::p_type>
cluster(InputIterator first, InputIterator last, double R, double power) {
  return cluster(first, last, R, power,
    [](const typename PseudoJet::p_type& p){ return true; });
}

} // end namespace

#endif
