// Written by Ivan Pogrebnyak @ MSU

#ifndef META_JET_HH
#define META_JET_HH

#include <vector>
#include <list>
#include <cmath>
#include <iterator>
#include <limits>

namespace metajet {

constexpr double twopi   = M_PI*2;
constexpr double max_rap = 1e5; // From FastJet
constexpr double max_pt  = 1e100;

template <typename T> inline T sq(T x) noexcept { return x*x; }
template <typename T, typename... TT>
inline T sq(T x, TT... xx) noexcept { return sq(x)+sq(xx...); }

#define small_fcn(name, expr) \
  template <typename T> inline auto \
    name(const T& p) noexcept -> decltype(expr) \
    { return expr; }

small_fcn( px, p[0] )
small_fcn( py, p[1] )
small_fcn( pz, p[2] )
small_fcn( E , p[3] )

small_fcn( pt2, sq(px(p),py(p)) )

#undef small_fcn

template <typename T>
inline double phi(const T& p) noexcept {
  double _phi = atan2(py(p),px(p));
  if (__builtin_expect(_phi < 0.,0)) _phi += twopi;
  else if (__builtin_expect(_phi >= twopi,0)) _phi -= twopi;
  return _phi;
}
template <typename T>
inline double phi(const T& p, double pt2) noexcept {
  return ( __builtin_expect(pt2==0.,0) ? 0. : phi(p) );
}

template <typename T>
double rap(const T& p, double pt2) noexcept {
  // Replicated FastJet rapidity calculation
  // for compatibility in cases of unphysical 4-momenta

  double _rap;
  const double abs_pz = abs(pz(p));

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
    _rap = 0.5*log(m2_pt2/sq(E(p)+abs_pz));
    if (pz(p) > 0.) _rap = -_rap;
  }
  return rap;
}

// --------------------------------------------------------

template <typename Iter>
using deref_t = std::iterator_traits<Iter>::value_type

template <typename P> struct pseudo_jet {
  P p; // particle
  double Rij_, diB_, dij_;

  pseudo_jet(P p): p(p) { }

  inline void operator+=(const pseudo_jet<P>& other) noexcept {
    p += other.p;
  }
  double Rij(const pseudo_jet<P>& j) noexcept {
    const double pt2i = pt2(p);
    const double pt2j = pt2(j.p);
    return Rij_ = sq( rap(p,pt2i)-rap(p,pt2j), phi(p,pt2i)-phi(p,pt2j) );
  }
  inline double diB(double power) noexcept {
    // p = -1 : anti-kt
    // p =  0 : Cambridge/Aachen
    // p =  1 : kt
    return diB_ = pow(pt2(p),power);
  }
  inline double dij() noexcept {
    return dij_ = diB_ * Rij_;
  }
};

// --------------------------------------------------------

template <typename InputIterator>
std::vector<deref_t<InputIterator>>
cluster(InputIterator first, InputIterator last, double power) {
  using ptype = deref_t<InputIterator>;

  std::list<pseudo_jet<ptype>> pp; // particles and pseudo-jets
  std::vector<ptype> jj; // complete jets

  double min_diB = std::numeric_limits<double>::max();
  double min_dij = std::numeric_limits<double>::max();

  for (auto it=first; it!=last; ++it) {
    pp.emplace_back(*it);
    double diB = pp.back().diB(power);
    if (diB < min_diB) min_diB = diB;
  }

  // TODO: https://github.com/ivankp/n2jet_java/blob/master/ClusterSequence.java

}

} // end namespace

#endif
