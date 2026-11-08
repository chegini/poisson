#ifndef POISSON_GEOMETRY_HH
#define POISSON_GEOMETRY_HH
#include <array>
#include <limits>
#include <stdexcept>
namespace PoissonGeometry {
inline int subdomainCount(int n, int dimension) {
  if (n<1 || (dimension!=2 && dimension!=3)) throw std::invalid_argument("Poisson requires positive subres and dimension 2 or 3");
  int count=1;
  for (int d=0; d<dimension; ++d) {
    if (count>std::numeric_limits<int>::max()/n) throw std::invalid_argument("Too many Poisson subdomains");
    count*=n;
  }
  return count;
}
template<int dim> std::array<int,dim> coordinates(int id,int n) {
  std::array<int,dim> result{};
  for (int d=dim-1; d>=0; --d) { result[d]=id%n; id/=n; }
  return result;
}
template<int dim,class Position> double source(Position const& z) {
  double x=z[0], y=z[1];
  if constexpr (dim==2)
    return -(2-6*x)*y*y*(1-y)-x*x*(1-x)*(2-6*y); // preserve historical 2D arithmetic
  else {
    static_assert(dim==3);
    double q=z[2];
    return -(2-6*x)*y*y*(1-y)*q*q*(1-q)
           -x*x*(1-x)*(2-6*y)*q*q*(1-q)
           -x*x*(1-x)*y*y*(1-y)*(2-6*q);
  }
}
template<int dim,class Position> bool outerBoundary(Position const& x) {
  for (int d=0; d<dim; ++d) if (x[d]<=1e-10 || x[d]>=1-1e-9) return true;
  return false;
}
template<int dim,class Position> bool localBoundary(Position const& x,Position const& lo,Position const& hi) {
  for (int d=0; d<dim; ++d) if (x[d]<=lo[d]+1e-10 || x[d]>=hi[d]-1e-10) return true;
  return false;
}
template<int dim,class Position> bool positionLess(Position const& x,Position const& y) {
  for (int d=0; d<dim; ++d) {
    if (x[d]+1e-10<y[d]) return true;
    if (x[d]>y[d]+1e-10) return false;
  }
  return false;
}
}
#endif
