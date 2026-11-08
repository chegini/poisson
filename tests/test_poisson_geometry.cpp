#include "poisson_geometry.hh"
#include <array>
#include <cmath>
#include <iostream>
#include <set>

void check(bool condition) { if (!condition) throw std::runtime_error("Poisson geometry regression"); }
int main() {
  using namespace PoissonGeometry;
  check(subdomainCount(2,3)==8 && subdomainCount(4,2)==16);
  bool rejected=false;try {subdomainCount(50000,3);} catch(std::invalid_argument const&) {rejected=true;}
  check(rejected);
  std::set<std::array<int,3>> ids;
  for(int i=0;i<27;++i) ids.insert(coordinates<3>(i,3));
  check(ids.size()==27 && coordinates<3>(5,2)==std::array<int,3>{1,0,1});
  check(coordinates<2>(6,4)==std::array<int,2>{1,2});
  using Point=std::array<double,3>;
  Point middle{.5,.5,.5},lo{0,0,0},hi{1,1,1};
  check(std::abs(source<3>(middle)-3./64)<1e-15);
  check(std::abs(source<2>(middle)-.25)<1e-15);
  auto u=[](Point const& x) {double p=1;for(double a:x) p*=a*a*(1-a);return p;};
  Point x{.27,.41,.63}; double laplace=0,h=1e-4;
  for(int d=0;d<3;++d) {auto a=x,b=x;a[d]+=h;b[d]-=h;laplace+=(u(a)-2*u(x)+u(b))/(h*h);}
  check(std::abs(source<3>(x)+laplace)<1e-8);
  check(!outerBoundary<3>(middle) && !localBoundary<3>(middle,lo,hi));
  for(int d=0;d<3;++d) for(double edge:{0.,1.}) {
    auto q=middle;q[d]=edge;
    check(outerBoundary<3>(q) && localBoundary<3>(q,lo,hi));
    check(u(q)==0);
  }
  auto z=middle;z[2]=.75;
  check(positionLess<3>(middle,z) && !positionLess<3>(z,middle));
  check(!positionLess<2>(middle,z));
  std::cout << "Poisson 2D/3D geometry: source, all boundary planes, z ordering, subdomain indexing and overflow passed\n";
}
