#include <array>
#include <cassert>
#include <iostream>
#include <limits>
#include "../poisson_validation.hh"

struct Vector : std::vector<std::array<double,1>>
{
  using std::vector<std::array<double,1>>::vector;
};
struct Matrix
{
  void mv(Vector const& x, Vector& y) const
  {
    y[0][0]=2*x[0][0]-x[1][0];
    y[1][0]=2*x[1][0]-x[0][0];
  }
};
struct LocalDof { int s,i; };
struct Subdomain
{
  Vector u;
  Vector getSolution() const { return u; }
};

int main()
{
  using namespace PoissonValidation;
  auto ids=makeMapping(std::vector<int>{2,2},
    std::vector<std::vector<LocalDof>>{{{0,1},{1,0}}});
  assert((ids==Mapping{{0,1},{1,2}}));
  Matrix a;
  Vector b0{{{0}},{{3}}}, b1{{{1}},{{4}}};
  std::vector<Matrix const*> matrices{&a,&a};
  std::vector<Vector const*> rhs{&b0,&b1};
  Monitor<Matrix,Vector> monitor(matrices,rhs,ids,{true,true},[](auto&) {});
  assert((monitor.globalRhs()==std::vector<double>{0,4,4}));
  auto initial=monitor.evaluateGlobal({0,0,0});
  assert(initial.relative==1 && !initial.converged(1e-7));
  std::vector<Subdomain> storage{{Vector{{{1}},{{2}}}}, {Vector{{{2}},{{3}}}}};
  auto exact=monitor.evaluate(storage);
  assert(exact.absolute==0 && exact.relativeJump==0 && exact.converged(1e-7));
  // Averaging alone could hide an interface discontinuity. It must not pass.
  storage[0].u[1][0]=1;
  storage[1].u[0][0]=3;
  auto discontinuous=monitor.evaluate(storage);
  assert(discontinuous.absolute==0 && discontinuous.relativeJump>0);
  assert(!discontinuous.converged(1e-7));

  // Emulate the other owner's contributions to each real reduction.
  int phase=0;
  Monitor<Matrix,Vector> ownerMonitor(matrices,rhs,ids,{true,false},
    [&](std::vector<double>& values)
    {
      if (phase==0 || phase==3) { values[1]+=1; values[2]+=4; }
      else if (phase==1) { values[1]+=2; values[2]+=3; }
      else assert(phase==2 && values.size()==1 && values[0]==0);
      ++phase;
    });
  storage[0].u[1][0]=2;
  auto owned=ownerMonitor.evaluate(storage);
  assert(phase==4 && owned.converged(1e-7));
  Vector zero(2);
  Monitor<Matrix,Vector> zeroMonitor(matrices,{&zero,&zero},ids,{true,true},[](auto&) {});
  auto zeroRhs=zeroMonitor.evaluateGlobal({1e-9,0,0});
  assert(zeroRhs.rhsNorm==0 && zeroRhs.relative==zeroRhs.absolute);
  assert(zeroRhs.converged(1e-7));
  auto bad=monitor.evaluateGlobal({std::numeric_limits<double>::quiet_NaN(),0,0});
  assert(!bad.converged(1e-7));
  std::cout << "Poisson global residual, shared-DOF assembly, owner reductions and continuity checks passed.\n";
}
