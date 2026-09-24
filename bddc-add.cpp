#include "dune/common/config.h"
#include "dune/grid/config.h"
#include "dune/grid/uggrid.hh"

#include "fem/assemble.hh"
#include "fem/functional_aux.hh"
#include "fem/gridBasics.hh"
#include "fem/lagrangespace.hh"
#include "fem/spaces.hh"
#include "io/vtk.hh"
#include "mg/add.hh"
#include "mg/bddc.hpp"
#include "utilities/gridGeneration.hh"
#include "utilities/kaskopt.hh"
#include "utilities/timing.hh"

#include <iostream>

using namespace Kaskade;





template <class VarSet>
class HeatFunctional: public Kaskade::FunctionalBase<Kaskade::VariationalFunctional>
{
public:
  using Scalar = double;
  using OriginVars = VarSet;
  using AnsatzVars = VarSet;
  using TestVars = VarSet;
  using Grid = typename AnsatzVars::Grid;
  using GridView = typename AnsatzVars::GridView;
  using Cell = Kaskade::Cell<Grid>;

  static constexpr int dim = Grid::dimension;
  static constexpr int uIdx = 0;
  static constexpr int uSpaceIdx = spaceIndex<AnsatzVars,uIdx>;

  // The domain cache defines F as well as its first and second directional derivative,
  // evaluated at the current iterate u.
  class DomainCache : public Kaskade::CacheBase<HeatFunctional,DomainCache>
  {
  public:
    DomainCache(HeatFunctional const&,
                typename AnsatzVars::VariableSet const& vars_,
                int flags=7):
      vars(vars_)
    {}

    void moveTo(Cell const& c)
    {
      cell = c;
    }

    template <class Position, class Evaluators>
    void evaluateAt(Position const& xi, Evaluators const& evaluators)
    {
      u  = vars.template value<0>(evaluators);
      du = vars.template derivative<0>(evaluators);
      auto x = cell.geometry().global(xi);
      f = x[0] > std::abs(x[1]-0.5)+0.51 && x[1]>=0.5? 2: 0;
      kappa = 1;
    }

    Scalar d0() const
    {
      return du[0]*kappa*du[0] / 2 - f*u;
    }

    template<int row>
    Scalar d1_impl(Kaskade::VariationalArg<Scalar,dim,TestVars::template Components<row>::m> const& arg) const
    {
      return du[0]*kappa*arg.derivative[0] - f*arg.value;
    }

    template<int row, int col>
    Scalar d2_impl(Kaskade::VariationalArg<Scalar,dim,TestVars::template Components<row>::m> const &arg1,
                   Kaskade::VariationalArg<Scalar,dim,AnsatzVars::template Components<row>::m> const &arg2) const
    {
      return arg1.derivative[0]*kappa*arg2.derivative[0];
    }

  private:
    double kappa;
    typename AnsatzVars::VariableSet const& vars;
    Dune::FieldVector<Scalar,1> u, f;
    Dune::FieldMatrix<Scalar,1,dim> du;
    Cell cell;
  };


  class BoundaryCache
  {
  public:
    BoundaryCache(HeatFunctional<AnsatzVars> const& functional,
                  typename AnsatzVars::VariableSet const& vars_,
                  int flags=7):
      vars(vars_)
    {}

    template <class FaceIterator>
    void moveTo(FaceIterator const& fi)
    {
      face = *fi;
    }

    template <class Position, class Evaluators>
    void evaluateAt(Position const& xi, Evaluators const& evaluators)
    {
      u = vars.template value<0>(evaluators);
      auto x = face.geometry().global(xi);
      penalty = (   x[0]<=1e-10 || x[0]>=1-1e-10                // domain boundary: Dirichlet penalty
                 || x[1]<=1e-10 || x[1]>=1-1e-10) ? 1e2: 0.0;
    }

    Scalar d0() const
    {
      return penalty*(u-dirichlet)*(u-dirichlet)/2;
    }

    template<int row>
    Dune::FieldVector<Scalar,1> d1(Kaskade::VariationalArg<Scalar,dim> const& argT) const
    {
      return penalty*(u-dirichlet)*argT.value;
    }

    template<int row, int col>
    Dune::FieldMatrix<Scalar,1,1> d2(Kaskade::VariationalArg<Scalar,dim> const &argT,
                                     Kaskade::VariationalArg<Scalar,dim> const &argA) const
    {
      return penalty*argT.value*argA.value;
    }

  private:
    typename AnsatzVars::VariableSet const& vars;
    Scalar penalty = 1e9;
    double u;
    double dirichlet = 0;
    Face<GridView> face;
  };

  // Declare static properties of the right hand side (i.e. the first derivative of the
  // variational functional. The default properties are fine.
  template <int row>
  struct D1: public Kaskade::FunctionalBase<Kaskade::VariationalFunctional>::D1<row>
  { };

  // Declare static properties of the stiffness matrix (i.e. the second derivative of the
  // variational functional. The default properties are fine.
  template <int row, int col>
  struct D2: public Kaskade::FunctionalBase<Kaskade::VariationalFunctional>::D2<row,col>
  { };



  template <class Cell>
  int integrationOrder(Cell const& /* cell */, int shapeFunctionOrder, bool boundary) const
  {
    if (boundary)
      return 2*shapeFunctionOrder;
    else
    {
      int stiffnessMatrixIntegrationOrder = 2*(shapeFunctionOrder-1);
      int sourceTermIntegrationOrder = shapeFunctionOrder;        // as rhs f is constant, i.e. of order 0

      return std::max(stiffnessMatrixIntegrationOrder,sourceTermIntegrationOrder);
    }
  }
};


class MySubdomain
{
  constexpr static int dim = 2;
  using Grid = Dune::UGGrid<dim>;
  using GridView = Grid::LeafGridView;
  using GridMan = Kaskade::GridManager<Grid>;

public:

  using Matrix = Kaskade::NumaBCRSMatrix<Dune::FieldMatrix<double,1,1>>;
  using Vector = Dune::BlockVector<Dune::FieldVector<double,1>>;

  MySubdomain(Dune::FieldVector<double,2> bl_,
              Dune::FieldVector<double,2> tr_,
              double dx, bool symmetricGrid)
  : bl(bl_)
  , tr(tr_)
  , gridman(new GridMan(createCuboid<Grid>(bl,tr-bl,dx,symmetricGrid)))
  { }

  void init()
  {
    using namespace Kaskade::BDDC;
    gridman->enforceConcurrentReads(true);
    h1Space = std::make_unique<H1Space<Grid>>(*gridman,gridman->grid().leafGridView(),1);

    auto spaces = makeSpaceList(h1Space.get());
    auto varSetDesc = makeVariableSetDescription(spaces,boost::fusion::make_vector(
                                                 Variable<SpaceIndex<0>,Components<1>>("sub")));
    using VarSetDesc = decltype(varSetDesc);
    using Functional = HeatFunctional<VarSetDesc>;
    using Assembler = VariationalFunctionalAssembler<LinearizationAt<Functional>>;

    // Compute local stiffness matrix and right hand side.
    Functional J;
    Assembler assembler(spaces);
    auto u = varSetDesc.variableSet();
    assembler.assemble(linearization(J,u));
    f = component<0>(assembler.template rhs<>());
    f *= -1; // Newton sign
    A = assembler.template get<Matrix>(false);

    // Find boundary nodes and their indices.

    auto gridView = h1Space->gridView();
    for (auto const& v: Dune::vertices(gridView))
    {
      auto x = v.geometry().center();
      if (   x[0]<=bl[0]+1e-10 || x[0]>=tr[0]-1e-10                 // Check if vertex is on the boundary of
          || x[1]<=bl[1]+1e-10 || x[1]>=tr[1]-1e-10 )               // the rectangle (given by [bl,tr]).
        bdryVertices.push_back({gridView.indexSet().index(v),x});   // Find index based on P1 assumption:
    }                                                               // dof indexing equals vertex indexing.
  }

  Matrix const& matrix() const
  {
    return A;
  }

  Vector const& rhs() const
  {
    return f;
  }

  void write(Vector const& u_, std::string const& filename) const
  {
    auto u = h1Space->template element<1>();
    u.coefficients() = u_;
    writeVTK(u,filename,IoOptions(),"sol");
  }
  
  std::vector<std::pair<int,Dune::FieldVector<double,2>>> const& boundaryNodes() const
  {
    return bdryVertices;
  }


private:
  Dune::FieldVector<double,2> bl, tr;
  std::unique_ptr<GridMan> gridman;
  std::unique_ptr<H1Space<Grid>> h1Space;
  Matrix A;
  Vector f;
  std::vector<std::pair<int,Dune::FieldVector<double,2>>> bdryVertices;
};



                      

int main(int argc, char** argv)
{
  using namespace Kaskade::BDDC;
  int const dim = 2;


  int iter, interfaceTypes, n, m, threads, refinements;
  bool symmetricGrid, timing, vtk;
  double decompTol, tolerance;
  std::string prefix;
  if (getKaskadeOptions(argc,argv,Options
  ("refine", refinements,           6,"uniform mesh refinements")
  ("prefix",      prefix,         ".","output path for vtu files")
  ("timing",     timing,         true,"whether to write timing info")
  ("vtk",           vtk,         true,"write VTK output")
  ("threads", threads,             -1,"number of threads to use in the thread pool. 0=sequential, <0=default")
  ("subres",        n,              2,"number of subdomains")
  ("elmres",        m,              2,"unidirectional element number")
  ("iter",      iter,              16,"number of BDDC iterations")
  ("tol",       tolerance,       1e-10,"BDDC residual tolerance")
  ("crisscross", symmetricGrid,  true,"meshing with criss-cross grid")
  ("interfacetypes",interfaceTypes, 7,"bit flags for coarse interfaces to include: 1 corner 2 edge 3 face")
  ("decompTol",   decompTol,   1000.0,"threshold for relative weak diagonal dominance violation / eps")
  )) return 0;

  if (threads >= 0)
    NumaThreadPool::instance(threads);

  using Entry = Dune::FieldMatrix<double,1,1>;
  using Matrix = NumaBCRSMatrix<Entry>;
  
  auto& timer = Timings::instance();

  using Grid = Dune::UGGrid<dim>;
  using GridView = Grid::LeafGridView;
  using GridMan = Kaskade::GridManager<Grid>;

  timer.start("FE grid creation");
  GridMan gridman(createCuboid<Grid>(GlobalPosition<Grid>({0,0}),GlobalPosition<Grid>({1,1}),
                                     GlobalPosition<Grid>({1.0/m,1.0/m}),symmetricGrid));
  gridman.enforceConcurrentReads(true);
  H1Space<Grid> h1Space(gridman,gridman.grid().leafGridView(),1);

  auto spaces = makeSpaceList(&h1Space);
  auto varSetDesc = makeVariableSetDescription(spaces,boost::fusion::make_vector(
                                               Variable<SpaceIndex<0>,Components<1>>("u")));
  using VarSetDesc = decltype(varSetDesc);
  using Functional = HeatFunctional<VarSetDesc>;
  using Assembler = VariationalFunctionalAssembler<LinearizationAt<Functional>>;

  // Compute local stiffness matrix and right hand side.
  Functional J;
  Assembler assembler(spaces);
  auto u = varSetDesc.variableSet();
  assembler.assemble(linearization(J,u));
  auto f = component<0>(assembler.template rhs<>());
  f *= -1; // Newton sign
  auto A = assembler.template get<Matrix>(false);
  timer.stop("FE grid creation");

  auto uDirect = varSetDesc.variableSet();
  timer.start("whole problem inversion");
  timer.start("factorization");
  auto Ainv = directInverseOperator(makeAssembledGalerkinOperator(assembler));
  timer.stop("factorization");
  timer.start("solution");
  auto du = varSetDesc.zeroCoefficientVector();
  Ainv.apply(assembler.template rhs<>(),du);
  timer.stop("solution");
  timer.stop("whole problem inversion");
  if (vtk)
  {
    timer.start("output");
    component<0>(uDirect) -= component<0>(du);
    std::cout << "norm solution = " << component<0>(uDirect).coefficients().two_norm() << "\n";
    writeVTK(uDirect,"solution",IoOptions());
    timer.stop("output");
  }


  timer.start("subdom creation");

  auto [As,Fs,subdomIndices,sharedDofs] = matrixDomainDecommposition(A,f,n,decompTol*std::numeric_limits<double>::epsilon());
  std::vector<int> subdomSize(n);
  for (int i=0; i<n; ++i)
    subdomSize[i] = As[i].N();
  InterfaceAverages<1,int> ifa(sharedDofs,subdomSize,interfaceTypes);

  timer.stop("subdom creation");



  timer.start("BDDC");

//
  timer.start("alg subdom creation");
  std::vector<std::unique_ptr<Subdomain<1>>> subsptr(n);
  parallelFor(0,n,[&](int i)
  {
    subsptr[i] = std::make_unique<Subdomain<1>>(i,As[i],ifa);
    subsptr[i]->setRhs(Fs[i]);
  });
  std::vector<Subdomain<1>> subs;
  for (auto& sp: subsptr)
    subs.push_back(std::move(*sp));
  timer.stop("alg subdom creation");

  timer.start("bddc creation");
  BDDCSolver<Subdomain<1>> bddcSolver(subs,ifa.coarseConstraints());
  timer.stop("bddc creation");

  std::vector<double> resNorm;
  for (int k=0; k<iter; ++k)
  {
    timer.start("BDDC solve");
    resNorm.push_back(bddcSolver.solve());
    timer.stop("BDDC solve");

    if (vtk)
    {
      ScopedTimingSection out("output");
      std::vector<std::string> files, filesCorr, filesRaw, filesRes, filesRestricted;
      auto u = h1Space.element<1>();
      auto du = h1Space.element<1>();
      for (int i=0; i<n; ++i)
      {
        auto ui = subs[i].getSolution();
        auto dui = ui; subs[i].getCorrection(dui);
        for (int j=0; j<ui.N(); ++j)
        {
          u.coefficients()[subdomIndices[i][j]] = ui[j];
          du.coefficients()[subdomIndices[i][j]] = dui[j];
        }
      }
      writeVTK(u,prefix+"/sol-"+paddedString(k),IoOptions(),"u");
      writeVTK(du,prefix+"/cor-"+paddedString(k),IoOptions(),"du");
      u -= component<0>(uDirect); u *= -1;
      writeVTK(u,prefix+"/err-"+paddedString(k),IoOptions(),"err");
    }

    if (resNorm.back() <= tolerance)
      break;
  }

  timer.stop("BDDC");




  if (resNorm.size() >= 2)
  {
    int lookback = std::min(10,static_cast<int>(resNorm.size())-1);
    double older = resNorm[resNorm.size()-lookback-1];
    if (older > 0)
    {
      double contraction = std::pow(resNorm.back()/older,1.0/lookback);
      std::cout << "Estimated contraction factor: " << contraction << ". (kappa ~ " << (1+contraction)/(1-contraction) << ").\n";
    }
    else
      std::cout << "Contraction factor unavailable: reference residual is zero.\n";
  }
  else
    std::cout << "Contraction factor unavailable: fewer than two iterations.\n";



  if (timing)
    std::cout << timer;
  return 0;
}
