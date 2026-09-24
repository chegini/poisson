#include "dune/common/config.h"
#include "dune/grid/config.h"
#include "dune/grid/uggrid.hh"

#include "fem/assemble.hh"
#include "fem/functional_aux.hh"
#include "fem/gridBasics.hh"
#include "fem/lagrangespace.hh"
#include "fem/spaces.hh"
#include "io/vtk.hh"
#include "utilities/gridGeneration.hh"
#include "utilities/kaskopt.hh"
#include "utilities/timing.hh"

#include <iostream>

using namespace Kaskade;


std::array<std::vector<int>,2> histogram{std::vector<int>(100),std::vector<int>(100)};

template <class Vector>
void addToHistogram(Vector const& x, int n)
{
  int N = histogram[n].size();

  auto [mn,mx] = std::minmax_element(x.begin(),x.end());
  double r = std::max(std::abs((*mn)[0]),std::abs((*mx)[0]));
  for (auto xi: x)
  {
    int i = r==0? 0: (int)( (N-1)*((xi[0]+r)/(2*r)) );
    ++histogram[n][i];
  }
}



#include "mg/bddc.hpp"

// ------------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------------



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
    DomainCache(HeatFunctional const& F,
                typename AnsatzVars::VariableSet const& vars_,
                int flags=7):
      vars(vars_)
    , sigma0(F.sigma)
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
      auto z = cell.geometry().global(xi);

      // Create manufactured solution u = x^2 * (1-x) * y^2 * (1-y).
      // - Laplace of u is -u_xx - u_yy as written below.
      double x = z[0];
      double y = z[1];
      f = - (2-6*x)*y*y*(1-y) - x*x*(1-x)*(2-6*y);

      if (x>=0.5 && y>=0.5)
        sigma = sigma0;
      else
        sigma = 1;
    }

    Scalar d0() const
    {
      return du[0]*du[0] / 2 - f*u;
    }

    template<int row>
    Scalar d1_impl(Kaskade::VariationalArg<Scalar,dim,TestVars::template Components<row>::m> const& arg) const
    {
      return du[0]*sigma*arg.derivative[0] - f*arg.value;
    }

    template<int row, int col>
    Scalar d2_impl(Kaskade::VariationalArg<Scalar,dim,TestVars::template Components<row>::m> const &arg1,
                   Kaskade::VariationalArg<Scalar,dim,AnsatzVars::template Components<row>::m> const &arg2) const
    {
      return arg1.derivative[0]*sigma*arg2.derivative[0];
    }

  private:
    typename AnsatzVars::VariableSet const& vars;
    Dune::FieldVector<Scalar,1> u, f;
    Dune::FieldMatrix<Scalar,1,dim> du;
    Cell cell;
    double sigma0, sigma;
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
      penalty = (   x[0]<=1e-10 || x[0]>=1-1e-9                // domain boundary: Dirichlet penalty
                 || x[1]<=1e-10 || x[1]>=1-1e-9) ? 1e8: 0.0;
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

  //
  // sigma: in to top-right corner [.5,1]x[.5,1] we use diffusion constant sigma instead of 1.
  HeatFunctional(double sigma_)
  : sigma(sigma_)
  {}

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

private:
  double sigma;
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
  {
  }

  //
  // sigma: in to top-right corner [.5,1]x[.5,1] we use diffusion constant sigma instead of 1.
  //
  void init(double sigma)
  {
    assert(sigma>0);

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
    Functional J(sigma);
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


void writeParallelVTK(std::ostream& out,
                      std::vector<std::string> const& files,
                      std::vector<std::tuple<std::string,std::string,int>> const& pointData)
{
  out << "<?xml version=\"1.0\"?>\n"
      << "<VTKFile type=\"PUnstructuredGrid\" >\n"
      << "  <PUnstructuredGrid GhostLevel=\"0\">\n"
      << "    <PPointData>\n";
  for (auto const& [type,name,comp]: pointData)
    out << "      <PDataArray type=\"" << type << "\" Name=\"" << name 
        << "\" NumberOfComponents=\"" << comp << "\" />\n";
  out << "    </PPointData>\n";
//     <PCellData>
//     </PCellData>
  out << "    <PPoints>\n"
      << "      <PDataArray type=\"Float32\" Name=\"Coordinates\" NumberOfComponents=\"3\"  />\n"
      << "    </PPoints>\n";
  for (auto const& f: files)
    out << "    <Piece Source=\"" << f << ".vtu\" />\n";
  out << "  </PUnstructuredGrid>\n"
      << "</VTKFile>\n";
}

// ------------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------------


int main(int argc, char** argv)
{
  using namespace Kaskade::BDDC;
  int const dim = 2;


  int iter, interfaceTypes, n, m, threads;
  bool symmetricGrid, timing, vtk;
  double sigma, tolerance;
  std::string prefix;
  if (getKaskadeOptions(argc,argv,Options
  ("prefix",      prefix,         ".","output path for vtu files")
  ("timing",     timing,         true,"whether to write timing info")
  ("vtk",           vtk,         true,"write VTK output")
  ("threads", threads,             -1,"number of threads to use in the thread pool. 0=sequential, <0=default")
  ("subres",        n,              2,"unidirectional subdomain number")
  ("elmres",        m,              2,"unidirectional element number per subdomain")
  ("iter",      iter,              16,"number of BDDC iterations")
  ("tol",       tolerance,       1e-10,"BDDC residual tolerance")
  ("crisscross", symmetricGrid,  true,"meshing with criss-cross grid")
  ("interfacetypes",interfaceTypes, 7,"bit flags for coarse interfaces to include: 1 corner 2 edge 4 face")
  ("sigma",      sigma,        1.0,   "diffusion constant to be used in top right corner [0.5,1]^2")
  )) return 0;

  if (threads >= 0)
    NumaThreadPool::instance(threads);

  auto& timer = Timings::instance();
  
  int const N = n*n; // number of subdomains

  // Subdivide the unit square into n x n square subdomains, sonsisting of m x m squares subdivided into triangles.
  // We subdivide grid creation (on construction) and space creation/assembly since UG grid creation is not
  // thread-safe. So grids are created sequentially (what is a pity, because that's excessively slow) and
  // the remaining Kaskade part of the initialization is done in parallel.
  timer.start("FE grid creation");
  std::vector<MySubdomain> subdomains;
  for (int i=0; i<n; ++i)
    for (int j=0; j<n; ++j)
    {
      Dune::FieldVector<double,dim> bl{(double)i/n,(double)j/n}, tr{(i+1.0)/n,(j+1.0)/n};
      double dx = 1.0/(n*m);

      subdomains.push_back(MySubdomain(bl,tr,dx,symmetricGrid));
    }
  timer.stop("FE grid creation");

  timer.start("subdom creation");
  parallelFor(0,N,[&](int i)
  {
    subdomains[i].init(sigma);
  });
  timer.stop("subdom creation");

  timer.start("BDDC");

  // Lexical comparison of floating point vectors - admitting some leeway for roundoff
  auto fpless = [](auto const& x, auto const& y)
  {
    if (x[0]+1e-10 < y[0]) return true;
    if (x[0] > y[0]+1e-10) return false;
    return x[1]+1e-10 < y[1];
  };
  
  timer.start("interface creation");
  // register local dofs by their spatial node position
  std::map<Dune::FieldVector<double,dim>,std::vector<BDDC::LocalDof>,
           decltype(fpless)> globalDofs(fpless);
  std::vector<int> subdomSize;
  for (int i=0; i<N; ++i)
  {
    for (auto const& [k,x]: subdomains[i].boundaryNodes())
      globalDofs[x].push_back({i,k});
    subdomSize.push_back(subdomains[i].rhs().N());
  }

  std::vector<std::vector<BDDC::LocalDof>> sharedDofs;
  for (auto const& [x,s]: globalDofs)
    if (s.size() > 1)                   // consider this only if more than one local dof is affected
      sharedDofs.push_back(s);          // - otherwise this is a boundary node
  InterfaceAverages<1,int> ifa(sharedDofs,subdomSize,interfaceTypes);
  timer.stop("interface creation");

  using TransmissionScalar = float;
  using BddcSubdomain = Subdomain<1,double,double,SpaceTransfer<1,double,TransmissionScalar>>;

  timer.start("alg subdom creation");
  std::vector<std::unique_ptr<BddcSubdomain>> subsptr(N);
  parallelFor(0,N,[&](int i)
  {
    auto const& s = subdomains[i];
    subsptr[i] = std::make_unique<BddcSubdomain>(i,s.matrix(),ifa);
    subsptr[i]->setRhs(s.rhs());
  });
  std::vector<BddcSubdomain> subs;
  for (auto& sp: subsptr)
    subs.push_back(std::move(*sp));
  timer.stop("alg subdom creation");

  timer.start("bddc creation");
  BDDCSolver<BddcSubdomain> bddcSolver(subs,ifa.coarseConstraints());
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
      for (int i=0; i<subs.size(); ++i)
      {
        auto u = subs[i].getSolution();
        files.push_back(prefix+"/sub-"+paddedString(i)+paddedString(k));
        subdomains[i].write(u,files.back());
        
        filesCorr.push_back(prefix+"/cor-"+paddedString(i)+paddedString(k));
        subs[i].getCorrection(u);
        subdomains[i].write(u,filesCorr.back());
        
        filesRaw.push_back(prefix+"/raw-"+paddedString(i)+paddedString(k));
        subs[i].getRawCorrection(u);
        subdomains[i].write(u,filesRaw.back());

        filesRes.push_back(prefix+"/res-"+paddedString(i)+paddedString(k));
        subs[i].getResidual(u);
        subdomains[i].write(u,filesRes.back());

        filesRestricted.push_back(prefix+"/restricted-"+paddedString(i)+paddedString(k));
        subdomains[i].write(subs[i].getRestrictedResidual(),filesRestricted.back());
      }
      std::ofstream outfile("sub-"+paddedString(k)+".pvtu");
      writeParallelVTK(outfile,files,
                       std::vector<std::tuple<std::string,std::string,int>>{{"Float64","sol",1}});
      std::ofstream outfileCorr("cor-"+paddedString(k)+".pvtu");
      writeParallelVTK(outfileCorr,filesCorr,
                       std::vector<std::tuple<std::string,std::string,int>>{{"Float64","sol",1}});
      std::ofstream outfileRaw("raw-"+paddedString(k)+".pvtu");
      writeParallelVTK(outfileRaw,filesRaw,
                       std::vector<std::tuple<std::string,std::string,int>>{{"Float64","sol",1}});
      std::ofstream outfileRes("res-"+paddedString(k)+".pvtu");
      writeParallelVTK(outfileRes,filesRes,
                       std::vector<std::tuple<std::string,std::string,int>>{{"Float64","sol",1}});
      std::ofstream outfileRestricted("restricted-"+paddedString(k)+".pvtu");
      writeParallelVTK(outfileRestricted,filesRestricted,
                       std::vector<std::tuple<std::string,std::string,int>>{{"Float64","sol",1}});
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

  auto [nIter, nTransfer, iB, rB, pB, cB] = bddcSolver.traffic();
  std::cout << "average amount of data exchanged per iteration: init=" << iB/(nIter*1024.0) << "kb "
            << "restriction=" << rB/(nIter*1024.0) << "kb "
            << "prolongation=" << pB/(nIter*1024.0) << "kb "
            << "coarseGrid=" << cB/(nIter*1024.0) << "kb \n"
            << "average message size: " << (iB+rB+pB)/(double)nTransfer << "b \n";

  if (timing)
    std::cout << timer;

std::ofstream hist("histogram.gnu");
for (int i=0; i<histogram[0].size(); ++i)
  hist << i << ' ' << histogram[0][i] << ' ' << histogram[1][i] << '\n';

  return 0;
}
