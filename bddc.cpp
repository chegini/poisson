#ifdef KASKADE_HAVE_MPI
#include <mpi.h>
#endif

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
#include <filesystem>

using namespace Kaskade;


#include "mg/bddc.hpp"
#include "poisson_bddc_options.hh"

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
  bool symmetricGrid, timing, vtk, histogram;
  PoissonBddcOptions bddcOptions;
  bool direct;
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
  ("mpi",         bddcOptions.mpi, false, "enable owner-only MPI BDDC")
  ("compression", bddcOptions.compression, false, "enable compressed BDDC transfers")
  ("compressionBits", bddcOptions.compressionBits, 16, "quantization bits for compression")
  ("restrictCompressionBits", bddcOptions.restrictCompressionBits, -1, "restriction bits: -1=inherit, 1..32=override")
  ("prolongateCompressionBits", bddcOptions.prolongateCompressionBits, -1, "prolongation bits: -1=inherit, 1..32=override")
  ("restrictionPasses", bddcOptions.restrictionPasses, 1, "restriction exchanges per update, 1..16")
  ("useCg", bddcOptions.useCg, false, "use CG instead of legacy steepest descent")
  ("verbose", bddcOptions.verbose, true, "print internal iteration table")
  ("relativeTolerance", bddcOptions.relativeTolerance, 0.0, "true global relative residual target; 0=legacy --tol")
  ("direct", direct, false, "solve assembled global problem with UMFPACK for validation")
  ("compressionMinInterfaceSize", bddcOptions.compressionMinInterfaceSize, 0, "minimum scalar interface size for compressed BDDC transfer; 0=all")
  ("compressionThresholdInterfaceTypes", bddcOptions.compressionThresholdInterfaceTypes, static_cast<int>(Kaskade::BDDC::FACE), "interface type mask for the size threshold: 2=edge, 4=face, 0=all")
  ("graphLifting", bddcOptions.graphLifting, false, "graph lifting before quantization (set graphLiftingAfterQuantization=0)")
  ("graphLiftingAfterQuantization", bddcOptions.graphLiftingAfterQuantization, true, "apply reversible integer graph lifting after quantization")
  ("bddcExcludeCornersFromCompression", bddcOptions.excludeCornersFromCompression, false, "exclude BDDC corner DOFs from compressed transfer payloads")
  ("huffman", bddcOptions.huffman, true, "use Huffman transfer encoding")
  ("bitlength", bddcOptions.bitlength, true, "use bit-length transfer encoding")
  ("zigzag", bddcOptions.zigzag, true, "use zigzag signed-symbol mapping")
  ("histogram", histogram, false, "write compression symbol histograms")
  )) return 0;

  validatePoissonOptions(bddcOptions);
  if (n<1 || m<1 || iter<1 || !std::isfinite(sigma) || sigma<=0 ||
      !std::isfinite(tolerance) || tolerance<=0)
    throw std::invalid_argument("subres, elmres, iter, sigma and tol must be positive");
  if (direct && bddcOptions.relativeTolerance<=0)
    throw std::invalid_argument("--direct requires --relativeTolerance > 0");
  if (bddcOptions.mpi && vtk)
    throw std::runtime_error("--vtk 1 is not supported with --mpi 1; use --vtk 0");
  PoissonMpiSession mpiSession(argc,argv,bddcOptions.mpi);
  if (poissonMpiRankZero(bddcOptions.mpi))
    std::cout << "Poisson compression order: "
              << (!bddcOptions.compression ? "uncompressed" :
                  bddcOptions.graphLiftingAfterQuantization && bddcOptions.huffman
                    ? "quantize-first" : bddcOptions.graphLifting &&
                       !bddcOptions.graphLiftingAfterQuantization
                    ? "transform-first" : "quantization-only") << '\n';


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
  auto const interfaceCount = ifa.interfaceCount();
  if (poissonMpiRankZero(bddcOptions.mpi))
    std::cout << "BDDC coarse interfaces: corners=" << interfaceCount[0]
            << ", edges=" << interfaceCount[1]
            << ", faces=" << interfaceCount[2]
            << " (interfacetypes=" << interfaceTypes << ")" << std::endl;
  timer.stop("interface creation");

  timer.start("alg subdom creation");
  using Matrix = NumaBCRSMatrix<Dune::FieldMatrix<double,1,1>>;
  std::vector<Matrix const*> matrices;
  matrices.reserve(N);
  for (auto const& subdomain: subdomains)
    matrices.push_back(&subdomain.matrix());
  timer.stop("alg subdom creation");

  using Vector = MySubdomain::Vector;
  using Monitor = PoissonValidation::Monitor<Matrix,Vector>;
  std::unique_ptr<Monitor> monitor;
  PoissonValidation::Mapping mapping;
  if (bddcOptions.relativeTolerance>0)
  {
    mapping=PoissonValidation::makeMapping(subdomSize,sharedDofs);
    std::vector<Vector const*> rhs;
    std::vector<bool> owned;
    MpiDomainLayout layout;
    layout.configure(bddcOptions.mpi,N);
    for (int i=0; i<N; ++i)
    {
      rhs.push_back(&subdomains[i].rhs());
      owned.push_back(layout.owns(i));
    }
    monitor=std::make_unique<Monitor>(matrices,rhs,mapping,owned,
      [&](auto& values) { poissonSum(values,bddcOptions.mpi); });
    if (poissonMpiRankZero(bddcOptions.mpi))
      std::cout << "Poisson global dofs: " << monitor->globalRhs().size()
                << "\nPoisson stopping rule: ||b-A*x||_2 / ||b||_2 <= "
                << bddcOptions.relativeTolerance << " (zero RHS: absolute residual); scaled interface jump <= same tolerance\n";
  }
  if (direct)
  {
    auto x=solvePoissonDirect<Matrix,Vector>(matrices,mapping,monitor->globalRhs(),bddcOptions.mpi);
    auto sample=monitor->evaluateGlobal(std::move(x));
    printPoissonResidual(sample,0,bddcOptions.mpi);
    savePoissonSolution(sample,prefix,bddcOptions.mpi);
    bool ok=sample.converged(bddcOptions.relativeTolerance);
    if (poissonMpiRankZero(bddcOptions.mpi))
      std::cout << "Poisson result: status=" << (ok ? "converged" : "not-converged")
                << " solver=direct calls=0 rel=" << sample.relative << '\n';
    return ok ? 0 : 2;
  }

  bool runConverged=false;
  auto run = [&]<class Transfer>()
  {
    auto started=std::chrono::steady_clock::now();
    using BddcSubdomain = Subdomain<1,double,double,Transfer>;
    auto subs = constructPoissonSubdomains<BddcSubdomain>(matrices,ifa,bddcOptions);
    for (int i=0; i<N; ++i)
      if (subs.has(i))
        subs[i].setRhs(subdomains[i].rhs());

    timer.start("bddc creation");
    std::vector<int> activeIds(N);
    std::iota(activeIds.begin(),activeIds.end(),0);
    BDDCSolver<BddcSubdomain,decltype(subs)> bddcSolver(subs,ifa.coarseConstraints(),
      activeIds,bddcOptions.useCg,bddcOptions.verbose,bddcOptions.mpi);
    bddcSolver.setRestrictionPasses(bddcOptions.restrictionPasses);
    timer.stop("bddc creation");
    if (histogram && poissonMpiRankZero(bddcOptions.mpi))
      std::filesystem::create_directories(prefix);
    bddcSolver.setHistogramOutput(histogram ? prefix : "");

    std::vector<double> resNorm;
    bool solverStopped = false;
    PoissonValidation::Sample sample;
    std::ofstream history;
    int accepted=0;
    if (monitor)
    {
      sample=monitor->evaluate(subs);
      runConverged=sample.converged(bddcOptions.relativeTolerance);
      printPoissonResidual(sample,0,bddcOptions.mpi);
      if (poissonMpiRankZero(bddcOptions.mpi))
      {
        std::filesystem::create_directories(prefix);
        history.open(prefix+"/convergence.csv");
        if (!history) throw std::runtime_error("Cannot write Poisson convergence.csv");
        history << "solver_calls,accepted_updates,absolute_residual,rhs_norm,relative_residual,interface_jump\n"
                << std::setprecision(17)
                << "0,0," << sample.absolute << ',' << sample.rhsNorm << ','
                << sample.relative << ',' << sample.relativeJump << '\n';
      }
    }
    for (int k=0; k<iter && !runConverged; ++k)
    {
      timer.start("BDDC solve");
      resNorm.push_back(bddcSolver.solve());
      timer.stop("BDDC solve");
      if (!bddcSolver.stopped()) ++accepted;
      if (monitor)
      {
        sample=monitor->evaluate(subs);
        runConverged=sample.converged(bddcOptions.relativeTolerance);
        printPoissonResidual(sample,k+1,bddcOptions.mpi);
        if (history.is_open())
          history << k+1 << ',' << accepted << ',' << sample.absolute << ','
                  << sample.rhsNorm << ',' << sample.relative << ',' << sample.relativeJump << '\n';
      }
      if (bddcSolver.stopped())
      {
        solverStopped = true;
        if (poissonMpiRankZero(bddcOptions.mpi))
        {
          std::cout << "BDDC stopped before applying the current correction: "
                    << (bddcSolver.stagnated() ? "stagnation" : "breakdown") << ".\n";
        }
        break;
      }
      if (vtk)
      {
        ScopedTimingSection out("output");
        std::vector<std::string> files, filesCorr, filesRaw, filesRes, filesRestricted;
        for (int i=0; i<N; ++i)
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
      }
      if (!monitor)
        runConverged=std::isfinite(resNorm.back()) && resNorm.back()<=tolerance;
    }
    if (poissonMpiRankZero(bddcOptions.mpi))
    {
      if (!solverStopped && resNorm.size() >= 2)
      {
        int lookback = std::min(10,static_cast<int>(resNorm.size())-1);
        double older = resNorm[resNorm.size()-lookback-1];
        if (older > 0)
        {
          double contraction = std::pow(resNorm.back()/older,1.0/lookback);
          std::cout << "Estimated contraction factor: " << contraction << ". (kappa ~ " << (1+contraction)/(1-contraction) << ").\n";
        }
      }
      auto [nIter, nTransfer, iB, rB, pB, cB] = bddcSolver.traffic();
      if (nIter > 0)
        std::cout << "average amount of data exchanged per iteration: init=" << iB/(nIter*1024.0) << "kb restriction=" << rB/(nIter*1024.0) << "kb prolongation=" << pB/(nIter*1024.0) << "kb coarseGrid=" << cB/(nIter*1024.0) << "kb\n";
    }
    if (bddcOptions.compression)
      printPoissonCompressionReport(bddcSolver.compressionTraffic(),bddcOptions.mpi);
    else
      printPoissonUncompressedReport(bddcSolver.traffic(),bddcOptions.mpi);
    if (monitor) savePoissonSolution(sample,prefix,bddcOptions.mpi);
    double seconds=poissonMax(std::chrono::duration<double>(
      std::chrono::steady_clock::now()-started).count(),bddcOptions.mpi);
    if (poissonMpiRankZero(bddcOptions.mpi))
      std::cout << "Poisson result: status="
                << (runConverged ? "converged" : bddcSolver.stopped() ?
                    (bddcSolver.stagnated() ? "stagnation" : "breakdown") : "iteration-limit")
                << " solver=" << (bddcOptions.useCg ? "cg" : "sd")
                << " calls=" << resNorm.size() << " accepted_updates=" << accepted
                << " rel=" << (monitor ? sample.relative : std::numeric_limits<double>::quiet_NaN())
                << "\nPoisson solver wall seconds: " << seconds << '\n';
    return resNorm;
  };

  std::vector<double> resNorm;
  if (bddcOptions.compression)
    resNorm = run.template operator()<SpaceTransferDataCompression<1,double,double,std::uint32_t,std::uint8_t>>();
  else
    resNorm = run.template operator()<SpaceTransfer<1,double,double>>();

  timer.stop("BDDC");

  if (timing && poissonMpiRankZero(bddcOptions.mpi))
    std::cout << timer;

  return runConverged ? 0 : 2;
}
