#ifndef POISSON_BDDC_OPTIONS_HH
#define POISSON_BDDC_OPTIONS_HH

#ifdef KASKADE_HAVE_MPI
#include <mpi.h>
#endif

#include "mg/bddc.hpp"
#include "linalg/direct.hh"
#include "utilities/kaskopt.hh"
#include "poisson_validation.hh"

#include <cstdint>
#include <cstddef>
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

struct PoissonBddcOptions
{
  bool mpi = false;
  bool compression = false;
  int compressionBits = 16;
  int restrictCompressionBits = -1;
  int prolongateCompressionBits = -1;
  int restrictionPasses = 1;
  bool useCg = false;
  bool verbose = true;
  double relativeTolerance = 0;
  int compressionMinInterfaceSize = 0;
  int compressionThresholdInterfaceTypes = Kaskade::BDDC::FACE;
  bool graphLifting = false;
  bool graphLiftingAfterQuantization = true;
  bool excludeCornersFromCompression = false;
  bool huffman = true;
  bool bitlength = true;
  bool zigzag = true;
};

inline void validatePoissonOptions(PoissonBddcOptions const& options)
{
  auto valid=[](int bits) { return bits>=1 && bits<=32; };
  if (!valid(options.compressionBits) ||
      (options.restrictCompressionBits!=-1 && !valid(options.restrictCompressionBits)) ||
      (options.prolongateCompressionBits!=-1 && !valid(options.prolongateCompressionBits)))
    throw std::invalid_argument("Compression bits must be 1..32; directional overrides may be -1 (inherit). 48/64-bit quantization is not supported.");
  if (options.restrictionPasses<1 || options.restrictionPasses>16 ||
      (options.restrictionPasses!=1 && !options.compression))
    throw std::invalid_argument("restrictionPasses must be 1..16; extra passes require --compression 1");
  if (!std::isfinite(options.relativeTolerance) || options.relativeTolerance<0)
    throw std::invalid_argument("relativeTolerance must be finite and nonnegative; 0 selects legacy stopping");
  // The legacy pre-quantization transform requantizes each neighbor separately,
  // while the published-residual accumulator tracks nodal quantization instead.
  // Feedback must track exactly what receivers reconstruct. Integer lifting
  // AFTER quantization preserves those symbols and their reconstruction.
  if (options.compression && options.graphLifting && !options.graphLiftingAfterQuantization)
    throw std::invalid_argument("Graph lifting before quantization is not supported by the validated feedback path; use --graphLifting 0 --graphLiftingAfterQuantization 1, or disable both transforms");
}

inline void configurePoissonTransfer(auto& transfer, PoissonBddcOptions const& options)
{
  if constexpr (requires {
                  transfer.setQuantizationBits(0);
                  transfer.setRestrictEncoding(true);
                  transfer.setProlongateEncoding(true);
                  transfer.setRestrictTransform(true);
                  transfer.setProlongateTransform(true);
                  transfer.setRestrictTransformAfterQuantization(true);
                  transfer.setProlongateTransformAfterQuantization(true);
                  transfer.setRestrictBitlengthEncoding(true);
                  transfer.setProlongateBitlengthEncoding(true);
                  transfer.setRestrictZigzagEncoding(true);
                  transfer.setProlongateZigzagEncoding(true);
                })
  {
    transfer.setQuantizationBits(options.compressionBits);
    if (options.restrictCompressionBits!=-1)
      transfer.setRestrictQuantizationBits(options.restrictCompressionBits);
    if (options.prolongateCompressionBits!=-1)
      transfer.setProlongateQuantizationBits(options.prolongateCompressionBits);
    if constexpr (requires { transfer.setCompressionMinInterfaceSize(0); })
      transfer.setCompressionMinInterfaceSize(static_cast<std::size_t>(std::max(0, options.compressionMinInterfaceSize)));
    if constexpr (requires { transfer.setCompressionThresholdInterfaceTypes(0); })
      transfer.setCompressionThresholdInterfaceTypes(options.compressionThresholdInterfaceTypes);
    transfer.enableTransform((options.graphLifting || options.graphLiftingAfterQuantization)
                               ? Kaskade::BDDC::TransformType::GRAPH_LIFTING
                               : Kaskade::BDDC::TransformType::NONE);
    transfer.setRestrictEncoding(options.huffman);
    transfer.setProlongateEncoding(options.huffman);
    transfer.setRestrictTransform(options.graphLifting);
    transfer.setProlongateTransform(options.graphLifting);
    transfer.setRestrictTransformAfterQuantization(options.graphLiftingAfterQuantization);
    transfer.setProlongateTransformAfterQuantization(options.graphLiftingAfterQuantization);
    transfer.setRestrictBitlengthEncoding(options.bitlength);
    transfer.setProlongateBitlengthEncoding(options.bitlength);
    transfer.setRestrictZigzagEncoding(options.zigzag);
    transfer.setProlongateZigzagEncoding(options.zigzag);
  }
}

class PoissonMpiSession
{
public:
  PoissonMpiSession(int& argc, char**& argv, bool enabled)
  {
    if (!enabled)
      return;
#ifdef KASKADE_HAVE_MPI
    int initialized = 0;
    MPI_Initialized(&initialized);
    if (!initialized)
    {
      MPI_Init(&argc,&argv);
      ownsSession = true;
    }
#else
    throw std::runtime_error("--mpi 1 requires an MPI-enabled build");
#endif
  }

  ~PoissonMpiSession()
  {
#ifdef KASKADE_HAVE_MPI
    if (ownsSession)
      MPI_Finalize();
#endif
  }

private:
  bool ownsSession = false;
};

inline bool poissonMpiRankZero(bool enabled)
{
#ifdef KASKADE_HAVE_MPI
  if (enabled)
  {
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    return rank == 0;
  }
#else
  (void)enabled;
#endif
  return true;
}

inline void poissonSum(std::vector<double>& values, bool enabled)
{
#ifdef KASKADE_HAVE_MPI
  if (enabled)
  {
    if (values.size()>static_cast<std::size_t>(std::numeric_limits<int>::max()))
      throw std::runtime_error("Poisson diagnostic exceeds MPI count range");
    MPI_Allreduce(MPI_IN_PLACE,values.data(),static_cast<int>(values.size()),
                  MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  }
#else
  (void)values;
  (void)enabled;
#endif
}

inline void printPoissonResidual(PoissonValidation::Sample const& sample,
                                 int calls, bool mpiEnabled)
{
  if (poissonMpiRankZero(mpiEnabled))
    std::cout << std::setprecision(17)
              << "Poisson true residual: calls=" << calls
              << " abs=" << sample.absolute << " rhs=" << sample.rhsNorm
              << " rel=" << sample.relative << " interface_jump=" << sample.relativeJump << '\n';
}

inline void savePoissonSolution(PoissonValidation::Sample const& sample,
                               std::string const& prefix, bool mpiEnabled)
{
  if (!poissonMpiRankZero(mpiEnabled)) return;
  std::filesystem::create_directories(prefix);
  std::ofstream out(prefix+"/solution.csv");
  if (!out) throw std::runtime_error("Cannot write Poisson solution.csv");
  out << "global_dof,value\n" << std::setprecision(17);
  for (std::size_t i=0; i<sample.solution.size(); ++i)
    out << i << ',' << sample.solution[i] << '\n';
}

inline void printPoissonUncompressedReport(std::array<int,6> const& traffic, bool mpiEnabled)
{
  // traffic() already contains MPI-global logical counts. Do not sum it again.
  if (!poissonMpiRankZero(mpiEnabled)) return;
  auto bytes=static_cast<std::int64_t>(traffic[3])+traffic[4];
  std::cout << "Compression report (BDDC, global):\n"
            << "  raw sent bytes: " << bytes << '\n'
            << "  transmitted sent bytes: " << bytes << '\n';
}

// Optional diagnostic reference: assemble exactly the same local operators and
// solve the continuous global system on rank zero with UMFPACK. This is intended
// for small validation problems, not the scaling campaign.
template <class Matrix, class Vector>
std::vector<double> solvePoissonDirect(std::vector<Matrix const*> const& matrices,
                                      PoissonValidation::Mapping const& mapping,
                                      std::vector<double> const& rhs, bool mpiEnabled)
{
  if (rhs.size()>static_cast<std::size_t>(std::numeric_limits<int>::max()))
    throw std::runtime_error("Poisson direct reference exceeds integer index range");
  std::vector<double> solution(rhs.size(),0);
  int failed=0;
  if (poissonMpiRankZero(mpiEnabled))
    try
    {
      Kaskade::MatrixAsTriplet<double> global;
      for (std::size_t s=0; s<matrices.size(); ++s)
      {
        Kaskade::MatrixAsTriplet<double> local(*matrices[s]);
        for (std::size_t j=0; j<local.data.size(); ++j)
          global.addEntry(static_cast<int>(mapping[s][local.ridx[j]]),
                          static_cast<int>(mapping[s][local.cidx[j]]),local.data[j]);
      }
      Vector b(rhs.size()), x(rhs.size());
      for (std::size_t i=0; i<rhs.size(); ++i) b[i][0]=rhs[i];
      Kaskade::DirectSolver<Vector,Vector> inverse(global,::DirectType::UMFPACK);
      inverse.apply(x,b);
      for (std::size_t i=0; i<rhs.size(); ++i) solution[i]=x[i][0];
    }
    catch (std::exception const& error)
    {
      failed=1;
      std::cerr << "Poisson direct reference failed: " << error.what() << '\n';
    }
#ifdef KASKADE_HAVE_MPI
  if (mpiEnabled)
  {
    MPI_Bcast(&failed,1,MPI_INT,0,MPI_COMM_WORLD);
    if (!failed)
      MPI_Bcast(solution.data(),static_cast<int>(solution.size()),MPI_DOUBLE,0,MPI_COMM_WORLD);
  }
#endif
  if (failed) throw std::runtime_error("Poisson direct reference factorization failed");
  return solution;
}

inline double poissonMax(double value, bool enabled)
{
#ifdef KASKADE_HAVE_MPI
  if (enabled) MPI_Allreduce(MPI_IN_PLACE,&value,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
#else
  (void)enabled;
#endif
  return value;
}

template <class Traffic>
inline void printPoissonCompressionReport(Traffic const& local, bool mpiEnabled)
{
  auto print = [](Traffic const& traffic)
  {
    auto const savedBytes = static_cast<std::int64_t>(traffic.originalSentBytes)
                          - static_cast<std::int64_t>(traffic.sentBytes);
    auto const savedBits = 8 * savedBytes;
    auto const savings = traffic.originalSentBytes > 0
                           ? 100.0 * static_cast<double>(savedBytes)
                             / static_cast<double>(traffic.originalSentBytes)
                           : 0.0;
    auto const ratio = traffic.sentBytes > 0
                         ? static_cast<double>(traffic.originalSentBytes)
                           / static_cast<double>(traffic.sentBytes)
                         : 0.0;
    std::cout << "Compression report (BDDC, global):\n"
              << "  raw sent bytes: " << traffic.originalSentBytes << "\n"
              << "  transmitted sent bytes: " << traffic.sentBytes << "\n"
              << "  bits saved: " << savedBits << "\n"
              << "  savings: " << savings << "%\n"
              << "  ratio: " << ratio << ":1\n"
              << "  shared codebook estimate: " << traffic.codebookBits << " bits\n"
              << "  raw received bytes: " << traffic.originalReceivedBytes << "\n"
              << "  received bytes: " << traffic.receivedBytes << "\n"
              << "  codec-only timings: unavailable (transfer timers are not instrumented)\n"
              << "  byte scope: logical restriction/prolongation payloads, including intra-rank delivery; excludes diagnostic reductions, MPI framing and codebook training\n";
  };

#ifdef KASKADE_HAVE_MPI
  int initialized = 0;
  MPI_Initialized(&initialized);
  if (mpiEnabled && initialized)
  {
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    unsigned long long localBytes[4] = {
      static_cast<unsigned long long>(local.originalSentBytes),
      static_cast<unsigned long long>(local.originalReceivedBytes),
      static_cast<unsigned long long>(local.sentBytes),
      static_cast<unsigned long long>(local.receivedBytes)};
    unsigned long long globalBytes[4] = {0,0,0,0};
    MPI_Allreduce(localBytes,globalBytes,4,MPI_UNSIGNED_LONG_LONG,MPI_SUM,MPI_COMM_WORLD);

    double localTimes[2] = {local.encodeTimeMs,local.decodeTimeMs};
    double globalTimes[2] = {0.0,0.0};
    MPI_Allreduce(localTimes,globalTimes,2,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

    unsigned long long localCodebook = static_cast<unsigned long long>(local.codebookBits);
    unsigned long long globalCodebook = 0;
    MPI_Allreduce(&localCodebook,&globalCodebook,1,MPI_UNSIGNED_LONG_LONG,MPI_MAX,MPI_COMM_WORLD);

    if (rank != 0)
      return;

    Traffic global = local;
    global.originalSentBytes = globalBytes[0];
    global.originalReceivedBytes = globalBytes[1];
    global.sentBytes = globalBytes[2];
    global.receivedBytes = globalBytes[3];
    global.encodeTimeMs = globalTimes[0];
    global.decodeTimeMs = globalTimes[1];
    global.codebookBits = globalCodebook;
    print(global);
    return;
  }
#else
  (void)mpiEnabled;
#endif
  print(local);
}

template <class Subdomain, class MatrixPointer, class Interfaces>
auto constructPoissonSubdomains(std::vector<MatrixPointer> const& matrixPointers,
                                Interfaces const& interfaces,
                                PoissonBddcOptions const& options)
{
  Kaskade::BDDC::MpiDomainLayout layout;
  layout.configure(options.mpi,static_cast<int>(matrixPointers.size()));
  Kaskade::BDDC::OwnedSubdomainStorage<Subdomain> storage(matrixPointers.size());
  for (size_t id = 0; id < matrixPointers.size(); ++id)
    if (layout.owns(static_cast<int>(id)))
    {
      auto const& matrix = *matrixPointers[id];
      storage.emplace(id,static_cast<int>(id),matrix,interfaces,
                      options.excludeCornersFromCompression);
      if constexpr (requires(Subdomain& subdomain, decltype(matrix) matrixObject) {
                      subdomain.transfer().setGraphLiftingFromMatrix(matrixObject);
                    })
        if (options.graphLifting || options.graphLiftingAfterQuantization)
          storage[id].transfer().setGraphLiftingFromMatrix(matrix);
      configurePoissonTransfer(storage[id].transfer(),options);
    }
  return storage;
}

#endif
