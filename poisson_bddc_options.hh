#ifndef POISSON_BDDC_OPTIONS_HH
#define POISSON_BDDC_OPTIONS_HH

#ifdef KASKADE_HAVE_MPI
#include <mpi.h>
#endif

#include "mg/bddc.hpp"
#include "utilities/kaskopt.hh"

#include <cstdint>
#include <cstddef>
#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <type_traits>
#include <vector>

struct PoissonBddcOptions
{
  bool mpi = false;
  bool compression = false;
  int compressionBits = 16;
  int compressionMinInterfaceSize = 0;
  int compressionThresholdInterfaceTypes = Kaskade::BDDC::FACE;
  bool graphLifting = true;
  bool graphLiftingAfterQuantization = false;
  bool excludeCornersFromCompression = false;
  bool huffman = true;
  bool bitlength = true;
  bool zigzag = true;
};

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
              << "  encode time (rank-summed): " << traffic.encodeTimeMs << " ms\n"
              << "  decode time (rank-summed): " << traffic.decodeTimeMs << " ms\n";
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
        if (options.graphLifting)
          storage[id].transfer().setGraphLiftingFromMatrix(matrix);
      configurePoissonTransfer(storage[id].transfer(),options);
    }
  return storage;
}

#endif
