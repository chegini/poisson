#ifndef POISSON_BDDC_OPTIONS_HH
#define POISSON_BDDC_OPTIONS_HH

#include "mg/bddc.hpp"
#include "utilities/kaskopt.hh"

#include <cstdint>
#include <cstddef>
#include <stdexcept>
#include <type_traits>
#include <vector>

#ifdef KASKADE_HAVE_MPI
#include <mpi.h>
#endif

struct PoissonBddcOptions
{
  bool mpi = false;
  bool compression = false;
  int compressionBits = 16;
  bool graphLifting = true;
  bool huffman = true;
  bool bitlength = true;
};

inline void configurePoissonTransfer(auto& transfer, PoissonBddcOptions const& options)
{
  if constexpr (requires {
                  transfer.setQuantizationBits(0);
                  transfer.setRestrictEncoding(true);
                  transfer.setProlongateEncoding(true);
                  transfer.setRestrictTransform(true);
                  transfer.setProlongateTransform(true);
                  transfer.setRestrictBitlengthEncoding(true);
                  transfer.setProlongateBitlengthEncoding(true);
                })
  {
    transfer.setQuantizationBits(options.compressionBits);
    transfer.enableTransform(options.graphLifting
                               ? Kaskade::BDDC::TransformType::GRAPH_LIFTING
                               : Kaskade::BDDC::TransformType::NONE);
    transfer.setRestrictEncoding(options.huffman);
    transfer.setProlongateEncoding(options.huffman);
    transfer.setRestrictTransform(options.graphLifting);
    transfer.setProlongateTransform(options.graphLifting);
    transfer.setRestrictBitlengthEncoding(options.bitlength);
    transfer.setProlongateBitlengthEncoding(options.bitlength);
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
      storage.emplace(id,static_cast<int>(id),matrix,interfaces);
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
