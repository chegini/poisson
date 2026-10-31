#ifndef POISSON_VALIDATION_HH
#define POISSON_VALIDATION_HH

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <numeric>
#include <stdexcept>
#include <utility>
#include <vector>

namespace PoissonValidation
{
  using Mapping = std::vector<std::vector<std::size_t>>;
  using Sum = std::function<void(std::vector<double>&)>;

  // Identify shared copies, retaining a distinct global ID for every interior DOF.
  // The input contains topological equivalences, not averaging weights.
  template <class SharedDofs>
  Mapping makeMapping(std::vector<int> const& sizes, SharedDofs const& shared)
  {
    std::vector<std::size_t> offsets(sizes.size()+1,0);
    for (std::size_t s=0; s<sizes.size(); ++s)
      offsets[s+1] = offsets[s]+sizes[s];
    std::vector<std::size_t> parent(offsets.back());
    std::iota(parent.begin(),parent.end(),0);
    auto root = [&](std::size_t i)
    {
      while (parent[i]!=i)
      {
        parent[i]=parent[parent[i]];
        i=parent[i];
      }
      return i;
    };
    for (auto const& copies: shared)
    {
      if (copies.empty()) continue;
      auto const& first=copies.front();
      auto a=root(offsets.at(first.s)+first.i);
      for (auto const& copy: copies)
        parent[root(offsets.at(copy.s)+copy.i)]=a;
    }
    Mapping result(sizes.size());
    std::vector<std::size_t> ids(parent.size(),parent.size());
    std::size_t next=0;
    for (std::size_t s=0; s<sizes.size(); ++s)
      for (int j=0; j<sizes[s]; ++j)
      {
        auto k=root(offsets[s]+j);
        if (ids[k]==parent.size()) ids[k]=next++;
        result[s].push_back(ids[k]);
      }
    return result;
  }

  inline double norm(std::vector<double> const& values)
  {
    double result=0;
    for (double v: values) result=std::hypot(result,v);
    return result;
  }

  struct Sample
  {
    double absolute=0, rhsNorm=0, relative=0, relativeJump=0;
    std::vector<double> solution;
    bool converged(double tolerance) const
    {
      // Do not let projection hide discontinuous subdomain solution copies.
      return std::isfinite(relative) && relative<=tolerance &&
             std::isfinite(relativeJump) && relativeJump<=tolerance;
    }
  };

  /** True residual of the assembled, continuous Poisson system.
   * Shared solution copies are averaged only for this read-only diagnostic.
   * Local matrix/RHS contributions are summed once per owner, then globally
   * reduced; interface RHS entries must not be counted as separate global DOFs.
   * This adds uncompressed diagnostic reductions and operator applications.
   */
  template <class Matrix, class Vector>
  class Monitor
  {
  public:
    Monitor(std::vector<Matrix const*> matrices_, std::vector<Vector const*> rhs_,
            Mapping mapping_, std::vector<bool> owned_, Sum sum_)
      : matrices(std::move(matrices_)), rhs(std::move(rhs_)),
        mapping(std::move(mapping_)), owned(std::move(owned_)), sum(std::move(sum_))
    {
      std::size_t size=0;
      for (auto const& ids: mapping)
        for (auto id: ids) size=std::max(size,id+1);
      multiplicity.assign(size,0);
      b.assign(size,0);
      for (std::size_t s=0; s<mapping.size(); ++s)
        for (std::size_t j=0; j<mapping[s].size(); ++j)
        {
          ++multiplicity[mapping[s][j]];
          if (owned[s]) b[mapping[s][j]]+=(*rhs[s])[j][0];
        }
      if (std::find(multiplicity.begin(),multiplicity.end(),0)!=multiplicity.end())
        throw std::invalid_argument("Poisson global DOF IDs must be contiguous");
      sum(b);
    }

    Sample evaluateGlobal(std::vector<double> x) const
    {
      std::vector<double> ax(b.size(),0);
      for (std::size_t s=0; s<mapping.size(); ++s)
        if (owned[s])
        {
          Vector local(mapping[s].size()), product(mapping[s].size());
          for (std::size_t j=0; j<mapping[s].size(); ++j)
            local[j][0]=x[mapping[s][j]];
          matrices[s]->mv(local,product);
          for (std::size_t j=0; j<mapping[s].size(); ++j)
            ax[mapping[s][j]]+=product[j][0];
        }
      sum(ax);
      for (std::size_t i=0; i<b.size(); ++i) ax[i]=b[i]-ax[i];
      Sample result;
      result.absolute=norm(ax);
      result.rhsNorm=norm(b);
      result.relative=result.rhsNorm>0 ? result.absolute/result.rhsNorm : result.absolute;
      result.solution=std::move(x);
      return result;
    }

    template <class Storage>
    Sample evaluate(Storage const& storage) const
    {
      std::vector<double> x(b.size(),0);
      for (std::size_t s=0; s<mapping.size(); ++s)
        if (owned[s])
        {
          auto local=storage[s].getSolution();
          for (std::size_t j=0; j<mapping[s].size(); ++j) x[mapping[s][j]]+=local[j][0];
        }
      sum(x);
      for (std::size_t i=0; i<x.size(); ++i) x[i]/=multiplicity[i];
      std::vector<double> defect(1,0);
      for (std::size_t s=0; s<mapping.size(); ++s)
        if (owned[s])
        {
          auto local=storage[s].getSolution();
          for (std::size_t j=0; j<mapping[s].size(); ++j)
          {
            double d=local[j][0]-x[mapping[s][j]];
            defect[0]+=d*d;
          }
        }
      sum(defect);
      auto result=evaluateGlobal(std::move(x));
      result.relativeJump=std::sqrt(defect[0])/std::max(1.0,norm(result.solution));
      return result;
    }

    std::vector<double> const& globalRhs() const { return b; }

  private:
    std::vector<Matrix const*> matrices;
    std::vector<Vector const*> rhs;
    Mapping mapping;
    std::vector<bool> owned;
    Sum sum;
    std::vector<double> multiplicity, b;
  };
}

#endif
