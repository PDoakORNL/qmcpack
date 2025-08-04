//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2025 QMCPACK developers.
//
// File developed by: Peter W. Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include "RandomGeneratorPool.hpp"
#include <RandomBase.h>
#include <StdRandom.h>

namespace qmcplusplus::testing
{

template<class RNGGEN>
RandomNumberGeneratorPool<RNGGEN>::RandomNumberGeneratorPool(const int num) : rng_pool(num)
{}

template<class RNGGEN>
RandomNumberGeneratorPool<RNGGEN>::RandomNumberGeneratorPool(const int num,
                                                             typename RNGGEN::uint_type seed,
                                                             unsigned long long discard)
    : rng_pool(num, seed)
{
  int i = 0;
  for (auto& rng : rng_pool)
  {
    rng.discard(i * discard);
    ++i;
  }
}

template<class RNGGEN>
RandomNumberGeneratorPool<RNGGEN>::~RandomNumberGeneratorPool() = default;

template<class RNGGEN>
RefVector<RNGGEN> RandomNumberGeneratorPool<RNGGEN>::getRngRefs()
{
  RefVector<RNGGEN> rng_refs;
  rng_refs.reserve(rng_pool.size());
  for (auto& rng : rng_pool)
    rng_refs.push_back(rng);
  return rng_refs;
}

template class RandomNumberGeneratorPool<StdRandom<double>>;

} // namespace qmcplusplus::testing
