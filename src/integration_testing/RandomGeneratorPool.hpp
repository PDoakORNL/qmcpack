//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2025 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//////////////////////////////////////////////////////////////////////////////////////

/** \file
 *  For many integration tests we need a random generator pool that is
 *  not static and global so we are not subject to sequencing
 *  indeterminancy from the order test are run in and other
 *  nightmares.
 *  IMHO RandomNumberControl should be replaced by something like this
 *  in production as well.  With the Pool owned at the QMCMain level
 *  instead of by the C++ runtime.
 */

#include "type_traits/template_types.hpp"
#include "Configuration.h"
#include <RandomBase.h>
#include <StdRandom.h>

namespace qmcplusplus::testing
{

template<class RNGGEN>
class RandomNumberGeneratorPool
{
public:
  using FullPrecReal = QMCTraits::FullPrecRealType;

  RandomNumberGeneratorPool(int num);
  RandomNumberGeneratorPool(int num, typename RNGGEN::uint_type, unsigned long long discard = 0x10000ULL);
  ~RandomNumberGeneratorPool();
  RandomNumberGeneratorPool(const RandomNumberGeneratorPool&)            = delete;
  RandomNumberGeneratorPool(RandomNumberGeneratorPool&&)                 = delete;
  RandomNumberGeneratorPool& operator=(const RandomNumberGeneratorPool&) = delete;
  RandomNumberGeneratorPool& operator=(RandomNumberGeneratorPool&&)      = delete;

  RefVector<RNGGEN> getRngRefs();

private:
  std::vector<RNGGEN> rng_pool;
};

extern template class RandomNumberGeneratorPool<StdRandom<double>>;

} // namespace qmcplusplus::testing
