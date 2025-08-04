//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2025 QMCPACK developers.
//
// File developed by: Peter W. Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_PARTICLESETEXTERN_HPP
#define QMCPLUSPLUS_PARTICLESETEXTERN_HPP

#include "ParticleSet.h"
#include <RandomBase.h>
#include <StdRandom.h>
#include <RandomGenerator.h>

namespace qmcplusplus
{

extern template void ParticleSet::randomizeFromSourceWithEngine<RandomBase<double>>(const ParticleSet& pset_source,
                                                                                    RandomBase<double>& rng);

extern template void ParticleSet::randomizeFromSourceWithEngine<RNGThreadSafe<StdRandom<double>>>(
    const ParticleSet& pset_source,
    RNGThreadSafe<StdRandom<double>>& rng);

} // namespace qmcplusplus
#endif
