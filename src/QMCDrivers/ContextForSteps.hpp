//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2025 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File refactored from: QMCDrinverNew.h
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_CONTEXTFORSTEPS_HPP
#define QMCPLUSPLUS_CONTEXTFORSTEPS_HPP

#include "RandomBase.h"

namespace qmcplusplus
{
/// a collection of driver-specific objects needed per batch
template<typename FULLPRECREAL>
class ContextForStepsT
{
public:
  ContextForStepsT(RandomBase<FULLPRECREAL>& random_gen) : random_gen_(random_gen) {}
  RandomBase<FULLPRECREAL>& get_random_gen() { return random_gen_; }

protected:
  RandomBase<FULLPRECREAL>& random_gen_;
};

extern template class ContextForStepsT<double>;

using ContextForSteps = ContextForStepsT<double>;

} // namespace qmcplusplus
#endif
