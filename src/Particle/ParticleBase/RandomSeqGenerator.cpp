//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2025 QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////

#include "RandomSeqGenerator.h"
#include "OhmmsPETE/TinyVector.h"
#include "config.h"
#include <RandomGenerator.h>

namespace qmcplusplus
{
template void makeGaussRandomWithEngine<OHMMS_PRECISION_FULL, OHMMS_DIM, StdRandom<OHMMS_PRECISION_FULL>>(
    std::vector<TinyVector<OHMMS_PRECISION_FULL, OHMMS_DIM>>& a,
    StdRandom<OHMMS_PRECISION_FULL>& rng);
}
