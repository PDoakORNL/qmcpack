//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Peter W. Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// Refactored from:  QMCHamiltonians/ReferencePoints.h
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_REFERENCE_POINTS_H
#define QMCPLUSPLUS_REFERENCE_POINTS_H

#include <Configuration.h>
#include "OhmmsData/OhmmsElementBase.h"
#include "OhmmsPETE/TinyVector.h"
#include "Particle/ParticleSet.h"
#include "QMCHamiltonians/ObservableHelper.h"
#include "OhmmsPETE/Tensor.h"
#include "ReferencePointsInput.h"

namespace qmcplusplus
{

/** Class representing input points used with SpaceGrid.s
 *  Currently these areused only by EnergyDensityEstimators but general purpose seeming
 */
class ReferencePoints
{
public:
  using Real = QMCTraits::RealType;
  using Points = std::map<std::string, TinyVector<Real, OHMMS_DIM>>;

  ReferencePoints(ReferencePointsInput&& input, ParticleSet& pset, RefVector<ParticleSet> prefs);

  void write_description(std::ostream& os, std::string& indent);
  void save(std::vector<ObservableHelper>& h5desc, hid_t gid) const;
private:
  Points points_;
  Tensor<Real, OHMMS_DIM> axes_;
};


} // namespace qmcplusplus


#endif
