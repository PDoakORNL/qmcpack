//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2023 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "EnergyDensityEstimator.h"
#include "EstimatorTesting.h"
#include "ValidEnergyDensityInput.h"
#include "Particle/tests/MinimalParticlePool.h"
#include <iostream>

namespace qmcplusplus
{

TEST_CASE("EnergyDensityEstimator::Constructor", "[estimators]")
{
  using Input = testing::ValidEnergyDensityInput;
  Communicate* comm;
  comm = OHMMS::Controller;

  ParticleSetPool particle_pool{MinimalParticlePool::make_diamondC_1x1x1(comm)};

  ParticleSet pset_elec{*(particle_pool.getParticleSet("e"))};
  ParticleSet pset_ions{*(particle_pool.getParticleSet("ion"))};

  pset_elec.R =
      ParticleSet::ParticlePos{{1.451870349, 1.381521229, 1.165202269}, {1.244515371, 1.382273176, 1.21105285},
                               {0.000459944, 1.329603408, 1.265030556}, {0.748660329, 1.63420622, 1.393637791},
                               {0.033228526, 1.391869137, 0.654413566}, {1.114198787, 1.654334594, 0.231075822},
                               {1.657151589, 0.883870516, 1.201243939}, {0.97317591, 1.245644974, 0.284564732}};

  Libxml2Document doc;
  bool okay       = doc.parseFromString(Input::xml[Input::valid::ION]);
  xmlNodePtr node = doc.getRoot();
  UPtr<EnergyDensityInput> edein;
  edein            = std::make_unique<EnergyDensityInput>(node);
  EnergyDensityEstimator e_den_est(*edein, particle_pool.getPool());
}


} // namespace qmcplusplus
