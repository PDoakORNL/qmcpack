//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2023 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//////////////////////////////////////////////////////////////////////////////////////

#include "catch.hpp"

#include "NESpaceGrid.h"
#include "SpaceGridInput.h"
#include "NEReferencePoints.h"
#include "ReferencePointsInput.h"
#include "ValidReferencePointsInput.h"
#include "ValidSpaceGridInput.h"
#include "OhmmsData/Libxml2Doc.h"
#include "EstimatorTesting.h"
#include "Particle/tests/MinimalParticlePool.h"

/** \file
 *  This is a postfacto unit testing written for NESpaceGrid during porting of EnergyDensity
 *  to the batched version of the estimators.
 */

namespace qmcplusplus
{
constexpr bool generate_test_data = false;

TEST_CASE("ReferencePoints::Construction", "[estimators]")
{
  Communicate* comm;
  comm = OHMMS::Controller;

  using Input = testing::ValidSpaceGridInput;
  Libxml2Document doc;
  bool okay       = doc.parseFromString(Input::xml[Input::valid::ORIGIN]);
  xmlNodePtr node = doc.getRoot();
  SpaceGridInput sgi(node);

  ReferencePointsInput rpi;

  auto particle_pool = MinimalParticlePool::make_diamondC_1x1x1(comm);
  auto& pset         = *(particle_pool.getParticleSet("e"));
  auto& pset_ions    = *(particle_pool.getParticleSet("ion"));
  // Setup particleset
  pset.R = ParticleSet::ParticlePos{{1.751870349, 4.381521229, 2.865202269}, {3.244515371, 4.382273176, 4.21105285},
                                    {3.000459944, 3.329603408, 4.265030556}, {3.748660329, 3.63420622, 5.393637791},
                                    {3.033228526, 3.391869137, 4.654413566}, {3.114198787, 2.654334594, 5.231075822},
                                    {3.657151589, 4.883870516, 4.201243939}, {2.97317591, 4.245644974, 4.284564732}};


  RefVector<ParticleSet> ref_psets;
  ref_psets.push_back(pset_ions);
  NEReferencePoints ref_points(std::move(rpi), pset, ref_psets);

  // EnergyDensityEstimator gets this from an enum giving indexes into each offset of AOS buffer.
  // It is a smell.
  NESpaceGrid space_grid(sgi, ref_points.get_points(), 3 , false);
  
}

}
