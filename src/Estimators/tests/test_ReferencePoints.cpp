//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2023 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//////////////////////////////////////////////////////////////////////////////////////

#include "catch.hpp"

#include "NEReferencePoints.h"
#include "ReferencePointsInput.h"
#include "ValidReferencePointsInput.h"
#include "OhmmsData/Libxml2Doc.h"
#include "EstimatorTesting.h"
#include "Particle/tests/MinimalParticlePool.h"

/** \file
 *  This is a postfacto unit testing written for reference points during porting of EnergyDensity
 *  to the batched version of the estimators.
 */

namespace qmcplusplus
{
constexpr bool generate_test_data = false;

template<typename T1, typename T2, unsigned D>
bool approxEquality(const TinyVector<T1, D>& val_a, const TinyVector<T2, D>& val_b)
{
  for (int i = 0; i < D; ++i)
    if (val_a[i] != Approx(val_b[i]))
      return false;
  return true;
}
  
TEST_CASE("ReferencePoints::DefaultConstruction", "[estimators]")
{
  using Input = testing::ValidReferencePointsInputs;
  Libxml2Document doc;
  bool okay       = doc.parseFromString(Input::xml[Input::valid::CELL]);
  xmlNodePtr node = doc.getRoot();
  ReferencePointsInput rpi(node);

  auto lattice = testing::makeTestLattice();
  Communicate* comm;
  comm = OHMMS::Controller;
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
  NEReferencePoints ref_points(rpi, pset, ref_psets);

  if (generate_test_data)
  {
    testing::TestableNEReferencePoints tref_points(ref_points);
    std::cout << "expected_reference_points" << tref_points;
  }
  typename NEReferencePoints::Points expected_reference_points{
      {"a1", {3.37316107749939, 3.37316107749939, 0}},
      {"a2", {0, 3.37316107749939, 3.37316107749939}},
      {"a3", {3.37316107749939, 0, 3.37316107749939}},
      {"cmmm", {-3.37316107749939, -3.37316107749939, -3.37316107749939}},
      {"cmmp", {0, -3.37316107749939, 0}},
      {"cmpm", {-3.37316107749939, 0, 0}},
      {"cmpp", {0, 0, 3.37316107749939}},
      {"cpmm", {0, 0, -3.37316107749939}},
      {"cpmp", {3.37316107749939, 0, 0}},
      {"cppm", {0, 3.37316107749939, 0}},
      {"cppp", {3.37316107749939, 3.37316107749939, 3.37316107749939}},
      {"f1m", {-1.686580538749695, -1.686580538749695, 0}},
      {"f1p", {1.686580538749695, 1.686580538749695, 0}},
      {"f2m", {0, -1.686580538749695, -1.686580538749695}},
      {"f2p", {0, 1.686580538749695, 1.686580538749695}},
      {"f3m", {-1.686580538749695, 0, -1.686580538749695}},
      {"f3p", {1.686580538749695, 0, 1.686580538749695}},
      {"ion1", {0, 0, 0}},
      {"ion2", {1.686580538749695, 1.686580538749695, 1.686580538749695}},
      {"r1", {3.37316107749939, 3.37316107749939, 0}},
      {"r2", {0, 3.37316107749939, 3.37316107749939}},
      {"r3", {3.37316107749939, 0, 3.37316107749939}},
      {"zero", {0, 0, 0}},
  };

  for (auto& [key, value] : ref_points.get_points())
  {
    bool coords_match = approxEquality(expected_reference_points[key], value);
    CHECK(coords_match);
  }
}

TEST_CASE("ReferencePoints::Construction", "[estimators]")
{
  using Input = testing::ValidReferencePointsInputs;
  Libxml2Document doc;
  bool okay       = doc.parseFromString(Input::xml[Input::valid::CELL]);
  xmlNodePtr node = doc.getRoot();
  ReferencePointsInput rpi(node);

  auto lattice = testing::makeTestLattice();
  Communicate* comm;
  comm = OHMMS::Controller;
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

  if (generate_test_data)
  {
    testing::TestableNEReferencePoints tref_points(ref_points);
    std::cout << "expected_reference_points" << tref_points;
  }
  typename NEReferencePoints::Points expected_reference_points{
      {"a1", {3.37316107749939, 3.37316107749939, 0}},
      {"a2", {0, 3.37316107749939, 3.37316107749939}},
      {"a3", {3.37316107749939, 0, 3.37316107749939}},
      {"cmmm", {-3.37316107749939, -3.37316107749939, -3.37316107749939}},
      {"cmmp", {0, -3.37316107749939, 0}},
      {"cmpm", {-3.37316107749939, 0, 0}},
      {"cmpp", {0, 0, 3.37316107749939}},
      {"cpmm", {0, 0, -3.37316107749939}},
      {"cpmp", {3.37316107749939, 0, 0}},
      {"cppm", {0, 3.37316107749939, 0}},
      {"cppp", {3.37316107749939, 3.37316107749939, 3.37316107749939}},
      {"f1m", {-1.686580538749695, -1.686580538749695, 0}},
      {"f1p", {1.686580538749695, 1.686580538749695, 0}},
      {"f2m", {0, -1.686580538749695, -1.686580538749695}},
      {"f2p", {0, 1.686580538749695, 1.686580538749695}},
      {"f3m", {-1.686580538749695, 0, -1.686580538749695}},
      {"f3p", {1.686580538749695, 0, 1.686580538749695}},
      {"ion1", {0, 0, 0}},
      {"ion2", {1.686580538749695, 1.686580538749695, 1.686580538749695}},
      {"r1", {3.37316107749939, 3.37316107749939, 0}},
      {"r2", {0, 3.37316107749939, 3.37316107749939}},
      {"r3", {3.37316107749939, 0, 3.37316107749939}},
      {"zero", {0, 0, 0}},
  };

  for (auto& [key, value] : ref_points.get_points())
  {
    bool coords_match = approxEquality(expected_reference_points[key], value);
    CHECK(coords_match);
  }
}


} // namespace qmcplusplus
