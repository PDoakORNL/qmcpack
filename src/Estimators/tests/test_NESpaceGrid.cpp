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

namespace testing
{
template<typename T>
class NESpaceGridTests
{
public:
  static int getBufferStart(const NESpaceGrid& nesg) { return nesg.buffer_start_; }
  static int getBufferEnd(const NESpaceGrid& nesg) { return nesg.buffer_end_; }
  static auto* getOdu(const NESpaceGrid& nesg) { return nesg.odu_; }
};

template<ValidSpaceGridInput::valid VALID>
class SpaceGridEnv
{
public:
  using Input = ValidSpaceGridInput;
  SpaceGridEnv(Communicate* comm)
      : particle_pool_(MinimalParticlePool::make_diamondC_1x1x1(comm)),
        pset_elec_(*(particle_pool_.getParticleSet("e"))),
        pset_ions_(*(particle_pool_.getParticleSet("ion")))
  {
    // Setup particleset
    // particle positions must be inside the unit cell
    pset_elec_.R =
        ParticleSet::ParticlePos{{1.451870349, 1.381521229, 1.165202269}, {1.244515371, 1.382273176, 1.21105285},
                                 {0.000459944, 1.329603408, 1.265030556}, {0.748660329, 1.63420622, 1.393637791},
                                 {0.033228526, 1.391869137, 0.654413566}, {1.114198787, 1.654334594, 0.231075822},
                                 {1.657151589, 0.883870516, 1.201243939}, {0.97317591, 1.245644974, 0.284564732}};

    Libxml2Document doc;
    bool okay       = doc.parseFromString(Input::xml[VALID]);
    xmlNodePtr node = doc.getRoot();
    sgi_            = std::make_unique<SpaceGridInput>(node);
    ReferencePointsInput rpi;
    ref_psets_.push_back(pset_ions_);
    ref_points_ = std::make_unique<NEReferencePoints>(std::move(rpi), pset_elec_, ref_psets_);
  }

  UPtr<SpaceGridInput> sgi_;
  UPtr<NEReferencePoints> ref_points_;
  RefVector<ParticleSet> ref_psets_;
  ParticleSetPool particle_pool_;
  ParticleSet pset_elec_;
  ParticleSet pset_ions_;
};
} // namespace testing

constexpr bool generate_test_data = false;
using Real                        = double;

TEST_CASE("SpaceGrid::Construction", "[estimators]")
{
  using Input = testing::ValidSpaceGridInput;
  Communicate* comm;
  comm = OHMMS::Controller;

  testing::SpaceGridEnv<Input::valid::ORIGIN> sge(comm);

  // EnergyDensityEstimator gets this from an enum giving indexes into each offset of AOS buffer.
  // It is a smell.
  NESpaceGrid space_grid(*(sge.sgi_), sge.ref_points_->get_points(), 1, false);
  PooledData<Real> data_pool;
  space_grid.allocate_buffer_space(data_pool);
  using NES         = testing::NESpaceGridTests<double>;
  auto buffer_start = NES::getBufferStart(space_grid);
  auto buffer_end   = NES::getBufferEnd(space_grid);
  space_grid.write_description(std::cout, std::string(""));
  auto& sgi = *(sge.sgi_);
  auto& agr = sgi.get_axis_grids();
  for (int id = 0; id < OHMMS_DIM; ++id)
    CHECK(NES::getOdu(space_grid)[id] == agr[id].odu);

  CHECK(buffer_start == 0);
  CHECK(buffer_end == 7999);
}

TEST_CASE("SpaceGrid::Basic", "[estimators]")
{
  using Input = testing::ValidSpaceGridInput;
  Communicate* comm;
  comm = OHMMS::Controller;
  testing::SpaceGridEnv<Input::valid::ORIGIN> sge(comm);
  int num_values = 3;
  NESpaceGrid space_grid(*(sge.sgi_), sge.ref_points_->get_points(), num_values, false);
  PooledData<Real> data_pool;
  space_grid.allocate_buffer_space(data_pool);
  using NES         = testing::NESpaceGridTests<double>;
  auto buffer_start = NES::getBufferStart(space_grid);
  auto buffer_end   = NES::getBufferEnd(space_grid);
  space_grid.write_description(std::cout, std::string(""));
  auto& sgi = *(sge.sgi_);
  auto& agr = sgi.get_axis_grids();
  for (int id = 0; id < OHMMS_DIM; ++id)
    CHECK(NES::getOdu(space_grid)[id] == agr[id].odu);
  CHECK(buffer_start == 0);
  CHECK(buffer_end == 23999);

  Matrix<Real> values;
  values.resize(sge.pset_elec_.getTotalNum(), num_values);

  for (int ip = 0; ip < sge.pset_elec_.getTotalNum(); ++ip)
    for (int iv = 0; iv < num_values; ++iv)
      values(ip, iv) = ip + 0.1 * iv;

  const int ei_tid = sge.pset_elec_.addTable(sge.pset_ions_);
  sge.pset_elec_.update();
  sge.pset_ions_.update();

  std::vector<bool> p_outside(8, false);
  space_grid.accumulate(sge.pset_elec_.R, values, data_pool, p_outside, sge.pset_elec_.getDistTableAB(ei_tid));

  // check that what's in data_pool is what is expected.

  // new pset R's
  // check again
}

TEST_CASE("SpaceGrid::BadPeriodic", "[estimators]")
{
  using Input = testing::ValidSpaceGridInput;
  Communicate* comm;
  comm = OHMMS::Controller;
  testing::SpaceGridEnv<Input::valid::ORIGIN> sge(comm);
  int num_values = 3;
  NESpaceGrid space_grid(*(sge.sgi_), sge.ref_points_->get_points(), num_values, false);
  PooledData<Real> data_pool;
  space_grid.allocate_buffer_space(data_pool);
  using NES         = testing::NESpaceGridTests<double>;
  auto buffer_start = NES::getBufferStart(space_grid);
  auto buffer_end   = NES::getBufferEnd(space_grid);
  space_grid.write_description(std::cout, std::string(""));
  auto& sgi = *(sge.sgi_);
  auto& agr = sgi.get_axis_grids();
  for (int id = 0; id < OHMMS_DIM; ++id)
    CHECK(NES::getOdu(space_grid)[id] == agr[id].odu);
  CHECK(buffer_start == 0);
  CHECK(buffer_end == 23999);

  Matrix<Real> values;
  values.resize(sge.pset_elec_.getTotalNum(), num_values);

  for (int ip = 0; ip < sge.pset_elec_.getTotalNum(); ++ip)
    for (int iv = 0; iv < num_values; ++iv)
      values(ip, iv) = ip + 0.1 * iv;

  const int ei_tid = sge.pset_elec_.addTable(sge.pset_ions_);
  sge.pset_elec_.update();
  sge.pset_ions_.update();

  // set a position outside of the cell
  sge.pset_elec_.R[2] = {1.451870349, 4.381521229, 1.165202269};
  
  std::vector<bool> p_outside(8, false);
  
  CHECK_THROWS_AS(space_grid.accumulate(sge.pset_elec_.R, values, data_pool, p_outside, sge.pset_elec_.getDistTableAB(ei_tid)), std::runtime_error);
}

  
} // namespace qmcplusplus
