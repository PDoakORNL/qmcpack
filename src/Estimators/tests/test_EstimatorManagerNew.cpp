//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File refactored from: Refactored from test_manager.cpp
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "Message/Communicate.h"
#include "OhmmsData/Libxml2Doc.h"
#include "QMCHamiltonians/QMCHamiltonian.h"
#include "Particle/tests/MinimalParticlePool.h"
#include "QMCWaveFunctions/tests/MinimalWaveFunctionPool.h"
#include "Estimators/EstimatorManagerInput.h"
#include "SpinDensityInput.h"
#include "MomentumDistributionInput.h"
#include "OneBodyDensityMatricesInput.h"
#include "test_EstimatorManagerInput.h"
#include "Estimators/EstimatorManagerNew.h"
#include "Estimators/ScalarEstimatorBase.h"
#include "Estimators/tests/EstimatorManagerNewTest.h"

#include <stdio.h>
#include <sstream>

namespace qmcplusplus
{

TEST_CASE("EstimatorManagerNew", "[estimators]")
{
  Communicate* comm = OHMMS::Controller;

  testing::EstimatorManagerNewTest embt(comm, 1);

  REQUIRE(embt.em.size() == 0);

  REQUIRE(embt.testAddGetEstimator());
}

TEST_CASE("EstimatorManagerNew::ConstructFromEstimatorManagerInput", "[estimators]")
{
  Communicate* comm = OHMMS::Controller;
  using namespace testing;
  Libxml2Document estimators_doc = createEstimatorManagerNewInputXML();
  EstimatorManagerInput emi(estimators_doc.getRoot());

  MinimalParticlePool mpp;
  ParticleSetPool particle_pool = mpp(comm);
  MinimalWaveFunctionPool wfp;
  WaveFunctionPool wavefunction_pool = wfp(comm, particle_pool);
  wavefunction_pool.setPrimary(wavefunction_pool.getWaveFunction("psi0"));
  auto& pset_target = *(particle_pool.getParticleSet("e"));
  auto& species_set = pset_target.getSpeciesSet();
  auto& wf_factory  = *(wavefunction_pool.getWaveFunctionFactory("wavefunction"));

  EstimatorManagerNew emn(comm, std::move(emi), pset_target, *(wavefunction_pool.getPrimary()), wf_factory);
}

TEST_CASE("EstimatorManagerNew::collectScalarEstimators", "[estimators]")
{
  Communicate* comm = OHMMS::Controller;

  testing::EstimatorManagerNewTest embt(comm, 1);
  // by design we have done no averaging here
  // the division by total weight happens only when a block is over and the
  // accumulated data has been reduced down.  So here there should just be simple sums.

  embt.fakeSomeScalarSamples();
  embt.collectScalarEstimators();
  double correct_value = 5.0;  
  REQUIRE(embt.em.get_AverageCache()[0] == Approx(correct_value));
  correct_value = 8.0;  
  REQUIRE(embt.em.get_AverageCache()[1] == Approx(correct_value));
  correct_value = 11.0;  
  REQUIRE(embt.em.get_AverageCache()[2] == Approx(correct_value));
  correct_value = 14.0;  
  REQUIRE(embt.em.get_AverageCache()[3] == Approx(correct_value));

}


TEST_CASE("EstimatorManagerNew::collectOperatorEstimators", "[estimators]")
{
  Communicate* comm = OHMMS::Controller;

  testing::EstimatorManagerNewTest embt(comm, 1);
  // by design we have done no averaging here
  // the division by total weight happens only when a block is over and the
  // accumulated data has been reduced down.  So here there should just be simple sums.

}


TEST_CASE("EstimatorManagerNew adhoc addVector operator", "[estimators]")
{
  int num_scalars = 3;
  std::vector<double> vec_a{1.0, 2.0, 3.0};
  std::vector<double> vec_b{2.0, 3.0, 4.0};
  std::vector<std::vector<double>> est{vec_a, vec_b};
  auto addVectors = [](const auto& vec_a, const auto& vec_b) {
    std::vector<ScalarEstimatorBase::RealType> result_vector(vec_a.size(), 0.0);
    for (int i = 0; i < vec_a.size(); ++i)
      result_vector[i] = vec_a[i] + vec_b[i];
    return result_vector;
  };

  std::vector<double> reduced_scalars(num_scalars);
  reduced_scalars = std::accumulate(est.begin(), est.end(), std::vector<double>(num_scalars, 0.0), addVectors);
  std::vector<double> correct{3.0, 5.0, 7.0};
  REQUIRE(reduced_scalars == correct);
}

} // namespace qmcplusplus
