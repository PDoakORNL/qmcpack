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
#include <iostream>
#include "EstimatorTesting.h"
#include "ValidEnergyDensityInput.h"
#include "Particle/Walker.h"
#include "Particle/tests/MinimalParticlePool.h"
#include "QMCWaveFunctions/tests/MinimalWaveFunctionPool.h"
#include "QMCHamiltonians/tests/MinimalHamiltonianPool.h"
#include "Utilities/for_testing/NativeInitializerPrint.hpp"
#include "Utilities/ProjectData.h"
#include "Utilities/ResourceCollection.h"
#include "Utilities/StdRandom.h"

constexpr bool generate_test_data = false;

namespace qmcplusplus
{

TEST_CASE("NEEnergyDensityEstimator::Constructor", "[estimators]")
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
  bool okay       = doc.parseFromString(Input::xml[Input::valid::CELL]);
  xmlNodePtr node = doc.getRoot();
  EnergyDensityInput edein{node};
  {
    NEEnergyDensityEstimator e_den_est(edein, particle_pool.getPool());
  }
}

TEST_CASE("NEEnergyDensityEstimator::AccumulateIntegration", "[estimators]")
{
  using Input = testing::ValidEnergyDensityInput;
  Communicate* comm;
  comm = OHMMS::Controller;

  ProjectData test_project("test", ProjectData::DriverVersion::BATCH);

  auto particle_pool = MinimalParticlePool::make_diamondC_1x1x1(comm);
  auto wavefunction_pool =
      MinimalWaveFunctionPool::make_diamondC_1x1x1(test_project.getRuntimeOptions(), comm, particle_pool);
  auto hamiltonian_pool = MinimalHamiltonianPool::make_hamWithEE(comm, particle_pool, wavefunction_pool);
  auto& twf             = *(wavefunction_pool.getWaveFunction("wavefunction"));
  auto& ham             = *(hamiltonian_pool.getPrimary());
  ParticleSet pset_elec{*(particle_pool.getParticleSet("e"))};
  ParticleSet pset_ions{*(particle_pool.getParticleSet("ion"))};

  pset_elec.R =
      ParticleSet::ParticlePos{{1.451870349, 1.381521229, 1.165202269}, {1.244515371, 1.382273176, 1.21105285},
                               {0.000459944, 1.329603408, 1.265030556}, {0.748660329, 1.63420622, 1.393637791},
                               {0.033228526, 1.391869137, 0.654413566}, {1.114198787, 1.654334594, 0.231075822},
                               {1.657151589, 0.883870516, 1.201243939}, {0.97317591, 1.245644974, 0.284564732}};

  auto& pset_target        = *(particle_pool.getParticleSet("e"));
  auto& trial_wavefunction = *(wavefunction_pool.getPrimary());

  Libxml2Document doc;
  bool okay       = doc.parseFromString(Input::xml[Input::valid::CELL]);
  xmlNodePtr node = doc.getRoot();
  UPtr<EnergyDensityInput> edein;
  edein = std::make_unique<EnergyDensityInput>(node);
  NEEnergyDensityEstimator e_den_est(*edein, particle_pool.getPool());

  UPtrVector<QMCHamiltonian> hams;
  UPtrVector<TrialWaveFunction> twfs;
  std::vector<ParticleSet> psets;
  hamiltonian_pool.getPrimary()->informOperatorsOfListener();
  int num_walkers   = 4;
  int num_electrons = particle_pool.getParticleSet("e")->getTotalNum();
  int num_ions      = particle_pool.getParticleSet("ion")->getTotalNum();
  for (int iw = 0; iw < num_walkers; ++iw)
  {
    psets.emplace_back(pset_target);
    psets.back().randomizeFromSource(*particle_pool.getParticleSet("ion"));
    twfs.emplace_back(trial_wavefunction.makeClone(psets.back()));
    hams.emplace_back(hamiltonian_pool.getPrimary()->makeClone(psets.back(), *twfs.back()));
  }

  using MCPWalker = typename OperatorEstBase::MCPWalker;
  std::vector<OperatorEstBase::MCPWalker> walkers(num_walkers, MCPWalker(pset_target.getTotalNum()));

  auto walker_refs = makeRefVector<MCPWalker>(walkers);
  RefVectorWithLeader<MCPWalker> walker_list{walker_refs[0], walker_refs};

  RefVector<QMCHamiltonian> ham_refs = convertUPtrToRefVector(hams);

  RefVectorWithLeader<QMCHamiltonian> ham_list{ham_refs[0], ham_refs};
  ResourceCollection ham_res("test_ham_res");
  ham_list.getLeader().createResource(ham_res);
  ResourceCollectionTeamLock<QMCHamiltonian> ham_lock(ham_res, ham_list);
  e_den_est.registerListeners(ham_list.getLeader());

  auto p_refs = makeRefVector<ParticleSet>(psets);
  RefVectorWithLeader<ParticleSet> p_list{p_refs[0], p_refs};

  ResourceCollection pset_res("test_pset_res");
  p_list.getLeader().createResource(pset_res);
  ResourceCollectionTeamLock<ParticleSet> pset_lock(pset_res, p_list);

  auto twf_refs = convertUPtrToRefVector(twfs);
  RefVectorWithLeader<TrialWaveFunction> twf_list{twf_refs[0], twf_refs};

  ResourceCollection wfc_res("test_wfc_res");
  twf_list.getLeader().createResource(wfc_res);
  ResourceCollectionTeamLock<TrialWaveFunction> mw_wfc_lock(wfc_res, twf_list);


  p_refs[0].get().L[0] = 1.0;
  p_refs[1].get().L[1] = 1.0;
  p_refs[2].get().L[2] = 1.0;
  p_refs[3].get().L[3] = 1.0;

  if constexpr (generate_test_data)
  {
    std::cout << "QMCHamiltonian-registerListeners initialize psets with:\n{";

    for (int iw = 0; iw < num_walkers; ++iw)
    {
      //psets.emplace_back(pset_target);
      std::cout << "{";
      for (auto r : p_refs[iw].get().R)
        std::cout << NativePrint(r) << ",";
      std::cout << "},\n";
    }
    std::cout << "}\n";
  }
  else
  {
    std::vector<ParticleSet::ParticlePos> deterministic_rs = {
        {
            {0.515677886, 0.9486072745, -1.17423246},
            {-0.3166678423, 1.376550506, 1.445290031},
            {1.96071365, 2.47265689, 1.051449486},
            {0.745853269, 0.5551359072, 4.080774681},
            {-0.3515016103, -0.5192222523, 0.9941510909},
            {-0.8354426872, 0.7071638258, -0.3409843552},
            {0.4386044751, 1.237378731, 2.331874152},
            {2.125850717, 0.3221067321, 0.5825731561},
        },
        {
            {-0.4633736785, 0.06318772224, -0.8580153742},
            {-1.174926354, -0.6276503679, 0.07458759314},
            {1.327618206, 2.085829379, 1.415749862},
            {0.9114727103, 0.1789183931, -0.08135540251},
            {-2.267908723, 0.802928773, 0.9522812957},
            {1.502715257, -1.84493529, 0.2805620469},
            {3.168934617, 0.1348337978, 1.371092768},
            {0.8310229518, 1.070827168, 1.18016733},
        },
        {
            {-0.04847732172, -1.201739871, -1.700527771},
            {0.1589259538, -0.3096047065, -2.066626415},
            {2.255976232, 1.629132391, -0.8024446773},
            {2.534792993, 3.121092901, 1.609780703},
            {-0.2892376071, -0.152022511, -2.727613712},
            {0.2477154804, 0.5039232765, 2.995702733},
            {3.679345099, 3.037770313, 2.808899306},
            {0.6418578532, 1.935944544, 1.515637954},
        },
        {
            {0.91126951, 0.0234699242, 1.442297821},
            {-0.9240061217, -0.1014997844, 0.9081020061},
            {1.887794866, 2.210192703, 2.209118551},
            {2.758945014, -1.21140421, 1.3337907},
            {0.376540703, 0.3485486555, 0.9056881595},
            {-0.3512770187, -0.4056820917, -2.068499576},
            {0.5358460986, 2.720153363, 1.41999706},
            {2.284020089, 1.173071915, 1.044597715},
        },
    };
    for (int iw = 0; iw < num_walkers; ++iw)
    {
      p_refs[iw].get().R = deterministic_rs[iw];
    }
  }

  ParticleSet::mw_update(p_list);
  TrialWaveFunction::mw_evaluateLog(twf_list, p_list);
  QMCHamiltonian::mw_evaluate(ham_list, twf_list, p_list);

  StdRandom<double> rng;
  rng.init(101);

  e_den_est.accumulate(walker_refs, p_list, twf_list, rng);
}


} // namespace qmcplusplus
