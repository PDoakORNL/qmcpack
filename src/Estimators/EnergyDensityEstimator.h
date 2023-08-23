//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2023 QMCPACK developers.
//
// File developed by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Peter W. Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// Some code refactored from: QMCHamiltonian/EnergyDensityEstimator.h
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_ESTIMATOR_ENERGY_DENSITY_H
#define QMCPLUSPLUS_ESTIMATOR_ENERGY_DENSITY_H

#include "OperatorEstBase.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "NEReferencePoints.h"
#include "NESpaceGrid.h"
#include "EnergyDensityInput.h"

namespace qmcplusplus
{
class NEEnergyDensityEstimator : public OperatorEstBase
{
public:
  using Real        = QMCTraits::RealType;
  using FullPrecReal = QMCTraits::FullPrecRealType;
  using Point       = typename NEReferencePoints::Point;
  using Points      = typename NEReferencePoints::Points;
  using POLT        = PtclOnLatticeTraits;
  using ParticlePos = POLT::ParticlePos;
  using PSPool      = std::map<std::string, const UPtr<ParticleSet>>;
  
  NEEnergyDensityEstimator(const EnergyDensityInput& input, const PSPool& PSP, DataLocality dl = DataLocality::crowd);

  NEEnergyDensityEstimator(const NEEnergyDensityEstimator& ede, const DataLocality dl);

private:
  /** shared construction code
   *  stateful setup from legacy.
   * removed redundant turnOnPerParticleSK this is handled via the EstimatorManagerNew if listeners are detected through
   * CoulombPBCAA{AB} which are the actual operators that need it.
   * Although it is possible that our copies of the particle sets will need per particle structure factors I don't think
   * they do. I think it is needed for the actual particle sets involved in walking.
   */
  void constructToReferencePoints(ParticleSet& pset_dynamic, const std::optional<ParticleSet>& pset_static);

public:
  /** @ingroup OperatorEstBase Overrides
   *  should NESpaceGrid have anything to do with OperatorEstBase API?
   *  @{
   */
  ~NEEnergyDensityEstimator() override;

  /** Register listeners for energy values
   */
  void registerListeners(QMCHamiltonian& ham_leader) override;

  ListenerVector<Real>::ReportingFunction getListener();

  /** accumulate 1 or more walkers of EnergyDensity samples
   */
  void accumulate(const RefVector<MCPWalker>& walkers,
                  const RefVector<ParticleSet>& psets,
                  const RefVector<TrialWaveFunction>& wfns,
                  RandomBase<FullPrecReal>& rng) override;

  void evaluate(ParticleSet& pset, const MCPWalker& walker, const int walker_index);

  /** this allows the EstimatorManagerNew to reduce without needing to know the details
   *  of SpaceGrid's data.
   *
   *  can use base class default until crowd level MomentumDistribution
   *  estimators don't have a copy of the density grid.
   */
  void collect(const RefVector<OperatorEstBase>& type_erased_operator_estimators) override;
  UPtr<OperatorEstBase> spawnCrowdClone() const override;

  /** start block entry point
   */
  void startBlock(int steps) override;

  /// @}

  /** return lambda function to register as listener
   *  the purpose of this function is to factor out the production of the lambda for unit testing
   *  \param[out] values
   */
  ListenerVector<Real>::ReportingFunction getListener(CrowdEnergyValues<Real>& values);

  const ParticleSet& getParticleSet(const PSPool& psetpool, const std::string& psname) const;

  void registerOperatorEstimator(hdf_archive& file) override;

  size_t getFullDataSize();

  RefVector<NESpaceGrid> getSpaceGrids();
private:

  auto extractIonPositionsAndCharge(const ParticleSet& pset);

  EnergyDensityInput input_;

  CrowdEnergyValues<Real> kinetic_values_;
  CrowdEnergyValues<Real> local_pot_values_;
  CrowdEnergyValues<Real> local_ion_pot_values_;

  /// working space for reduced_values_;
  std::vector<Vector<Real>> reduced_local_pot_values_;
  std::vector<Vector<Real>> reduced_local_kinetic_values_;
  std::vector<Vector<Real>> reduced_local_ion_pot_values_;

  ParticleSet pset_dynamic_;
  std::optional<ParticleSet> pset_static_;
  int dtable_index_ = -1;

  int n_particles_;
  int n_ions_;

  UPtr<NEReferencePoints> ref_points_;

  /// unboxed SpaceGridInputs for child and clone objects to refrence
  std::vector<SpaceGridInput> spacegrid_inputs_;
  
  /** EnergyDenstity quantities
   *  legacy style enum into vector not ideal.
   */
  enum
  {
    W = 0,
    T,
    V,
    N_EDVALS
  };

  /** @ingroup outside particles
   *  for the EnergyDensity of particles falling outside any spacegrid
   *  @{ */
  /** in legacy the starting index into the collectibles buffer.
   *  Maintained to use more legacy code without modificaiton in the short term.
   *  In the long term its possible the entire way the grid data is structured in memory should be redesigned.
   */
  int outside_buffer_offset{0};
  std::vector<bool> particles_outside_;

  /// @}

  /// spacegrids are used to find which cell domain contains the Energy information of particles
  UPtrVector<NESpaceGrid> spacegrids_;

  /** @ingroup working variables reused for each walker
   *  Do not use for persistent state.
   *  @{ */

  //particle positions
  ParticlePos r_work_;
  Matrix<Real> r_ion_work_;

  // values, seems like these should just be return values
  Matrix<Real> ed_values_;
  Matrix<Real> ed_ion_values_;
  /// @}


  //number of samples accumulated
  int nsamples;

  DataLocality data_locality_{DataLocality::crowd};

  // //needed (temporarily) for chempot
  // //ParticleSet should carry Zptcl so it doesn't have
  // // to be computed everywhere from species
  // std::vector<Real> Zptcl;
  // ParticlePos Rptcl;
};
} // namespace qmcplusplus

#endif
