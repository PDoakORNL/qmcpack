//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2023 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Peter W. Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// Some code refactored from: QMCHamiltonian/EnergyDensityEstimator.cpp
//////////////////////////////////////////////////////////////////////////////////////

#include "EnergyDensityEstimator.h"

namespace qmcplusplus
{

struct PosCharge
{
  EnergyDensityEstimator::ParticlePos r_ptcls;
  std::vector<EnergyDensityEstimator::Real> z_ptcls;
};

auto EnergyDensityEstimator::extractIonPositionsAndCharge(const ParticleSet& pset)
{
  const SpeciesSet& species(pset.getSpeciesSet());
  int charge_index = species.findAttribute("charge");
  int nspecies     = species.TotalNum;
  int nps          = pset.getTotalNum();
  std::vector<Real> z_spec;
  z_spec.resize(nspecies);
  std::vector<Real> z_ptcls;
  z_ptcls.resize(nps);
  for (int spec = 0; spec < nspecies; spec++)
    z_spec[spec] = species(charge_index, spec);
  for (int i = 0; i < nps; i++)
    z_ptcls[i] = z_spec[pset.GroupID[i]];
  ParticlePos r_ptcls;
  r_ptcls.resize(pset.R.size());
  for (int i = 0; i < pset.R.size(); i++)
    r_ptcls[i] = pset.R[i];
  if (pset.getLattice().SuperCellEnum != SUPERCELL_OPEN)
    pset.applyMinimumImage(r_ptcls);
  return PosCharge{r_ptcls, z_ptcls};
}

EnergyDensityEstimator::EnergyDensityEstimator(const EnergyDensityInput& input,
                                               const PSPool& pset_pool,
                                               DataLocality data_locality)
    : OperatorEstBase(data_locality), input_(input), pset_dynamic_(getParticleSet(pset_pool, input.get_dynamic()))
{
  // Bringing a bunch of adhoc setup from legacy.
  // removed redundant turnOnPerParticleSK this is handled via the EstimatorManagerNew if listeners are detected through
  // CoulombPBCAA{AB} which are the actual operators that need it.
  // Although it is possible that our copies of the particle sets will need per particle structure factors I don't think
  // they do.
  RefVector<ParticleSet> pset_refs;
  if (input.get_static().empty())
  {
    dtable_index_ = -1;
  }
  else
  {
    pset_static_.emplace(getParticleSet(pset_pool, input.get_static()));
    dtable_index_ = pset_dynamic_.addTable(pset_static_.value());
    pset_refs.push_back(pset_static_.value());
    if (!input.get_ion_points())
      n_particles_ += pset_static_->getTotalNum();
  }
  r_work_.resize(n_particles_);
  ed_values_.resize(n_particles_, N_EDVALS);
  if (input.get_ion_points())
  {
    n_ions_ = pset_static_->getTotalNum();
    ed_ion_values_.resize(n_ions_, N_EDVALS);
    r_ion_work_.resize(n_ions_, OHMMS_DIM);
    for (int i = 0; i < n_ions_; ++i)
      for (int d = 0; d < OHMMS_DIM; ++d)
        r_ion_work_(i, d) = pset_static_->R[i][d];
  }
  particles_outside_.resize(n_particles_, true);
  ref_points_ = std::make_unique<NEReferencePoints>(input.get_ref_points_input(), pset_dynamic_, pset_refs);

  bool periodic = pset_dynamic_.getLattice().SuperCellEnum != SUPERCELL_OPEN;

  auto space_grid_inputs = input.get_space_grid_inputs();
  spacegrids_.reserve(space_grid_inputs.size());
  for (int ig = 0; ig < spacegrids_.size(); ++ig)
  {
    if (pset_static_)
    {
      auto [r_ptcls, z_ptcls] = extractIonPositionsAndCharge(*pset_static_);
      spacegrids_[ig]         = std::make_unique<NESpaceGrid>(space_grid_inputs[ig], ref_points_->get_points(), r_ptcls,
                                                      z_ptcls, pset_dynamic_.getTotalNum(), N_EDVALS, periodic);
    }
    else
      spacegrids_[ig] =
          std::make_unique<NESpaceGrid>(space_grid_inputs[ig], ref_points_->get_points(), N_EDVALS, periodic);
  }
}

EnergyDensityEstimator::EnergyDensityEstimator(const EnergyDensityEstimator& ede, const DataLocality dl)
    : OperatorEstBase(dl), input_(ede.input_), pset_dynamic_(ede.pset_dynamic_)
{
  data_locality_ = dl;
  ref_points_    = std::make_unique<NEReferencePoints>(*ref_points_);
}

EnergyDensityEstimator::~EnergyDensityEstimator() {}

void EnergyDensityEstimator::registerListeners(QMCHamiltonian& ham_leader)
{
  ListenerVector<Real> kinetic_listener("kinetic", getListener(kinetic_values_));
  QMCHamiltonian::mw_registerKineticListener(ham_leader, kinetic_listener);
  ListenerVector<Real> potential_listener("potential", getListener(local_pot_values_));
  QMCHamiltonian::mw_registerLocalPotentialListener(ham_leader, potential_listener);
  ListenerVector<Real> ion_potential_listener("potential", getListener(local_ion_pot_values_));
  QMCHamiltonian::mw_registerLocalIonPotentialListener(ham_leader, ion_potential_listener);
}

/** This function collects the per particle energies.
 *  right now these are indentified by a string for each type.  This could be optimized but
 *  could also be an insiginificant cost versus the frequently large number of values handled.
 *  The values themselves are a vector of size particle_num.
 */
ListenerVector<QMCTraits::RealType>::ReportingFunction EnergyDensityEstimator::getListener(
    CrowdEnergyValues<Real>& values)
{
  auto& local_values = values;
  return [&local_values](const int walker_index, const std::string& name, const Vector<Real>& inputV) {
    if (walker_index >= local_values[name].size())
      local_values[name].resize(walker_index + 1);
    local_values[name][walker_index] = inputV;
  };
}

const ParticleSet& EnergyDensityEstimator::getParticleSet(const PSPool& psetpool, const std::string& psname) const
{
  auto pset_iter(psetpool.find(psname));
  if (pset_iter == psetpool.end())
  {
    throw UniformCommunicateError("Particle set pool does not contain \"" + psname +
                                  "so EnergyDensityEstimator::get_particleset fails!");
  }
  return *(pset_iter->second.get());
}

void EnergyDensityEstimator::accumulate(const RefVector<MCPWalker>& walkers,
                                        const RefVector<ParticleSet>& psets,
                                        const RefVector<TrialWaveFunction>& wfns,
                                        RandomGenerator& rng)
{
  combinePerParticleEnergies(local_pot_values_, reduced_local_pot_values_);
  combinePerParticleEnergies(kinetic_values_, reduced_local_kinetic_values_);
  if (pset_static_)
    combinePerParticleEnergies(local_ion_pot_values_, reduced_local_ion_pot_values_);
  for (int iw = 0; iw < walkers.size(); ++iw)
  {
    walkers_weight_ += walkers[iw].get().Weight;
    evaluate(psets[iw], walkers[iw], iw);
  }
}

void EnergyDensityEstimator::evaluate(ParticleSet& pset, const MCPWalker& walker, const int walker_index)
{
  //Collect positions from ParticleSets
  int p_count = 0;
  {
    const ParticlePos& Rs = pset.R;
    for (int i = 0; i < Rs.size(); i++)
    {
      r_work_[p_count] = Rs[i];
      p_count++;
    }
  }
  if (pset_static_ && !input_.get_ion_points())
  {
    const ParticlePos& Rs = pset_static_->R;
    for (int i = 0; i < Rs.size(); i++)
    {
      r_work_[p_count] = Rs[i];
      p_count++;
    }
  }
  if (pset.getLattice().SuperCellEnum != SUPERCELL_OPEN)
    pset.applyMinimumImage(r_work_);
  //Convert information accumulated in ParticleSets into EnergyDensity quantities
  Real weight = walker.Weight;
  p_count     = 0;
  {
    const auto& Ts = reduced_local_kinetic_values_[walker_index];
    const auto& Vd = reduced_local_pot_values_[walker_index];
    for (int i = 0; i < pset.getTotalNum(); i++)
    {
      ed_values_(p_count, W) = weight;
      ed_values_(p_count, T) = weight * Ts[i];
      ed_values_(p_count, V) = weight * Vd[i];
      p_count++;
    }
  }
  if (pset_static_)
  {
    const ParticleSet& Ps = *pset_static_;
    const auto& Vs        = reduced_local_ion_pot_values_[walker_index];
    if (!input_.get_ion_points())
      for (int i = 0; i < Ps.getTotalNum(); i++)
      {
        ed_values_(p_count, W) = weight;
        ed_values_(p_count, T) = 0.0;
        ed_values_(p_count, V) = weight * Vs[i];
        p_count++;
      }
    else
      for (int i = 0; i < Ps.getTotalNum(); i++)
      {
        ed_ion_values_(i, W) = weight;
        ed_ion_values_(i, T) = 0.0;
        ed_ion_values_(i, V) = weight * Vs[i];
      }
  }
  //Accumulate energy density in spacegrids
  const auto& dtab(pset.getDistTableAB(dtable_index_));
  fill(particles_outside_.begin(), particles_outside_.end(), true);
  for (int i = 0; i < spacegrids_.size(); i++)
  {
    NESpaceGrid& sg = *(spacegrids_[i]);
    sg.accumulate(r_work_, ed_values_, particles_outside_, dtab);
  }

  //Accumulate energy density of particles outside any spacegrid
  int bi, v;
  const int bimax = outside_buffer_offset + N_EDVALS;
  for (int p = 0; p < particles_outside_.size(); p++)
  {
    if (particles_outside_[p])
    {
      for (bi = outside_buffer_offset, v = 0; bi < bimax; bi++, v++)
      {
        data_[bi] += ed_values_(p, v);
      }
    }
  }
  if (input_.get_ion_points())
  {
    // Accumulate energy density for ions at a point field
    bi = outside_buffer_offset + N_EDVALS;
    for (int i = 0; i < n_ions_; i++)
      for (v = 0; v < N_EDVALS; v++, bi++)
      {
        data_[bi] += ed_ion_values_(i, v);
      }
  }
  nsamples++;
}

void EnergyDensityEstimator::collect(const RefVector<OperatorEstBase>& type_erased_operator_estimators)
{
  int num_crowds = type_erased_operator_estimators.size();
  for (int ig = 0; ig < spacegrids_.size(); ++ig)
  {
    RefVector<NESpaceGrid> crowd_grids;
    crowd_grids.reserve(num_crowds);
    for (OperatorEstBase& crowd_oeb : type_erased_operator_estimators)
    {
      EnergyDensityEstimator& crowd_ede = dynamic_cast<EnergyDensityEstimator&>(crowd_oeb);
      NESpaceGrid& grid_ref             = *(crowd_ede.spacegrids_[ig]);
      crowd_grids.push_back(grid_ref);
    }
    NESpaceGrid::collect(*(spacegrids_[ig]), crowd_grids);
  }
  OperatorEstBase::collect(type_erased_operator_estimators);
}

void EnergyDensityEstimator::registerOperatorEstimator(hdf_archive& file)
{
  hdf_path hdf_name{my_name_};
  h5desc_.emplace_back(hdf_name / "variables");
  auto& oh = h5desc_.back();
  oh.addProperty(const_cast<int&>(n_particles_), "nparticles", file);
  int nspacegrids = spacegrids_.size();
  oh.addProperty(const_cast<int&>(nspacegrids), "nspacegrids", file);
  oh.addProperty(const_cast<int&>(nsamples), "nsamples", file);

  if (input_.get_ion_points())
  {
    oh.addProperty(const_cast<int&>(n_ions_), "nions", file);
    oh.addProperty(const_cast<Matrix<Real>&>(r_ion_work_), "ion_positions", file);
  }

  ref_points_->save(h5desc_, file);

  h5desc_.emplace_back(hdf_name / "outside");
  auto& ohOutside = h5desc_.back();
  std::vector<int> ng(1);
  ng[0] = N_EDVALS;
  ohOutside.set_dimensions(ng, outside_buffer_offset);
  for (int i = 0; i < spacegrids_.size(); i++)
  {
    spacegrids_[i]->registerGrid(file, h5desc_, i);
  }
  if (input_.get_ion_points())
  {
    std::vector<int> ng2(2);
    ng2[0] = n_ions_;
    ng2[1] = N_EDVALS;

    h5desc_.emplace_back(hdf_name / "ions");
    auto& ohIons = h5desc_.back();
    ohIons.set_dimensions(ng2, 0);
  }
}

std::unique_ptr<OperatorEstBase> EnergyDensityEstimator::spawnCrowdClone() const
{
  auto spawn_data_locality = data_locality_;
  auto data_size           = this->data_.size();
  UPtr<EnergyDensityEstimator> spawn(std::make_unique<EnergyDensityEstimator>(*this, spawn_data_locality));
  spawn->get_data().resize(data_size);
  return spawn;
}

void EnergyDensityEstimator::startBlock(int steps) {}

} // namespace qmcplusplus
