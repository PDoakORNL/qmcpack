//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File refactored from: EnergyDensityEstimator.h
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_ENERGY_DENSITY_ESTIMATOR_NEW_H
#define QMCPLUSPLUS_ENERGY_DENSITY_ESTIMATOR_NEW_H
#include "OperatorEstBase.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "QMCWaveFunctions/WaveFunctionTypes.hpp"
#include "QMCHamiltonians/ReferencePoints.h"
#include "QMCHamiltonians/SpaceGrid.h"
#include <map>
#include <vector>

namespace qmcplusplus
{
  
class EnergyDensityEstimatorNew : public OperatorEstBase, public PtclOnLatticeTraits
{
public:
  using WFT = WaveFunctionTypes<QMCTraits::RealType, QMCTraits::FullPrecRealType>;
  using Real = WFT::Real;
  
  using Point  = ReferencePoints::Point;
  using PSPool = std::map<std::string, const std::unique_ptr<ParticleSet>>;

  EnergyDensityEstimator(EnergyDensityInput&& edi,
			 const ParticleSet& pset_target,
			 const QMCHamiltonian& ham);
  
  EnergyDensityEstimator(const PSPool& PSP, const std::string& defaultKE);
  ~EnergyDensityEstimator() override;

  void accumulate(const RefVector<MCPWalker>& walkers,
                  const RefVector<ParticleSet>& psets,
                  const RefVector<TrialWaveFunction>& wfns,
                  RandomGenerator& rng) override;
  
  void addObservables(PropertySetType& plist) {}
  void addObservables(PropertySetType& plist, BufferType& olist) override;
  void registerCollectables(std::vector<ObservableHelper>& h5desc, hid_t gid) const override;
  void setObservables(PropertySetType& plist) override;
  void setParticlePropertyList(PropertySetType& plist, int offset) override;
  bool put(xmlNodePtr cur) override;
  bool put(xmlNodePtr cur, ParticleSet& P);
  bool get(std::ostream& os) const override;
  std::unique_ptr<OperatorBase> makeClone(ParticleSet& qp, TrialWaveFunction& psi) final;

  void write_description(std::ostream& os);

  void write_EDValues(void);
  void write_nonzero_domains(const ParticleSet& P);
  void write_Collectables(std::string& label, int& cnt, ParticleSet& P);
  void registerListeners(QMCHamiltonian& hamiltonian);
  
private:
  //system information
  std::string defKE;
  const PSPool& psetpool;
  ParticleSet* Pdynamic;
  ParticleSet* Pstatic;
  ParticleSet* get_particleset(std::string& psname);
  int dtable_index;
  int nparticles;
  bool ion_points;
  int nions;
  int ion_buffer_offset;
  Matrix<RealType> Rion;
  //collection of points from which to build spacegrid origin and axes
  ReferencePoints ref;
  //EnergyDenstity quantities
  enum
  {
    W = 0,
    T,
    V,
    nEDValues
  };
  Matrix<RealType> EDValues;
  Matrix<RealType> EDIonValues;
  //for EnergyDensity of particles falling outside any spacegrid
  int outside_buffer_offset;
  std::vector<bool> particles_outside;
  //spacegrids are used to find which cell domain
  //  contains the Energy information of particles
  std::vector<SpaceGrid*> spacegrids;
  //particle positions
  ParticlePos R;
  //number of samples accumulated
  int nsamples;

  //needed (temporarily) for chempot
  //ParticleSet should carry Zptcl so it doesn't have
  // to be computed everywhere from species
  std::vector<RealType> Zptcl;
  ParticlePos Rptcl;
  void set_ptcl(void);
  void unset_ptcl(void);

  /// Should be a per walker quantity
  std::vector<std::vector<Real>> weight_samples_;

  /// per walker per particle kinetic energy
  std::vector<std::vector<Real>> kinetic_samples_;

  std::vector<Vector<Read>> local_potential_samples_;
  
  std::vector<Real> position_samples_;
  std::vector<Vector<Real>> vd_samples_;
  std::vector<Vector<Real>> vs_samples_;  
};


} // namespace qmcplusplus
#endif
