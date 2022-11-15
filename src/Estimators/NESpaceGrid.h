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
// Some code refactored from: QMCHamiltonian/SpaceGrid.h
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_NESPACEGRID_H
#define QMCPLUSPLUS_NESPACEGRID_H

#include <Configuration.h>
#include "OhmmsPETE/Tensor.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "Pools/PooledData.h"
#include "QMCHamiltonians/ObservableHelper.h"
#include "Particle/DistanceTable.h"

namespace qmcplusplus
{
/** SpaceGrid refactored for use with batched estimator design
 *  NE should be dropped when QMCHamiltonian/SpaceGrid has been deleted.
 */
class NESpaceGrid
{
public:
  using Real   = QMCTraits::RealType;
  using Points = typename NEReferencePoints::Points;
  using BufferType = PooledData<Real>;
  using Matrix_t   = Matrix<RealType>;
  using POLT    = PtclOnLatticeTraits;
  using ParticlePos = POLT::ParticlePos;

  enum class ReferenceEnergy
  {
    vacuum,
    neutral,
    noref
  };

  
  NESpaceGrid(int& nvalues);
  NESpaceGrid(SpaceGridInput& sgi,
           NEReferencePoints::Points& points,
           int ndp,
	      bool is_periodic);

  bool initialize_voronoi(std::map<std::string, Point>& points);
  void write_description(std::ostream& os, std::string& indent);
  int allocate_buffer_space(BufferType& buf);
  void registerCollectables(std::vector<ObservableHelper>& h5desc, hid_t gid, int grid_index) const;
  void evaluate(const ParticlePos& R,
                const Matrix<Real>& values,
                BufferType& buf,
                std::vector<bool>& particles_outside,
                const DistanceTableAB& dtab);

  bool check_grid(void);
  inline int nDomains(void) { return ndomains; }

  void sum(const BufferType& buf, Real* vals);

private:
  bool initialize_rectilinear(SpaceGridInput& input, Points& points);

  SpaceGridInput& input_;
  int ndparticles_;
  bool periodic_;
  NEReferencePoints& points_;
  
  // _NO_
  int buffer_start_;
  int buffer_end_;

  //private:

  int buffer_offset_;
  int ndomains_;
  int nvalues_per_domain_;
  Matrix<Real> domain_volumes_;
  Matrix<Real> domain_centers_;

  //in use if sorting by particle count
  bool chempot_;
  int npmin_, npmax_;
  int npvalues_;
  Matrix<Real> cellsamples_;

  std::vector<int> reference_count_;

  //really only used for cartesian-like grids
  Point origin_;
  Tensor<Real, OHMMS_DIM> axes_;
  Tensor<Real, OHMMS_DIM> axinv_;
  Real volume_;
  Matrix<Real> domain_uwidths_;
  std::string axlabel_[OHMMS_DIM];
  std::vector<int> gmap_[OHMMS_DIM];
  Real odu_[OHMMS_DIM];
  Real umin_[OHMMS_DIM];
  Real umax_[OHMMS_DIM];
  int dimensions_[OHMMS_DIM];
  int dm_[OHMMS_DIM];
  ReferenceEnergy reference_energy_;
  
  struct IRPair
  {
    Real r;
    int i;
  };
  std::vector<IRPair> nearcell_;

  //used only in evaluate
  // Then why are they here and not on the stack in evaluate!
  Point u_, ub_;
};


} // namespace qmcplusplus

#endif
