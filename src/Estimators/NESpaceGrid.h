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
// Some code refactored from: QMCHamiltonian/SpaceGrid.h
//////////////////////////////////////////////////////////////////////////////////////

/** \file
 *  This is the port of QMCHamiltonian/SpaceGrid to the new Estimator design.
 *  Name clashes were undesirable as the legacy implementation needed to remain.
 *  \todo rename to more obvious Spacegrid once QMCHamiltonian/SpaceGrid is removed
 */
#ifndef QMCPLUSPLUS_NESPACEGRID_H
#define QMCPLUSPLUS_NESPACEGRID_H

#include <Configuration.h>

#include "SpaceGridInput.h"
#include "OhmmsPETE/Tensor.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "Pools/PooledData.h"
#include "QMCHamiltonians/ObservableHelper.h"
#include "Particle/DistanceTable.h"
#include "NEReferencePoints.h"

namespace qmcplusplus
{
/** SpaceGrid refactored for use with batched estimator design
 *  NE should be dropped when QMCHamiltonian/SpaceGrid has been deleted.
 *
 *  This class has more going on than just representing a spacial grid
 *  I'm still working out how much of this just because of the Voronoi code that shouldn't be
 *  part of the same class as the simpler grid and how much is particleset contimination etc.
 */
class NESpaceGrid
{
public:
  using Real        = QMCTraits::RealType;
  using Point       = typename NEReferencePoints::Point;
  using Points      = typename NEReferencePoints::Points;
  using BufferType  = PooledData<Real>;
  using POLT        = PtclOnLatticeTraits;
  using ParticlePos = POLT::ParticlePos;
  using AxTensor    = Tensor<Real, OHMMS_DIM>;
  enum class ReferenceEnergy
  {
    vacuum,
    neutral,
    noref
  };

  /** This constructor is used for electron only grids
   * \param[in]  sgi            input object for space grid.
   * \param[in]  reference      reference points from which origin and on other reference points referenced in input object are to be found
   * \param[in]  nvalues        number of fields the owning class wants for each grid point.
   * \param[in]  is_period      properly names is what is says
   */
  NESpaceGrid(SpaceGridInput& sgi, const Points& points, const int nvalues, const bool is_periodic);

  /** This is the general constructor
   * \param[in]  sgi            input object for space grid.
   * \param[in]  reference      reference points from which origin and on other reference points referenced in input object are to be found
   * \param[in]  ndp            number of ions that can move
   * \param[in]  nvalues        number of fields the owning class wants for each grid point.
   * \param[in]  is_period      properly names is what is says
   */
  NESpaceGrid(SpaceGridInput& sgi, const Points& points, const int ndp, const int nvalues, const bool is_periodic);

  void write_description(std::ostream& os, std::string& indent);
  int allocate_buffer_space(BufferType& buf);
  void registerCollectables(std::vector<ObservableHelper>& h5desc, hid_t gid, int grid_index) const;
  void evaluate(const ParticlePos& R,
                const Matrix<Real>& values,
                BufferType& buf,
                std::vector<bool>& particles_outside,
                const DistanceTableAB& dtab);

  bool check_grid(void);
  int nDomains(void) const { return ndomains_; }

  void sum(const BufferType& buf, Real* vals);

private:
  // The following funciton are static to provide some discipline and actual
  // visibility of there data dependence and effect on the object state when
  // called.
  /** Initialize NESpaceGrid for rectilinear grid
   *  \param[in]  input     SpaceGridInput object
   *  \param[in]  points    ReferencePoints object for grid
   *  Causes side effects updating
   *    origin_    fixed up origin for grid
   *    axes_      axes with scaling applied to it.
   *    axinv_     the inverse of the axes with scaling applied   
   */
  bool initializeRectilinear(const SpaceGridInput& input, const Points& points);
  /** Deal with some Axes issues
   *  \param[in]  input   input object
   *  \param[out] axes    axes with scaling applied to it.
   *  \param[out] axinv   the inverse of the axes with scaling applied
   *  could just be used for reporting.
   */

  /** Another function to cut scopes to sort of manageable size.
   *  does nothing but create many side effects
   */
  void someMoreAxisGridStuff();
    
  static void processAxis(const SpaceGridInput& input_, AxTensor& axes, AxTensor& axinv);

  static bool checkAxisGridValues(const SpaceGridInput& input_, const AxTensor& axes);
  
  SpaceGridInput& input_;
  int ndparticles_;
  bool is_periodic_;
  const Points& points_;

  // _NO_
  int buffer_start_;
  int buffer_end_;

  //private:

  int buffer_offset_; /// Assuming this is just written into a shared buffer of type Real
  int ndomains_;
  int nvalues_per_domain_;
  Matrix<Real> domain_volumes_;
  Matrix<Real> domain_centers_;

  std::vector<int> reference_count_;

  //really only used for cartesian-like grids
  Point origin_;
  AxTensor axes_;
  AxTensor axinv_;
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
