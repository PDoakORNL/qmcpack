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
// Some code refactored from: QMCHamiltonian/SpaceGrid.cpp
//////////////////////////////////////////////////////////////////////////////////////

#include "NESpaceGrid.h"
#include "OhmmsData/AttributeSet.h"
#include "Utilities/string_utils.h"
#include <cmath>
#include "OhmmsPETE/OhmmsArray.h"
#include "NEReferencePoints.h"
#include "Concurrency/OpenMP.h"

namespace qmcplusplus
{
using std::acos;
using std::atan2;
using std::cos;
using std::floor;
using std::sin;
using std::sqrt;

NESpaceGrid::NESpaceGrid(SpaceGridInput& sgi,
                         const NEReferencePoints::Points& points,
                         const int nvalues,
                         const bool is_periodic)
    : NESpaceGrid(sgi, points, 0, nvalues, is_periodic)
{}

NESpaceGrid::NESpaceGrid(SpaceGridInput& sgi,
                         const NEReferencePoints::Points& points,
                         const int ndp,
                         const int nvalues,
                         const bool is_periodic)
    : input_(sgi), ndparticles_(ndp), is_periodic_(is_periodic), points_(points), nvalues_per_domain_(nvalues)
{
  ndomains_         = -1;
  bool init_success = initializeRectilinear(input_, points_);
  if (!init_success)
    throw std::runtime_error("NESpaceGrid initialization failed");
}

void NESpaceGrid::processAxis(const SpaceGridInput& input, AxTensor& axes, AxTensor& axinv)
{
  auto& axis_labels = input.get_axis_labels();
  auto& axis_p1s    = input.get_axis_p1s();
  auto& axis_p2s    = input.get_axis_p2s();
  auto& axis_scales = input.get_axis_scales();
  auto& axis_grids  = input.get_axis_grids();

  // SpaceGridInput check particular validity insures the number of axis == OHMMS_DIM
  for (int iaxis = 0; iaxis < OHMMS_DIM; ++iaxis)
  {
    Real frac = axis_scales[iaxis];
    for (int d = 0; d < OHMMS_DIM; d++)
      axes(d, iaxis) = frac * (axis_p1s[iaxis][d] - axis_p2s[iaxis][d]);
  }
  axinv = inverse(axes);
}

bool NESpaceGrid::initializeRectilinear(const SpaceGridInput& input, const Points& points_)
{
  // This code should be refactored to SpaceGridInput such that a simple map of
  // axis is available.
  using CoordForm = typename SpaceGridInput::CoordForm;
  std::map<std::string, int> cmap;
  auto coord_form = input.get_coord_form();

  // Belongs in input class?
  bool is_periodic = input.isPeriodic();

  const std::string& origin_p1 = input.get_origin_p1();
  if (origin_p1.size() > 0)
  {
    const std::string& origin_p2 = input.get_origin_p2();
    origin_ = points_.at(origin_p1) + input.get_origin_fraction() * (points_.at(origin_p2) - points_.at(origin_p1));
  }
  else
    origin_ = points_.at("zero");

  // variables for loop
  Real utol = 1e-5;
  std::string grid;
  std::vector<int> ndu_per_interval[OHMMS_DIM];

  // awful state variables for processAxes
  int iaxis = 0;
  int naxes = 0;
  processAxis(input, axes_, axinv_);
  bool succeeded = checkAxisGridValues(input, axes_);

  someMoreAxisGridStuff();
  return succeeded;
}

void NESpaceGrid::someMoreAxisGridStuff()
{
  auto& axis_grids = input_.get_axis_grids();
  // This dates back to the legacy implementation and I'm not sure why both code blocks are here.
  //set grid dimensions
  // C/Python style indexing
  dm_[0] = axis_grids[1].dimensions * axis_grids[2].dimensions;
  dm_[1] = axis_grids[2].dimensions;
  dm_[2] = 1;
  // Fortran/Matlab style indexing
  //dm[0] = 1;
  //dm[1] = dimensions[0];
  //dm[2] = dimensions[0]*dimensions[1];

  ndomains_ = axis_grids[0].dimensions * axis_grids[1].dimensions * axis_grids[2].dimensions;
  volume_   = std::abs(det(axes_)) * 8.0; //axes span only one octant
  //compute domain volumes, centers, and widths

  domain_volumes_.resize(ndomains_, 1);
  domain_centers_.resize(ndomains_, OHMMS_DIM);
  domain_uwidths_.resize(ndomains_, OHMMS_DIM);
  std::vector<Real> interval_centers[OHMMS_DIM];
  std::vector<Real> interval_widths[OHMMS_DIM];

  auto& agr = axis_grids;

  for (int d = 0; d < OHMMS_DIM; d++)
  {
    int nintervals = agr[d].ndu_per_interval.size();
    app_log() << "nintervals " << d << " " << nintervals << std::endl;
    interval_centers[d].resize(nintervals);
    interval_widths[d].resize(nintervals);
    interval_widths[d][0]  = agr[d].ndu_per_interval[0] / agr[d].odu;
    interval_centers[d][0] = interval_widths[d][0] / 2.0 + agr[d].umin;
    for (int i = 1; i < nintervals; i++)
    {
      interval_widths[d][i]  = agr[d].ndu_per_interval[i] / agr[d].odu;
      interval_centers[d][i] = interval_centers[d][i - 1] + .5 * (interval_widths[d][i] + interval_widths[d][i - 1]);
    }
    //app_log()<<"  interval widths"<< std::endl;
    //for(int i=0;i<nintervals;i++)
    //  app_log()<<"    "<<i<<" "<<interval_widths[d][i]<< std::endl;
    //app_log()<<"  interval centers"<< std::endl;
    //for(int i=0;i<nintervals;i++)
    //  app_log()<<"    "<<i<<" "<<interval_centers[d][i]<< std::endl;
  }
  Point du, uc, ubc, rc;
  Real vol        = 0.0;
  Real vscale     = std::abs(det(axes_));
  using CoordForm = SpaceGridInput::CoordForm;
  for (int i = 0; i < agr[0].dimensions; i++)
  {
    for (int j = 0; j < agr[1].dimensions; j++)
    {
      for (int k = 0; k < agr[2].dimensions; k++)
      {
        int idomain = dm_[0] * i + dm_[1] * j + dm_[2] * k;
        du[0]       = interval_widths[0][i];
        du[1]       = interval_widths[1][j];
        du[2]       = interval_widths[2][k];
        uc[0]       = interval_centers[0][i];
        uc[1]       = interval_centers[1][j];
        uc[2]       = interval_centers[2][k];
        switch (input_.get_coord_form())
        {
        case (CoordForm::CARTESIAN):
          vol = du[0] * du[1] * du[2];
          ubc = uc;
          break;
        case (CoordForm::CYLINDRICAL):
          uc[1]  = 2.0 * M_PI * uc[1] - M_PI;
          du[1]  = 2.0 * M_PI * du[1];
          vol    = uc[0] * du[0] * du[1] * du[2];
          ubc[0] = uc[0] * cos(uc[1]);
          ubc[1] = uc[0] * sin(uc[1]);
          ubc[2] = uc[2];
          break;
        case (CoordForm::SPHERICAL):
          uc[1] = 2.0 * M_PI * uc[1] - M_PI;
          du[1] = 2.0 * M_PI * du[1];
          uc[2] = M_PI * uc[2];
          du[2] = M_PI * du[2];
          vol   = (uc[0] * uc[0] + du[0] * du[0] / 12.0) * du[0] //r
              * du[1]                                            //theta
              * 2.0 * sin(uc[2]) * sin(.5 * du[2]);              //phi
          ubc[0] = uc[0] * sin(uc[2]) * cos(uc[1]);
          ubc[1] = uc[0] * sin(uc[2]) * sin(uc[1]);
          ubc[2] = uc[0] * cos(uc[2]);
          break;
        default:
          break;
        }
        vol *= vscale;
        rc = dot(axes_, ubc) + origin_;
        //app_log()<< std::endl;
        //app_log()<<"umin "<<uc-du/2<< std::endl;
        //app_log()<<"uc "<<uc<< std::endl;
        //app_log()<<"umax "<<uc+du/2<< std::endl;
        //app_log()<<"rc "<<rc<< std::endl;
        domain_volumes_(idomain, 0) = vol;
        for (int d = 0; d < OHMMS_DIM; d++)
        {
          domain_uwidths_(idomain, d) = du[d];
          domain_centers_(idomain, d) = rc[d];
        }
      }
    }
  }
  ////the following check is only valid if grid spans maximum amount
  ////check that the amount of space the grid takes up is correct
  //Real vfrac;
  //switch(coordinate){
  //case(cartesian):
  //  vfrac=1.0;
  //  break;
  //case(cylindrical):
  //  vfrac=M_PI/4.0;
  //  break;
  //case(spherical):
  //  vfrac=M_PI/6.0;
  //  break;
  //default:
  //  vfrac=vol_tot/volume;
  //}
  //if(std::abs(vol_tot/volume-vfrac)>1e-6){
  //  app_log()<<"  "<<coord<<" relative volume"<< std::endl;
  //  app_log()<<"  spacegrid volume fraction "<<vol_tot/volume<< std::endl;
  //  app_log()<<"                  should be "<<vfrac<< std::endl;
  //  succeeded=false;
  //}
  //find the actual volume of the grid
  for (int d = 0; d < OHMMS_DIM; d++)
  {
    du[d] = agr[d].umax - agr[d].umin;
    uc[d] = .5 * (agr[d].umax + agr[d].umin);
  }
  switch (input_.get_coord_form())
  {
  case CoordForm::CARTESIAN:
    vol = du[0] * du[1] * du[2];
    break;
  case CoordForm::CYLINDRICAL:
    uc[1] = 2.0 * M_PI * uc[1] - M_PI;
    du[1] = 2.0 * M_PI * du[1];
    vol   = uc[0] * du[0] * du[1] * du[2];
    break;
  case CoordForm::SPHERICAL:
    uc[1] = 2.0 * M_PI * uc[1] - M_PI;
    du[1] = 2.0 * M_PI * du[1];
    uc[2] = M_PI * uc[2];
    du[2] = M_PI * du[2];
    vol   = (uc[0] * uc[0] + du[0] * du[0] / 12.0) * du[0] //r
        * du[1]                                            //theta
        * 2.0 * sin(uc[2]) * sin(.5 * du[2]);              //phi
    break;
  default:
    break;
  }
  volume_ = vol * det(axes_);
  return;
}

bool NESpaceGrid::checkAxisGridValues(const SpaceGridInput& input, const AxTensor& axes)
{
  auto& axis_labels = input.get_axis_labels();
  auto& axis_grids  = input.get_axis_grids();

  //check that all axis grid values fall in the allowed intervals
  std::map<std::string, int> cartmap;
  for (int d = 0; d < OHMMS_DIM; d++)
  {
    cartmap[std::string(SpaceGridInput::ax_cartesian.at(d))] = d;
  }
  bool succeeded = true;
  for (int d = 0; d < OHMMS_DIM; d++)
  {
    if (cartmap.find(axis_labels[d]) != cartmap.end())
    {
      if (axis_grids[d].umin < -1.0 || axis_grids[d].umax > 1.0)
      {
        app_log() << "  grid values for " << axis_labels[d] << " must fall in [-1,1]" << std::endl;
        app_log() << "  interval provided: [" << axis_grids[d].umin << "," << axis_grids[d].umax << "]" << std::endl;
        succeeded = false;
      }
    }
    else if (axis_labels[d] == "phi")
    {
      if (std::abs(axis_grids[d].umin) + std::abs(axis_grids[d].umax) > 1.0)
      {
        app_log() << "  phi interval cannot be longer than 1" << std::endl;
        app_log() << "  interval length provided: " << std::abs(axis_grids[d].umin) + std::abs(axis_grids[d].umax)
                  << std::endl;
        succeeded = false;
      }
    }
    else
    {
      if (axis_grids[d].umin < 0.0 || axis_grids[d].umax > 1.0)
      {
        app_log() << "  grid values for " << axis_labels[d] << " must fall in [0,1]" << std::endl;
        app_log() << "  interval provided: [" << axis_grids[d].umin << "," << axis_grids[d].umax << "]" << std::endl;
        succeeded = false;
      }
    }
  }
  return succeeded;
}

void NESpaceGrid::write_description(std::ostream& os, std::string& indent)
{
  os << indent + "SpaceGrid" << std::endl;
  std::string s;
  using CoordForm = SpaceGridInput::CoordForm;
  switch (input_.get_coord_form())
  {
  case CoordForm::CARTESIAN:
    s = "cartesian";
    break;
  case CoordForm::CYLINDRICAL:
    s = "cylindrical";
    break;
  case CoordForm::SPHERICAL:
    s = "spherical";
    break;
  default:
    break;
  }
  auto& agr         = input_.get_axis_grids();
  auto& axis_labels = input_.get_axis_labels();
  os << indent + "  SpaceGridInput::lookup_input_ename_value(input_.get_coord_form()  = " + s << std::endl;
  os << indent + "  buffer_offset = " << buffer_offset_ << std::endl;
  os << indent + "  ndomains_      = " << ndomains_ << std::endl;
  os << indent + "  axes  = " << axes_ << std::endl;
  os << indent + "  axinv = " << axinv_ << std::endl;
  for (int d = 0; d < OHMMS_DIM; d++)
  {
    os << indent + "  axis " << axis_labels[d] << ":" << std::endl;
    os << indent + "    umin = " << agr[d].umin << std::endl;
    os << indent + "    umax = " << agr[d].umax << std::endl;
    os << indent + "    du   = " << 1.0 / agr[d].odu << std::endl;
    os << indent + "    dm   = " << dm_[d] << std::endl;
    os << indent + "    gmap = ";
    for (int g = 0; g < gmap_[d].size(); g++)
    {
      os << gmap_[d][g] << " ";
    }
    os << std::endl;
  }
  os << indent + "end NESpaceGrid" << std::endl;
}

int NESpaceGrid::allocate_buffer_space(BufferType& buf)
{
  buffer_offset_ = buf.size();
  std::vector<Real> tmp(nvalues_per_domain_ * ndomains_);
  buf.add(tmp.begin(), tmp.end());
  buffer_start_ = buffer_offset_;
  buffer_end_   = buffer_start_ + nvalues_per_domain_ * ndomains_ - 1;
  return buffer_offset_;
}

void NESpaceGrid::registerCollectables(hdf_archive& file, int grid_index) const
{
  using iMatrix = Matrix<int>;
  iMatrix imat;
  std::vector<int> ng(1);
  int cshift = 1;
  std::stringstream ss;
  ss << grid_index + cshift;
  std::string name = "spacegrid" + ss.str();
  h5desc.emplace_back(name);
  auto& oh = h5desc.back();
  ng[0] = nvalues_per_domain_ * ndomains_;
  oh.set_dimensions(ng, buffer_offset_);
  oh.open(gid);
  int coord = (int)coordinate;
  oh.addProperty(const_cast<int&>(coord), "coordinate");
  oh.addProperty(const_cast<int&>(ndomains_), "ndomains");
  oh.addProperty(const_cast<int&>(nvalues_per_domain_), "nvalues_per_domain_");
  oh.addProperty(const_cast<Real&>(volume), "volume");
  oh.addProperty(const_cast<Matrix<Real>&>(domain_volumes), "domain_volumes");
  oh.addProperty(const_cast<Matrix<Real>&>(domain_centers), "domain_centers");
  if (coordinate != voronoi)
  {
    oh.addProperty(const_cast<Point&>(origin), "origin");
    oh.addProperty(const_cast<Tensor_t&>(axes), "axes");
    oh.addProperty(const_cast<Tensor_t&>(axinv), "axinv");
    oh.addProperty(const_cast<Matrix<Real>&>(domain_uwidths), "domain_uwidths");
    //add dimensioned quantities
    std::map<std::string, int> axtmap;
    axtmap["x"]     = 0;
    axtmap["y"]     = 1;
    axtmap["z"]     = 2;
    axtmap["r"]     = 3;
    axtmap["phi"]   = 4;
    axtmap["theta"] = 5;
    int axtypes[DIM];
    for (int d = 0; d < OHMMS_DIM; d++)
    {
      axtypes[d] = axtmap[axis_labels[d]];
    }
    int n;
    const int ni = 3;
    int* ivar[ni];
    std::string iname[ni];
    n        = 0;
    ivar[n]  = (int*)axtypes;
    iname[n] = "axtypes";
    n++;
    ivar[n]  = (int*)dimensions;
    iname[n] = "dimensions";
    n++;
    ivar[n]  = (int*)dm;
    iname[n] = "dm";
    n++;
    const int nr = 3;
    Real* rvar[nr];
    std::string rname[nr];
    n        = 0;
    rvar[n]  = (Real*)odu;
    rname[n] = "odu";
    n++;
    rvar[n]  = (Real*)umin;
    rname[n] = "umin";
    n++;
    rvar[n]  = (Real*)umax;
    rname[n] = "umax";
    n++;
    imat.resize(DIM, 1);
    for (int i = 0; i < ni; i++)
    {
      for (int d = 0; d < OHMMS_DIM; d++)
        imat(d, 0) = ivar[i][d];
      oh.addProperty(const_cast<iMatrix&>(imat), iname[i]);
    }
    Matrix<Real> rmat;
    rmat.resize(DIM, 1);
    for (int i = 0; i < nr; i++)
    {
      for (int d = 0; d < OHMMS_DIM; d++)
        rmat(d, 0) = rvar[i][d];
      oh.addProperty(const_cast<Matrix<Real>&>(rmat), rname[i]);
    }
    for (int d = 0; d < OHMMS_DIM; d++)
    {
      int gsize = gmap[d].size();
      imat.resize(gsize, 1);
      for (int i = 0; i < gsize; i++)
      {
        int gval   = gmap[d][i];
        imat(i, 0) = gval;
      }
      int ival           = d + 1;
      std::string gmname = "gmap" + int2string(ival);
      oh.addProperty(const_cast<iMatrix&>(imat), gmname);
    }
  }

  return;
}


#define NESpaceGrid_CHECK


void NESpaceGrid::evaluate(const ParticlePos& R,
                           const Matrix<Real>& values,
                           BufferType& buf,
                           std::vector<bool>& particles_outside,
                           const DistanceTableAB& dtab)
{
  int p, v;
  int nparticles = values.size1();
  int nvalues    = values.size2();
  int iu[DIM];
  int buf_index;
  const Real o2pi = 1.0 / (2.0 * M_PI);
  if (!chempot)
  {
    switch (coordinate)
    {
    case cartesian:
      if (periodic)
      {
        for (p = 0; p < nparticles; p++)
        {
          particles_outside[p] = false;
          u                    = dot(axinv, (R[p] - origin));
          for (int d = 0; d < OHMMS_DIM; ++d)
            iu[d] = gmap[d][floor((u[d] - umin[d]) * odu[d])];
          buf_index = buffer_offset;
          for (int d = 0; d < OHMMS_DIM; ++d)
            buf_index += nvalues * dm[d] * iu[d];
          for (v = 0; v < nvalues; v++, buf_index++)
            buf[buf_index] += values(p, v);
        }
      }
      else
      {
        for (p = 0; p < nparticles; p++)
        {
          u = dot(axinv, (R[p] - origin));
          if (u[0] > umin[0] && u[0] < umax[0] && u[1] > umin[1] && u[1] < umax[1] && u[2] > umin[2] && u[2] < umax[2])
          {
            particles_outside[p] = false;
            iu[0]                = gmap[0][floor((u[0] - umin[0]) * odu[0])];
            iu[1]                = gmap[1][floor((u[1] - umin[1]) * odu[1])];
            iu[2]                = gmap[2][floor((u[2] - umin[2]) * odu[2])];
            buf_index            = buffer_offset + nvalues * (dm[0] * iu[0] + dm[1] * iu[1] + dm[2] * iu[2]);
            for (v = 0; v < nvalues; v++, buf_index++)
              buf[buf_index] += values(p, v);
          }
        }
      }
      break;
    case cylindrical:
      for (p = 0; p < nparticles; p++)
      {
        ub   = dot(axinv, (R[p] - origin));
        u[0] = sqrt(ub[0] * ub[0] + ub[1] * ub[1]);
        u[1] = atan2(ub[1], ub[0]) * o2pi + .5;
        u[2] = ub[2];
        if (u[0] > umin[0] && u[0] < umax[0] && u[1] > umin[1] && u[1] < umax[1] && u[2] > umin[2] && u[2] < umax[2])
        {
          particles_outside[p] = false;
          iu[0]                = gmap[0][floor((u[0] - umin[0]) * odu[0])];
          iu[1]                = gmap[1][floor((u[1] - umin[1]) * odu[1])];
          iu[2]                = gmap[2][floor((u[2] - umin[2]) * odu[2])];
          buf_index            = buffer_offset + nvalues * (dm[0] * iu[0] + dm[1] * iu[1] + dm[2] * iu[2]);
          for (v = 0; v < nvalues; v++, buf_index++)
            buf[buf_index] += values(p, v);
        }
      }
      break;
    case spherical:
      for (p = 0; p < nparticles; p++)
      {
        ub   = dot(axinv, (R[p] - origin));
        u[0] = sqrt(ub[0] * ub[0] + ub[1] * ub[1] + ub[2] * ub[2]);
        u[1] = atan2(ub[1], ub[0]) * o2pi + .5;
        u[2] = acos(ub[2] / u[0]) * o2pi * 2.0;
        if (u[0] > umin[0] && u[0] < umax[0] && u[1] > umin[1] && u[1] < umax[1] && u[2] > umin[2] && u[2] < umax[2])
        {
          particles_outside[p] = false;
          iu[0]                = gmap[0][floor((u[0] - umin[0]) * odu[0])];
          iu[1]                = gmap[1][floor((u[1] - umin[1]) * odu[1])];
          iu[2]                = gmap[2][floor((u[2] - umin[2]) * odu[2])];
          buf_index            = buffer_offset + nvalues * (dm[0] * iu[0] + dm[1] * iu[1] + dm[2] * iu[2]);
          for (v = 0; v < nvalues; v++, buf_index++)
            buf[buf_index] += values(p, v);
        }
      }
      break;
    case voronoi:
      //find cell center nearest to each dynamic particle
      int nd, nn;
      Real dist;
      for (p = 0; p < ndparticles_; p++)
      {
        const auto& dist = dtab.getDistRow(p);
        for (nd = 0; nd < ndomains_; nd++)
          if (dist[nd] < nearcell[p].r)
          {
            nearcell[p].r = dist[nd];
            nearcell[p].i = nd;
          }
      }
      //accumulate values for each dynamic particle
      for (p = 0; p < ndparticles_; p++)
      {
        buf_index = buffer_offset + nvalues * nearcell[p].i;
        for (v = 0; v < nvalues; v++, buf_index++)
          buf[buf_index] += values(p, v);
      }
      //accumulate values for static particles (static particles == cell centers)
      buf_index = buffer_offset;
      for (p = ndparticles_; p < nparticles; p++)
        for (v = 0; v < nvalues; v++, buf_index++)
          buf[buf_index] += values(p, v);
      //each particle belongs to some voronoi cell
      for (p = 0; p < nparticles; p++)
        particles_outside[p] = false;
      //reset distances
      for (p = 0; p < ndparticles_; p++)
        nearcell[p].r = std::numeric_limits<Real>::max();
      break;
    default:
      app_log() << "  coordinate type must be cartesian, cylindrical, spherical, or voronoi" << std::endl;
      APP_ABORT("SpaceGrid::evaluate");
    }
  }
  else
  //chempot: sort values by particle count in volumes
  {
    int cell_index;
    int nd;
    std::fill(cellsamples.begin(), cellsamples.end(), 0.0);
    switch (coordinate)
    {
    case cartesian:
      for (p = 0; p < nparticles; p++)
      {
        u = dot(axinv, (R[p] - origin));
        if (u[0] > umin[0] && u[0] < umax[0] && u[1] > umin[1] && u[1] < umax[1] && u[2] > umin[2] && u[2] < umax[2])
        {
          particles_outside[p] = false;
          iu[0]                = gmap[0][floor((u[0] - umin[0]) * odu[0])];
          iu[1]                = gmap[1][floor((u[1] - umin[1]) * odu[1])];
          iu[2]                = gmap[2][floor((u[2] - umin[2]) * odu[2])];
          cell_index           = dm[0] * iu[0] + dm[1] * iu[1] + dm[2] * iu[2];
          for (v = 0; v < nvalues; v++)
            cellsamples(cell_index, v) += values(p, v);
          cellsamples(cell_index, nvalues) += 1.0;
        }
      }
      break;
    case cylindrical:
      for (p = 0; p < nparticles; p++)
      {
        ub   = dot(axinv, (R[p] - origin));
        u[0] = sqrt(ub[0] * ub[0] + ub[1] * ub[1]);
        u[1] = atan2(ub[1], ub[0]) * o2pi + .5;
        u[2] = ub[2];
        if (u[0] > umin[0] && u[0] < umax[0] && u[1] > umin[1] && u[1] < umax[1] && u[2] > umin[2] && u[2] < umax[2])
        {
          particles_outside[p] = false;
          iu[0]                = gmap[0][floor((u[0] - umin[0]) * odu[0])];
          iu[1]                = gmap[1][floor((u[1] - umin[1]) * odu[1])];
          iu[2]                = gmap[2][floor((u[2] - umin[2]) * odu[2])];
          cell_index           = dm[0] * iu[0] + dm[1] * iu[1] + dm[2] * iu[2];
          for (v = 0; v < nvalues; v++)
            cellsamples(cell_index, v) += values(p, v);
          cellsamples(cell_index, nvalues) += 1.0;
        }
      }
      break;
    case spherical:
      for (p = 0; p < nparticles; p++)
      {
        ub   = dot(axinv, (R[p] - origin));
        u[0] = sqrt(ub[0] * ub[0] + ub[1] * ub[1] + ub[2] * ub[2]);
        u[1] = atan2(ub[1], ub[0]) * o2pi + .5;
        u[2] = acos(ub[2] / u[0]) * o2pi * 2.0;
        if (u[0] > umin[0] && u[0] < umax[0] && u[1] > umin[1] && u[1] < umax[1] && u[2] > umin[2] && u[2] < umax[2])
        {
          particles_outside[p] = false;
          iu[0]                = gmap[0][floor((u[0] - umin[0]) * odu[0])];
          iu[1]                = gmap[1][floor((u[1] - umin[1]) * odu[1])];
          iu[2]                = gmap[2][floor((u[2] - umin[2]) * odu[2])];
          cell_index           = dm[0] * iu[0] + dm[1] * iu[1] + dm[2] * iu[2];
          for (v = 0; v < nvalues; v++)
            cellsamples(cell_index, v) += values(p, v);
          cellsamples(cell_index, nvalues) += 1.0;
        }
      }
      break;
    case voronoi:
      //find cell center nearest to each dynamic particle
      int nn;
      Real dist;
      APP_ABORT("SoA transformation needed for Voronoi grids")
      //for (nd = 0; nd < ndomains_; nd++)
      //  for (nn = dtab.M[nd], p = 0; nn < dtab.M[nd + 1]; ++nn, ++p)
      //  {
      //    dist = dtab.r(nn);
      //    if (dist < nearcell[p].r)
      //    {
      //      nearcell[p].r = dist;
      //      nearcell[p].i = nd;
      //    }
      //  }
      //accumulate values for each dynamic particle
      for (p = 0; p < ndparticles_; p++)
      {
        cell_index = nearcell[p].i;
        for (v = 0; v < nvalues; v++)
          cellsamples(cell_index, v) += values(p, v);
        cellsamples(cell_index, nvalues) += 1.0;
      }
      //accumulate values for static particles (static particles == cell centers)
      for (p = ndparticles_, cell_index = 0; p < nparticles; p++, cell_index++)
      {
        for (v = 0; v < nvalues; v++)
          cellsamples(cell_index, v) += values(p, v);
        cellsamples(cell_index, nvalues) += 1.0;
      }
      //each particle belongs to some voronoi cell
      for (p = 0; p < nparticles; p++)
        particles_outside[p] = false;
      //reset distances
      for (p = 0; p < ndparticles_; p++)
        nearcell[p].r = std::numeric_limits<Real>::max();
      break;
    default:
      app_log() << "  coordinate type must be cartesian, cylindrical, spherical, or voronoi" << std::endl;
      APP_ABORT("SpaceGrid::evaluate");
    }
    //now place samples in the buffer according to how
    // many particles are in each cell
    int nincell;
    buf_index = buffer_offset;
    for (nd = 0; nd < ndomains_; nd++)
    {
      nincell = cellsamples(nd, nvalues) - reference_count[nd];
      if (nincell >= npmin && nincell <= npmax)
      {
        buf_index = buffer_offset + (nd * npvalues + nincell - npmin) * nvalues;
        for (v = 0; v < nvalues; v++, buf_index++)
          buf[buf_index] += cellsamples(nd, v);
      }
    }
  }
}


void NESpaceGrid::sum(const BufferType& buf, Real* vals)
{
  for (int v = 0; v < nvalues_per_domain_; v++)
  {
    vals[v] = 0.0;
  }
  for (int i = 0, n = buffer_offset; i < ndomains_; i++, n += nvalues_per_domain_)
  {
    for (int v = 0; v < nvalues_per_domain_; v++)
    {
      vals[v] += buf[n + v];
    }
  }
}


bool NESpaceGrid::check_grid(void)
{
  app_log() << "SpaceGrid::check_grid" << std::endl;
  const Real o2pi = 1.0 / (2.0 * M_PI);
  int iu[DIM];
  int idomain;
  bool ok = true;
  Point dc;
  for (int i = 0; i < ndomains_; i++)
  {
    for (int d = 0; d < OHMMS_DIM; d++)
      dc[d] = domain_centers(i, d);
    ub = dot(axinv, (dc - origin));
    switch (coordinate)
    {
    case cartesian:
      u = ub;
      break;
    case cylindrical:
      u[0] = sqrt(ub[0] * ub[0] + ub[1] * ub[1]);
      u[1] = atan2(ub[1], ub[0]) * o2pi + .5;
      u[2] = ub[2];
      break;
    case spherical:
      u[0] = sqrt(ub[0] * ub[0] + ub[1] * ub[1] + ub[2] * ub[2]);
      u[1] = atan2(ub[1], ub[0]) * o2pi + .5;
      u[2] = acos(ub[2] / u[0]) * o2pi * 2.0;
      break;
    default:
      break;
    }
    iu[0]   = gmap[0][floor((u[0] - umin[0]) * odu[0])];
    iu[1]   = gmap[1][floor((u[1] - umin[1]) * odu[1])];
    iu[2]   = gmap[2][floor((u[2] - umin[2]) * odu[2])];
    idomain = dm[0] * iu[0] + dm[1] * iu[1] + dm[2] * iu[2];
    if (idomain != i)
    {
      app_log() << "  cell mismatch " << i << " " << idomain << std::endl;
      ok = false;
    }
  }
  if (!ok)
  {
    app_log() << "  NESpaceGrid cells do not map onto themselves" << std::endl;
  }
  app_log() << "end NESpaceGrid::check_grid" << std::endl;
  return ok;
}

} // namespace qmcplusplus
