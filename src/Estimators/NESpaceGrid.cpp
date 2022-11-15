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


NESpaceGrid::NESpaceGrid(SpaceGridInput&& sgi, NERefernecePoints::Points& points, ) : input_(std::move(sgi)) {}

NESpaceGrid::NESpaceGrid(int& nvalues)
{
  nvalues_per_domain = nvalues;
  Rptcl              = 0;
  ndparticles        = 0;
  periodic           = false;
}

NESpaceGrid::NESpaceGrid(SpaceGridInput& sgi, NEReferencePoints::Points& points, int ndp, bool is_periodic)
    : input_(sgi), ndparticles_(ndp), periodic_(is_periodic),
{
  ndomains = -1;
  init_success = initialize_rectilinear(cur, input_.get_space_grid_coord(), points);
}


void SpaceGrid::processAxis()
{
    {
      naxes++;
      if (cmap.find(label) == cmap.end())
      {
        app_log() << "      grid label " << label << " is invalid for " << coord << " coordinates" << std::endl;
        app_log() << "      valid options are: ";
        for (int d = 0; d < DIM; d++)
          app_log() << axlabel_[d] << ", ";
        app_log() << std::endl;
        APP_ABORT("SpaceGrid::put");
      }
      iaxis = cmap[label];
      if (scale == undefined)
        frac = 1.0;
      else
        frac = string2real(scale);
      for (int d = 0; d < DIM; d++)
        axes(d, iaxis) = frac * (points[p1][d] - points[p2][d]);
      //read in the grid contents
      //  remove spaces inside of parentheses
      std::string gtmp;
      bool inparen = false;
      char gc;
      for (int i = 0; i < grid.size(); i++)
      {
        gc = grid[i];
        if (gc == '(')
        {
          inparen = true;
          gtmp += ' ';
        }
        if (!(inparen && gc == ' '))
          gtmp += gc;
        if (gc == ')')
        {
          inparen = false;
          gtmp += ' ';
        }
      }
      grid = gtmp;
      //  break into tokens
      std::vector<std::string> tokens = split(grid, " ");
      //  count the number of intervals
      int nintervals = 0;
      for (int i = 0; i < tokens.size(); i++)
        if (tokens[i][0] != '(')
          nintervals++;
      nintervals--;
      //  allocate temporary interval variables
      std::vector<int> ndom_int, ndu_int;
      std::vector<Real> du_int;
      ndom_int.resize(nintervals);
      du_int.resize(nintervals);
      ndu_int.resize(nintervals);
      //  determine number of domains in each interval and the width of each domain
      Real u1     = string2real(tokens[0]);
      umin[iaxis] = u1;
      if (std::abs(u1) > 1.0000001)
      {
        app_log() << "  interval endparticles cannot be greater than " << 1 << std::endl;
        app_log() << "  endpoint provided: " << u1 << std::endl;
        succeeded = false;
      }
      bool is_int        = false;
      bool has_paren_val = false;
      Real du_i;
      int ndom_i   = 1;
      int interval = -1;
      for (int i = 1; i < tokens.size(); i++)
      {
        if (tokens[i][0] != '(')
        {
          Real u2     = string2real(tokens[i]);
          umax[iaxis] = u2;
          if (!has_paren_val)
            du_i = u2 - u1;
          has_paren_val = false;
          interval++;
          if (u2 < u1)
          {
            app_log() << "  interval (" << u1 << "," << u2 << ") is negative" << std::endl;
            succeeded = false;
          }
          if (std::abs(u2) > 1.0000001)
          {
            app_log() << "  interval endparticles cannot be greater than " << 1 << std::endl;
            app_log() << "  endpoint provided: " << u2 << std::endl;
            succeeded = false;
          }
          if (is_int)
          {
            du_int[interval]   = (u2 - u1) / ndom_i;
            ndom_int[interval] = ndom_i;
          }
          else
          {
            du_int[interval]   = du_i;
            ndom_int[interval] = floor((u2 - u1) / du_i + .5);
            if (std::abs(u2 - u1 - du_i * ndom_int[interval]) > utol)
            {
              app_log() << "  interval (" << u1 << "," << u2 << ") not divisible by du=" << du_i << std::endl;
              succeeded = false;
            }
          }
          u1 = u2;
        }
        else
        {
          has_paren_val         = true;
          std::string paren_val = tokens[i].substr(1, tokens[i].length() - 2);
          is_int                = tokens[i].find(".") == std::string::npos;
          if (is_int)
          {
            ndom_i = string2int(paren_val);
            du_i   = -1.0;
          }
          else
          {
            ndom_i = 0;
            du_i   = string2real(paren_val);
          }
        }
      }
      // find the smallest domain width
      Real du_min = 1.0;
      for (int i = 0; i < du_int.size(); i++)
        du_min = std::min(du_min, du_int[i]);
      odu[iaxis] = 1.0 / du_min;
      // make sure it divides into all other domain widths
      for (int i = 0; i < du_int.size(); i++)
      {
        ndu_int[i] = floor(du_int[i] / du_min + .5);
        if (std::abs(du_int[i] - ndu_int[i] * du_min) > utol)
        {
          app_log() << "interval " << i + 1 << " of axis " << iaxis + 1 << " is not divisible by smallest subinterval "
                    << du_min << std::endl;
          succeeded = false;
        }
      }
      // set up the interval map such that gmap[u/du]==domain index
      gmap[iaxis].resize(floor((umax[iaxis] - umin[iaxis]) * odu[iaxis] + .5));
      int n  = 0;
      int nd = -1;
      for (int i = 0; i < ndom_int.size(); i++)
        for (int j = 0; j < ndom_int[i]; j++)
        {
          nd++;
          for (int k = 0; k < ndu_int[i]; k++)
          {
            //app_log()<<"        accessing gmap "<<iaxis<<" "<<n<< std::endl;
            gmap[iaxis][n] = nd;
            n++;
          }
        }
      dimensions[iaxis] = nd + 1;
      //end read in the grid contents
      //save interval width information
      int ndom_tot = 0;
      for (int i = 0; i < ndom_int.size(); i++)
        ndom_tot += ndom_int[i];
      ndu_per_interval[iaxis].resize(ndom_tot);
      int idom = 0;
      for (int i = 0; i < ndom_int.size(); i++)
        for (int ii = 0; ii < ndom_int[i]; ii++)
        {
          ndu_per_interval[iaxis][idom] = ndu_int[i];
          idom++;
        }
}

bool SpaceGrid::initialize_rectilinear(SpaceGridInput& input, Points& points)
{
  // This code should be refactored to SpaceGridInput such that a simple map of
  // axis is available.
  using CoordForm = typename SpaceGridInput::CoordForm;
  std::map<std::string, int> cmap;
  coord_form = input_.get_coordinate_form();
  if (coord_form == CoordForm::CYLINDRICAL || coord_form == CoordForm::SPHERICAL)
    periodic_ = false;

  int iaxis             = 0;
  int naxes             = 0;
  std::string& origin_p1 = input_.get_origin_p1();
  if (origin_p1.size() > 0)
  {
    std::string& origin_p2 = input_.get_origin_p2();
    origin = points_[p1] + input_.get_origin_fraction() * (points[origin_p2] - points[origin_p1]);
  }
  else
    origin_                = points_["zero"];
  
  // variables for loop
  Real utol = 1e-5;
  std::string grid;
  std::vector<int> ndu_per_interval[DIM];
  while (element != NULL)
  {
    std::string name      = (const char*)element->name;
    OhmmsAttributeSet* ea = new OhmmsAttributeSet;
    std::string p1        = undefined;
    std::string p2        = "zero";
    std::string fraction  = undefined;
    Real frac             = -1.0;
    std::string scale     = undefined;
    std::string label     = undefined;
    grid                  = "0 1";
  
  axinv = inverse(axes);
  //check that all axis grid values fall in the allowed intervals
  std::map<std::string, int> cartmap;
  for (int d = 0; d < DIM; d++)
  {
    cartmap[ax_cartesian[d]] = d;
  }
  for (int d = 0; d < DIM; d++)
  {
    if (cartmap.find(axlabel_[d]) != cartmap.end())
    {
      if (umin[d] < -1.0 || umax[d] > 1.0)
      {
        app_log() << "  grid values for " << axlabel_[d] << " must fall in [-1,1]" << std::endl;
        app_log() << "  interval provided: [" << umin[d] << "," << umax[d] << "]" << std::endl;
        succeeded = false;
      }
    }
    else if (axlabel_[d] == "phi")
    {
      if (std::abs(umin[d]) + std::abs(umax[d]) > 1.0)
      {
        app_log() << "  phi interval cannot be longer than 1" << std::endl;
        app_log() << "  interval length provided: " << std::abs(umin[d]) + std::abs(umax[d]) << std::endl;
        succeeded = false;
      }
    }
    else
    {
      if (umin[d] < 0.0 || umax[d] > 1.0)
      {
        app_log() << "  grid values for " << axlabel_[d] << " must fall in [0,1]" << std::endl;
        app_log() << "  interval provided: [" << umin[d] << "," << umax[d] << "]" << std::endl;
        succeeded = false;
      }
    }
  }
  //set grid dimensions
  // C/Python style indexing
  dm[0] = dimensions[1] * dimensions[2];
  dm[1] = dimensions[2];
  dm[2] = 1;
  // Fortran/Matlab style indexing
  //dm[0] = 1;
  //dm[1] = dimensions[0];
  //dm[2] = dimensions[0]*dimensions[1];
  ndomains = dimensions[0] * dimensions[1] * dimensions[2];
  volume   = std::abs(det(axes)) * 8.0; //axes span only one octant
  //compute domain volumes, centers, and widths
  domain_volumes.resize(ndomains, 1);
  domain_centers.resize(ndomains, DIM);
  domain_uwidths.resize(ndomains, DIM);
  std::vector<Real> interval_centers[DIM];
  std::vector<Real> interval_widths[DIM];
  for (int d = 0; d < DIM; d++)
  {
    int nintervals = ndu_per_interval[d].size();
    app_log() << "nintervals " << d << " " << nintervals << std::endl;
    interval_centers[d].resize(nintervals);
    interval_widths[d].resize(nintervals);
    interval_widths[d][0]  = ndu_per_interval[d][0] / odu[d];
    interval_centers[d][0] = interval_widths[d][0] / 2.0 + umin[d];
    for (int i = 1; i < nintervals; i++)
    {
      interval_widths[d][i]  = ndu_per_interval[d][i] / odu[d];
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
  Real vol    = 0.0;
  Real vscale = std::abs(det(axes));
  for (int i = 0; i < dimensions[0]; i++)
  {
    for (int j = 0; j < dimensions[1]; j++)
    {
      for (int k = 0; k < dimensions[2]; k++)
      {
        int idomain = dm[0] * i + dm[1] * j + dm[2] * k;
        du[0]       = interval_widths[0][i];
        du[1]       = interval_widths[1][j];
        du[2]       = interval_widths[2][k];
        uc[0]       = interval_centers[0][i];
        uc[1]       = interval_centers[1][j];
        uc[2]       = interval_centers[2][k];
        switch (coordinate)
        {
        case (cartesian):
          vol = du[0] * du[1] * du[2];
          ubc = uc;
          break;
        case (cylindrical):
          uc[1]  = 2.0 * M_PI * uc[1] - M_PI;
          du[1]  = 2.0 * M_PI * du[1];
          vol    = uc[0] * du[0] * du[1] * du[2];
          ubc[0] = uc[0] * cos(uc[1]);
          ubc[1] = uc[0] * sin(uc[1]);
          ubc[2] = uc[2];
          break;
        case (spherical):
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
        rc = dot(axes, ubc) + origin;
        //app_log()<< std::endl;
        //app_log()<<"umin "<<uc-du/2<< std::endl;
        //app_log()<<"uc "<<uc<< std::endl;
        //app_log()<<"umax "<<uc+du/2<< std::endl;
        //app_log()<<"rc "<<rc<< std::endl;
        domain_volumes(idomain, 0) = vol;
        for (int d = 0; d < DIM; d++)
        {
          domain_uwidths(idomain, d) = du[d];
          domain_centers(idomain, d) = rc[d];
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
  for (int d = 0; d < DIM; d++)
  {
    du[d] = umax[d] - umin[d];
    uc[d] = .5 * (umax[d] + umin[d]);
  }
  switch (coordinate)
  {
  case (cartesian):
    vol = du[0] * du[1] * du[2];
    break;
  case (cylindrical):
    uc[1] = 2.0 * M_PI * uc[1] - M_PI;
    du[1] = 2.0 * M_PI * du[1];
    vol   = uc[0] * du[0] * du[1] * du[2];
    break;
  case (spherical):
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
  volume = vol * det(axes);
  if (chempot)
  {
    reference_count.resize(ndomains);
    //rectilinear is referenced to vacuum for now
    fill(reference_count.begin(), reference_count.end(), 0.0);
  }
  succeeded = succeeded && check_grid();
  return succeeded;
}


void SpaceGrid::write_description(std::ostream& os, std::string& indent)
{
  os << indent + "SpaceGrid" << std::endl;
  std::string s;
  switch (coordinate)
  {
  case cartesian:
    s = "cartesian";
    break;
  case cylindrical:
    s = "cylindrical";
    break;
  case spherical:
    s = "spherical";
    break;
  default:
    break;
  }
  os << indent + "  coordinates   = " + s << std::endl;
  os << indent + "  buffer_offset = " << buffer_offset << std::endl;
  os << indent + "  ndomains      = " << ndomains << std::endl;
  os << indent + "  axes  = " << axes << std::endl;
  os << indent + "  axinv = " << axinv << std::endl;
  for (int d = 0; d < DIM; d++)
  {
    os << indent + "  axis " << axlabel_[d] << ":" << std::endl;
    os << indent + "    umin = " << umin[d] << std::endl;
    os << indent + "    umax = " << umax[d] << std::endl;
    os << indent + "    du   = " << 1.0 / odu[d] << std::endl;
    os << indent + "    dm   = " << dm[d] << std::endl;
    os << indent + "    gmap = ";
    for (int g = 0; g < gmap[d].size(); g++)
    {
      os << gmap[d][g] << " ";
    }
    os << std::endl;
  }
  os << indent + "end SpaceGrid" << std::endl;
}


int SpaceGrid::allocate_buffer_space(BufferType& buf)
{
  buffer_offset = buf.size();
  if (!chempot)
  {
    std::vector<Real> tmp(nvalues_per_domain * ndomains);
    buf.add(tmp.begin(), tmp.end());
    buffer_start = buffer_offset;
    buffer_end   = buffer_start + nvalues_per_domain * ndomains - 1;
  }
  else
  {
    std::vector<Real> tmp(nvalues_per_domain * npvalues * ndomains);
    buf.add(tmp.begin(), tmp.end());
    buffer_start = buffer_offset;
    buffer_end   = buffer_start + nvalues_per_domain * npvalues * ndomains - 1;
  }
  return buffer_offset;
}


void SpaceGrid::registerCollectables(std::vector<ObservableHelper>& h5desc, hid_t gid, int grid_index) const
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
  if (!chempot)
    ng[0] = nvalues_per_domain * ndomains;
  else
    ng[0] = nvalues_per_domain * npvalues * ndomains;
  oh.set_dimensions(ng, buffer_offset);
  oh.open(gid);
  int coord = (int)coordinate;
  oh.addProperty(const_cast<int&>(coord), "coordinate");
  oh.addProperty(const_cast<int&>(ndomains), "ndomains");
  oh.addProperty(const_cast<int&>(nvalues_per_domain), "nvalues_per_domain");
  oh.addProperty(const_cast<Real&>(volume), "volume");
  oh.addProperty(const_cast<Matrix_t&>(domain_volumes), "domain_volumes");
  oh.addProperty(const_cast<Matrix_t&>(domain_centers), "domain_centers");
  if (chempot)
  {
    oh.addProperty(const_cast<int&>(npmin), "min_part");
    oh.addProperty(const_cast<int&>(npmax), "max_part");
    int ref = (int)reference;
    oh.addProperty(const_cast<int&>(ref), "reference");
    imat.resize(reference_count.size(), 1);
    for (int i = 0; i < reference_count.size(); i++)
      imat(i, 0) = reference_count[i];
    oh.addProperty(const_cast<iMatrix&>(imat), "reference_count");
  }
  if (coordinate != voronoi)
  {
    oh.addProperty(const_cast<Point&>(origin), "origin");
    oh.addProperty(const_cast<Tensor_t&>(axes), "axes");
    oh.addProperty(const_cast<Tensor_t&>(axinv), "axinv");
    oh.addProperty(const_cast<Matrix_t&>(domain_uwidths), "domain_uwidths");
    //add dimensioned quantities
    std::map<std::string, int> axtmap;
    axtmap["x"]     = 0;
    axtmap["y"]     = 1;
    axtmap["z"]     = 2;
    axtmap["r"]     = 3;
    axtmap["phi"]   = 4;
    axtmap["theta"] = 5;
    int axtypes[DIM];
    for (int d = 0; d < DIM; d++)
    {
      axtypes[d] = axtmap[axlabel_[d]];
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
      for (int d = 0; d < DIM; d++)
        imat(d, 0) = ivar[i][d];
      oh.addProperty(const_cast<iMatrix&>(imat), iname[i]);
    }
    Matrix_t rmat;
    rmat.resize(DIM, 1);
    for (int i = 0; i < nr; i++)
    {
      for (int d = 0; d < DIM; d++)
        rmat(d, 0) = rvar[i][d];
      oh.addProperty(const_cast<Matrix_t&>(rmat), rname[i]);
    }
    for (int d = 0; d < DIM; d++)
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


#define SPACEGRID_CHECK


void SpaceGrid::evaluate(const ParticlePos& R,
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
          for (int d = 0; d < DIM; ++d)
            iu[d] = gmap[d][floor((u[d] - umin[d]) * odu[d])];
          buf_index = buffer_offset;
          for (int d = 0; d < DIM; ++d)
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
      for (p = 0; p < ndparticles; p++)
      {
        const auto& dist = dtab.getDistRow(p);
        for (nd = 0; nd < ndomains; nd++)
          if (dist[nd] < nearcell[p].r)
          {
            nearcell[p].r = dist[nd];
            nearcell[p].i = nd;
          }
      }
      //accumulate values for each dynamic particle
      for (p = 0; p < ndparticles; p++)
      {
        buf_index = buffer_offset + nvalues * nearcell[p].i;
        for (v = 0; v < nvalues; v++, buf_index++)
          buf[buf_index] += values(p, v);
      }
      //accumulate values for static particles (static particles == cell centers)
      buf_index = buffer_offset;
      for (p = ndparticles; p < nparticles; p++)
        for (v = 0; v < nvalues; v++, buf_index++)
          buf[buf_index] += values(p, v);
      //each particle belongs to some voronoi cell
      for (p = 0; p < nparticles; p++)
        particles_outside[p] = false;
      //reset distances
      for (p = 0; p < ndparticles; p++)
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
      //for (nd = 0; nd < ndomains; nd++)
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
      for (p = 0; p < ndparticles; p++)
      {
        cell_index = nearcell[p].i;
        for (v = 0; v < nvalues; v++)
          cellsamples(cell_index, v) += values(p, v);
        cellsamples(cell_index, nvalues) += 1.0;
      }
      //accumulate values for static particles (static particles == cell centers)
      for (p = ndparticles, cell_index = 0; p < nparticles; p++, cell_index++)
      {
        for (v = 0; v < nvalues; v++)
          cellsamples(cell_index, v) += values(p, v);
        cellsamples(cell_index, nvalues) += 1.0;
      }
      //each particle belongs to some voronoi cell
      for (p = 0; p < nparticles; p++)
        particles_outside[p] = false;
      //reset distances
      for (p = 0; p < ndparticles; p++)
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
    for (nd = 0; nd < ndomains; nd++)
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


void SpaceGrid::sum(const BufferType& buf, Real* vals)
{
  for (int v = 0; v < nvalues_per_domain; v++)
  {
    vals[v] = 0.0;
  }
  for (int i = 0, n = buffer_offset; i < ndomains; i++, n += nvalues_per_domain)
  {
    for (int v = 0; v < nvalues_per_domain; v++)
    {
      vals[v] += buf[n + v];
    }
  }
}


bool SpaceGrid::check_grid(void)
{
  app_log() << "SpaceGrid::check_grid" << std::endl;
  const Real o2pi = 1.0 / (2.0 * M_PI);
  int iu[DIM];
  int idomain;
  bool ok = true;
  Point dc;
  for (int i = 0; i < ndomains; i++)
  {
    for (int d = 0; d < DIM; d++)
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
    app_log() << "  SpaceGrid cells do not map onto themselves" << std::endl;
  }
  app_log() << "end SpaceGrid::check_grid" << std::endl;
  return ok;
}

} // namespace qmcplusplus
