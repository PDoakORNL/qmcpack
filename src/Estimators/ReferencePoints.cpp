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
// Refactored from:  QMCHamiltonians/ReferencePoints.cpp
//////////////////////////////////////////////////////////////////////////////////////


#include "ReferencePoints.h"

namespace qmcplusplus
{

ReferencePoints::ReferencePoints(ReferencePointsInput&& input, ParticleSet& P, std::vector<ParticleSet*>& Psets)
    : input_(input)
{
  using ReferencePointsInput::CoordForm;

  points_["zero"] = 0 * P.getLattice().a(0);
  points_["a1"]   = P.getLattice().a(0);
  points_["a2"]   = P.getLattice().a(1);
  points_["a3"]   = P.getLattice().a(2);
  //points_["center"]= .5*(P.getLattice().a(0)+P.getLattice().a(1)+P.Lattice.a(2))
  //set points_ on face centers
  points_["f1p"] = points_["zero"] + .5 * points_["a1"];
  points_["f1m"] = points_["zero"] - .5 * points_["a1"];
  points_["f2p"] = points_["zero"] + .5 * points_["a2"];
  points_["f2m"] = points_["zero"] - .5 * points_["a2"];
  points_["f3p"] = points_["zero"] + .5 * points_["a3"];
  points_["f3m"] = points_["zero"] - .5 * points_["a3"];
  //set points_ on cell corners
  points_["cmmm"] = points_["zero"] + .5 * (-1 * points_["a1"] - points_["a2"] - points_["a3"]);
  points_["cpmm"] = points_["zero"] + .5 * (points_["a1"] - points_["a2"] - points_["a3"]);
  points_["cmpm"] = points_["zero"] + .5 * (-1 * points_["a1"] + points_["a2"] - points_["a3"]);
  points_["cmmp"] = points_["zero"] + .5 * (-1 * points_["a1"] - points_["a2"] + points_["a3"]);
  points_["cmpp"] = points_["zero"] + .5 * (-1 * points_["a1"] + points_["a2"] + points_["a3"]);
  points_["cpmp"] = points_["zero"] + .5 * (points_["a1"] - points_["a2"] + points_["a3"]);
  points_["cppm"] = points_["zero"] + .5 * (points_["a1"] + points_["a2"] - points_["a3"]);
  points_["cppp"] = points_["zero"] + .5 * (points_["a1"] + points_["a2"] + points_["a3"]);
  //get points from requested particle sets
  int cshift = 1;
  for (int i = 0; i < Psets.size(); i++)
  {
    ParticleSet& PS = *Psets[i];
    for (int p = 0; p < PS.getTotalNum(); p++)
    {
      std::stringstream ss;
      ss << p + cshift;
      points_[PS.getName() + ss.str()] = PS.R[p];
    }
  }

  for (int i = 0; i < OHMMS_DIM; i++)
    for (int d = 0; d < OHMMS_DIM; d++)
      axes_(d, i) = P.getLattice().a(i)[d];
  Tensor<Real, OHMMSDIM> crd;
  switch (input_.get_coord())
  {
  case CoordForm::CELL:
    crd = axes_;
    break;
  case CoordForm::CARTESIAN:
    for (int i = 0; i < OHMMSDIM; i++)
      for (int d = 0; d < OHMMSDIM; d++)
        if (d == i)
          crd(i, i) = 1.0;
        else
          crd(d, i) = 0.0;
    break;
  }

  Points input_points = input_.get_points();
  for(auto it = input_points.begin(); it!= input_points.end; ++ it) {
    point& rp = it->second;
    rp                 = dot(crd, rp);
  }
  std::move(input_points.begin(), input_points.end(), std::back_inserter(points_));
}

void ReferencePoints::write_description(std::ostream& os, std::string& indent)
{
  os << indent + "reference_points" << std::endl;
  std::map<std::string, Point>::const_iterator it, end = points.end();
  for (it = points.begin(); it != end; ++it)
  {
    os << indent + "  " << it->first << ": " << it->second << std::endl;
  }
  os << indent + "end reference_points" << std::endl;
  return;
}

void ReferencePoints::save(std::vector<ObservableHelper>& h5desc, hid_t gid) const
{
  h5desc.emplace_back("reference_points");
  auto& oh = h5desc.back();
  oh.open(gid);
  std::map<std::string, Point>::const_iterator it;
  for (it = points.begin(); it != points.end(); ++it)
  {
    oh.addProperty(const_cast<Point&>(it->second), it->first);
  }
  return;
}


} // namespace qmcplusplus
