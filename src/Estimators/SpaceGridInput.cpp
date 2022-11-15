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
// Some code refactored from:  QMCHamiltonians/SpaceGrid.cpp
//////////////////////////////////////////////////////////////////////////////////////

#include "SpaceGridInput.h"

#include <sstream>

#include "EstimatorInput.h"
#include "Message/UniformCommunicateError.h"


namespace qmcplusplus
{

SpaceGridInput::SpaceGridInput(xmlNodePtr cur)
{
  input_section_.readXML(cur);
  auto setIfInInput = LAMBDA_setIfInInput;
  setIfInInput(coord_form_, "coord");
  // rip open the axes inputs we're guarateed they have the proper dimensionality already

  auto axes = input_section_.get<std::vector<std::any>>("axis");
  checkAxes(axes);
  for (int d = 0; d < OHMMS_DIM; d++)
  {
    auto* axis_input = std::any_cast<SpaceGridAxisInput>(&axes[d]);
    auto& input = axis_input->get_input();
    axis_label_[d] = input.get<std::string>("label");
    axis_p1_[d] = input.get<std::string>("p1");
    axis_grid_[d] = input.get<std::vector<double>>("grid");
    axis_scale_[d] = input.get<double>("scale");
  }
}


void SpaceGridInput::checkAxes(std::vector<std::any>& axes) {
  auto checkLabels = [&](const std::array<const std::string, 3>& ax_labels, auto& axes) {
    for (auto& axis : axes)
    {
      auto* axis_input = std::any_cast<SpaceGridAxisInput>(&axis);
      std::string axis_label = axis_input->get_input().template get<std::string>("label");
      auto result      = std::find(std::begin(ax_labels), std::end(ax_labels), axis_label);
      if (result == std::end(ax_labels))
        throw UniformCommunicateError(axis_label + " is not a valid label for coord form " +
                                      input_section_.get<std::string>("coord"));
    }
  };

  // many many checks belong here
  switch (coord_form_)
  {
  case (CoordForm::CARTESIAN):
    checkLabels(ax_cartesian, axes);
    break;
  case (CoordForm::CYLINDRICAL):
    checkLabels(ax_cylindrical, axes);
    break;
  case (CoordForm::SPHERICAL):
    checkLabels(ax_spherical, axes);
    break;
  case (CoordForm::VORONOI):
    throw UniformCommunicateError("VORONOI SpaceGrid coordinates are not yet supported!");
    break;
  }
}

std::any SpaceGridInput::SpaceGridInputSection::assignAnyEnum(const std::string& name) const
{
  return lookupAnyEnum(name, get<std::string>(name), lookup_input_enum_value);
}

std::any makeSpaceGridInput(xmlNodePtr cur, std::string& value_label)
{
  SpaceGridInput sgi{cur};
  value_label = "spacegrid";
  return sgi;
}

void SpaceGridInput::SpaceGridInputSection::checkParticularValidity()
{
  auto axes       = std::any_cast<std::vector<std::any>>(values_["axis"]);
  auto axis_count = axes.size();
  if (axis_count != OHMMS_DIM)
  {
    std::ostringstream error;
    error << "SpaceGrid input must contain " << OHMMS_DIM << " axes, " << axis_count << " found!";
    throw UniformCommunicateError(error.str());
  }
}

} // namespace qmcplusplus
