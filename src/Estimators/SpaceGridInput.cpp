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

#include "EstimatorInput.h"



SpaceGridInput::SpaceGridInput(xmlNodePtr cur)
{
  input_section_.readXML(cur);
  auto setIfInInput = LAMBDA_setIfInInput;
  setIfInInput(space_grid_coord_, "coord");
  setIfInInput(space_gird_

}

std::any SpaceGridInput::SpaceGridInputSection::assignAnyEnum(const std::string& name) const
{
  return lookupAnyEnum(name, get<std::string>(name), lookup_input_enum_value);
}

