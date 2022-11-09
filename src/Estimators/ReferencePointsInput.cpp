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
// Some code refactored from:  QMCHamiltonians/ReferencePoints.cpp
//////////////////////////////////////////////////////////////////////////////////////

#include "ReferencePointsInput.h"

#include <string_view>

#include "EstimatorInput.h"

namespace qmcplusplus
{

ReferencePointsInput::ReferencePointsInput(xmlNodePtr cur)
{
  input_section_.readXML(cur);
  auto setIfInInput = LAMBDA_setIfInInput;
  setIfInInput(coord_form_, "coord");
  if (!readXML(cur))
    throw UniformCommunicateError("Failure in parsing reference points node content");
}

ReferencePointsInput::ReferencePointsInput(const Points& points, const CoordForm coord_form)
    : points_(points), coord_form_(coord_form)
{}

bool ReferencePointsInput::readXML(xmlNodePtr cur)
{
  using modernstrutil::split;

  bool succeeded = true;
  //read in the point contents
      std::vector<std::string_view> lines = split(XMLNodeString{cur}, "\n");
      for (int i = 0; i < lines.size(); i++)
      {
        succeeded                            = true;
        std::vector<std::string_view> tokens = split(lines[i], " ");
        if (tokens.size() != OHMMS_DIM + 1)
        {
          app_log() << "  reference point must have 4 entries, given " << tokens.size() << ": " << lines[i]
                    << std::endl;
          succeeded = false;
          break;
        }
        else
        {
          Point rp;
          for (int d = 0; d < OHMMS_DIM; d++)
          {
            rp[d] = string2Real<Real>(tokens[d + 1]);
          }
	  std::string key{tokens[0]};
          points_[key] = rp;
        }
      }
if (succeeded)
  return true;
else
  return false;
}

std::any ReferencePointsInput::ReferencePointsInputSection::assignAnyEnum(const std::string& name) const
{
  return lookupAnyEnum(name, get<std::string>(name), lookup_input_enum_value);
}
} // namespace qmcplusplus
