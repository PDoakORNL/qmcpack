//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// Some code refactored from: ReferencePoints.h & ReferencePoints.cpp
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_REFERNCE_POINTS_INPUT_H
#define QMCPLUSPLUS_REFERNCE_POINTS_INPUT_H

#include "NESpaceGrid.h"

namespace qmcplusplus
{
namespace testing
{
template<typename T>
class EnergyDensityTests;

class ReferencePointsInput
{
public:
  using Consumer = SpaceGrid;
  class ReferencePointsInputSection : public InputSection
  {
  public:
    ReferencePointInputSection()
    {
      section_name = "ReferencePoints";
      attributes = {"name", "type", "coord", "min_part", "max_part", "reference", "periodic"};
    }
  };
private:
  Coordinate coordinate;
  
}
};
}

#endif
