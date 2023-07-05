//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2023 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_VALID_SPACEGRID_INPUT_H
#define QMCPLUSPLUS_VALID_SPACEGRID_INPUT_H

#include <array>
#include <string_view>

namespace qmcplusplus
{
namespace testing
{

struct ValidSpaceGridInput
{
  static constexpr std::array<std::string_view, 2> xml{
      R"XML(
           <spacegrid coord="cartesian">
             <axis p1="a1" scale=".5" label="x" grid="-1 (.1) 1"/>
             <axis p1="a2" scale=".5" label="y" grid="-1 (.1) 1"/>
             <axis p1="a3" scale=".5" label="z" grid="-1 (.1) 1"/>
           </spacegrid>
      )XML",
      R"XML(
           <spacegrid coord="cartesian">
             <origin p1="zero"/>
             <axis p1="a1" scale=".5" label="x" grid="-1 (.1) 1"/>
             <axis p1="a2" scale=".5" label="y" grid="-1 (.1) 1"/>
             <axis p1="a3" scale=".5" label="z" grid="-1 (.1) 1"/>
           </spacegrid>
      )XML"};
  enum valid
  {
    DEFAULT = 0,
    ORIGIN
  };
};

struct InvalidSpaceGridInput
{
  static constexpr std::array<std::string_view, 1> xml{
      R"XML(
           <spacegrid coord="cartesian">
             <axis p1="a1" scale=".5" label="x" grid="-1 (.1) 1"/>
             <axis p1="a2" scale=".5" label="y" grid="-1 (.1) 1"/>
             <axis p1="a3" scale=".5" label="z" grid="-1 (.1) 1"/>
           </spacegrid>
      )XML"};
};
} // namespace testing
} // namespace qmcplusplus

#endif
