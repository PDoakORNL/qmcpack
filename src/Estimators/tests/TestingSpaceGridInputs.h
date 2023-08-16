//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_TESTING_SPACEGRIDINPUTS_H
#define QMCPLUSPLUS_TESTING_SPACEGRIDINPUTS_H

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
  <spacegrid coord="spherical">
    <origin p1="ion1p"/>
    <axis p1="r1" scale="6.9" label="r"     grid="0 (0.1) 1"/>
    <axis p1="r2" scale="6.9" label="phi"   grid="0 (0.1) 1"/>
    <axis p1="r3" scale="6.9" label="theta" grid="0 (0.1) 1"/>
  </spacegrid>
)XML",
      R"XML(
  <spacegrid coord="spherical">
    <origin p1="ion2"/>
    <axis p1="r1" scale="6.9" label="r"     grid="0 1"/>
    <axis p1="r2" scale="6.9" label="phi"   grid="0 1"/>
    <axis p1="r3" scale="6.9" label="theta" grid="0 1"/>
  </spacegrid>
)XML"};

  enum valid
  {
    WITH_STEP = 0,
    WITHOUT_STEP
  };
};

struct InvalidSpaceGridInput
{
  static constexpr std::array<std::string_view, 3> xml{
      R"XML(
  <spacegrid coord="sphericalt">
    <origin p1="ion1p"/>
    <axis p1="r6" scale="6.9" label="r"     grid="0 (0.1) 1"/>
    <axis p1="r2" scale="6.9" label="phi"   grid="0 (0.1) 1"/>
    <axis p1="r3" scale="6.9" label="theta" grid="0 (0.1) 1"/>
  </spacegrid>
)XML",
      R"XML(
  <spacegrid coord="sphericalt">
    <origin p1="ion2"/>
    <axis p1="r1" scale="6.9" label="x"     grid="0 1"/>
    <axis p1="r2" scale="6.9" label="phi"   grid="0 1"/>
    <axis p1="r3" scale="6.9" label="theta" grid="0 1"/>
  </spacegrid>
)XML",
      R"XML(
  <spacegrid coord="spherical">
    <origin p1="ion2"/>
    <axis p1="r2" scale="6.9" label="phi"   grid="0 1"/>
    <axis p1="r3" scale="6.9" label="theta" grid="0 1"/>
  </spacegrid>
)XML"};
};

} // namespace testing
} // namespace qmcplusplus

#endif
