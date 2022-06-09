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

#ifndef QMCPLUSPLUS_VALID_REFERECEPOINTS_INPUT_H
#define QMCPLUSPLUS_VALID_REFERECEPOINTS_INPUT_H

#include <array>

namespace qmcplusplus
{
namespace testing
{


constexpr std::array<std::string_view, 2> valid_reference_points_input_sections{
    R"XML(
  <reference_points coord="cell">
    r1 1 0 0
    r2 0 1 0
    r3 0 0 1
  </reference_points>
)XML",
    R"XML(
  <reference_points coord="cartesian">
    r1 1 0 0
    r2 0 1 0
    r3 0 0 1
  </reference_points>
)XML"};

constexpr std::array<std::string_view, 1> invalid_reference_points_input_sections{
    R"XML(
  <reference_points coord="cartesian">
    r1 1 0 0
    r2 0 1
    r3 0 0 1
  </reference_points>
)XML"};
/* )XML", */
/*     R"XML( */
/*   <reference_points coord="cartesian"> */
/*     r1 1 0 0 */
/*     r2 0 1 0 */
/*   </reference_points> */
/* )XML"}; */

constexpr int valid_referencepointsinput_cell        = 0;
constexpr int valid_referencepointsinput_cartesian   = 1;
constexpr int invalid_referencepointsinput_point_param = 0;
/* constexpr int invalid_referencepointsinput_point_num = 1; */

} // namespace testing
} // namespace qmcplusplus

#endif
