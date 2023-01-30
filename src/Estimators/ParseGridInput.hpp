//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2023 QMCPACK developers.
//
// File developed by: Peter  W. Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Peter  W. Doak, doakpw@ornl.gov, Oak Ridge National Lab
//////////////////////////////////////////////////////////////////////////////////////

/** \file
 *  The AxisGrid data structure and the ParseGridInput factor parsing in a manner usable in acustom handler
 *  method of any InputSection subtype. Right now it is only required by
 *  SpaceGridInput::SpaceGridAxisInputSection::setFromStreamCustom.
 *  This takes the one dimensional grid input specified for EnergyDensityEstimator (described in the for
 *  the legacy energy density estimator). To the suprising complex input variables required by SpaceGrid.
 *  AxisGrid should be subjected to analysis at a later date for which of these data elements are actually
 *  required to successfully construct a NESpaceGrid instance.
 */

#ifndef QMCPLUSPLUS_PARSE_GRID_INPUT_HPP
#define QMCPLUSPLUS_PARSE_GRID_INPUT_HPP
#include <vector>
#include <string>
#include <sstream>
#include <iostream>

namespace qmcplusplus
{

template<typename REAL>
struct AxisGrid
{
  // I think the parse should just return something like this.
  // std::vector<REAL> points;
  // std::vector<REAL> p_to_p_du;
  // however the parse is hopelessly tangled up with a great deal of calculation.
  // which produces this
  std::vector<int> ndom_int;
  std::vector<int> ndu_int;
  std::vector<REAL> du_int;
  REAL umin;
  REAL umax;
  REAL odu;
  std::vector<int> gmap;
  std::vector<int> ndu_per_interval;
  int dimensions;

  /** equality operator, all values must be equal.
   *  In C++20 this will just be the defaulted operator==.
   */
  bool operator==(const AxisGrid& ag) const {
    return (ag.ndom_int == ndom_int && ag.ndu_int == ndu_int && ag.du_int == du_int
	    && ag.umin == umin && ag.umax == umax && ag.odu == odu && ag.gmap == gmap &&
	    ag.ndu_per_interval == ndu_per_interval && ag.dimensions == dimensions);
  }
};

/** Parses the one dimensional grid specification format originated for EnergyDensity
 *  This has be refactored from QMCHamiltonian/SpaceGrid.cpp but it is ridiculously hard to read.
 * My evidence for correctness is that it can pass the unit test I wrote.  By inspection I would have guessed it would be broken.
 */
template<typename REAL>
AxisGrid<REAL> parseGridInput(std::istringstream& grid_input_stream);

// explicit instantiation declaration.
extern template struct AxisGrid<float>;
extern template struct AxisGrid<double>;
extern template AxisGrid<double> parseGridInput<double>(std::istringstream& grid_input_stream);
extern template AxisGrid<float> parseGridInput<float>(std::istringstream& grid_input_stream);

} // namespace qmcplusplus

#endif
