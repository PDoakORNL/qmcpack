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


#include "catch.hpp"

#include "SpaceGridInput.h"
#include "TestingSpaceGridInputs.h"
#include "OhmmsData/Libxml2Doc.h"

namespace qmcplusplus
{

TEST_CASE("SpaceGridInputs::parseXML::valid", "[estimators]")
{
  for (auto input_xml : testing::valid_space_grid_input_sections)
  {
    Libxml2Document doc;
    bool okay       = doc.parseFromString(input_xml);
    xmlNodePtr node = doc.getRoot();

    // Will throw if input is invalid.
    SpaceGridInput spi(node);
  }
}

TEST_CASE("SpaceGridInputs::parseXML::invalid", "[estimators]")
{
  for (auto input_xml : testing::invalid_space_grid_input_sections)
  {
    Libxml2Document doc;
    bool okay = doc.parseFromString(input_xml);
    REQUIRE(okay);
    xmlNodePtr node = doc.getRoot();

    auto constructBadSpaceGrid = [](xmlNodePtr cur) { SpaceGridInput spi(cur); };
    CHECK_THROWS_AS(constructBadSpaceGrid(node), UniformCommunicateError);
  }
}

template<typename T>
bool checkVec(T& vec, T expected_vec) {
  return vec == expected_vec;
}
  
TEST_CASE("SpaceGridInputs::parseXML::axes", "[estimators]")
{
  auto& input_xml = testing::valid_space_grid_input_sections[testing::valid_space_grid_input];
  Libxml2Document doc;
  bool okay = doc.parseFromString(input_xml);
  REQUIRE(okay);
  xmlNodePtr node = doc.getRoot();

  SpaceGridInput sgi(node);

  auto& axis_p1 = sgi.get_axis_p1();
  CHECK(checkVec(axis_p1, {"r1","r2","r3"}));
  auto& axis_scale = sgi.get_axis_scale();
  CHECK(checkVec(axis_scale, {6.9,6.9,6.9}));
  auto& axis_label = sgi.get_axis_label();
  CHECK(checkVec(axis_label, {"r","phi","theta"}));
  auto& axis_grid = sgi.get_axis_grid();
  CHECK(checkVec(axis_grid, {{{0.0, 1.0},{0, 1},{0, 1}}}));
}
  
} // namespace qmcplusplus
