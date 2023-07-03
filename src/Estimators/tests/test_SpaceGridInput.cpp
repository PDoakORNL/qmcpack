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
  using input = qmcplusplus::testing::ValidSpaceGridInput;
  for (auto input_xml : input::xml)
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
  using input = qmcplusplus::testing::InvalidSpaceGridInput;
  for (auto input_xml : input::xml)
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
  using input = qmcplusplus::testing::ValidSpaceGridInput;
  auto& input_xml = input::xml[input::WITH_STEP];
  Libxml2Document doc;
  bool okay = doc.parseFromString(input_xml);
  REQUIRE(okay);
  xmlNodePtr node = doc.getRoot();

  SpaceGridInput sgi(node);

  auto& axis_p1 = sgi.get_axis_p1s();
  CHECK(checkVec(axis_p1, {"r1","r2","r3"}));
  auto& axis_scale = sgi.get_axis_scales();
  CHECK(checkVec(axis_scale, {6.9,6.9,6.9}));
  auto& axis_label = sgi.get_axis_labels();
  CHECK(checkVec(axis_label, {"r","phi","theta"}));
  auto& axis_grid = sgi.get_axis_grids();
  //CHECK(checkVec(axis_grid, {{{0.0, 1.0},{0, 1},{0, 1}}}));
}
  
} // namespace qmcplusplus
