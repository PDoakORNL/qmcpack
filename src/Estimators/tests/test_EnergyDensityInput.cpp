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


#include "catch.hpp"

#include "EnergyDensityInput.h"
#include "EstimatorTesting.h"
#include "ValidEnergyDensityInput.h"
#include <iostream>

namespace qmcplusplus
{

TEST_CASE("EnergyDensityInput::parseXML::valid", "[estimators]")
{
  using input = qmcplusplus::testing::ValidEnergyDensityInput;
  for (auto input_xml : input::xml)
    {
      Libxml2Document doc;
    bool okay       = doc.parseFromString(input_xml);
    xmlNodePtr node = doc.getRoot();
    EnergyDensityInput edi(node);
    }

}

TEST_CASE("EnergyDensityInput::parseXML::invalid", "[estimators]")
{
  using input = qmcplusplus::testing::InvalidEnergyDensityInput;
  for (auto input_xml : input::xml)
    {
      Libxml2Document doc;
    bool okay       = doc.parseFromString(input_xml);
    xmlNodePtr node = doc.getRoot();
    auto constructBadEnergyDensityInput = [](xmlNodePtr cur) { EnergyDensityInput edi(cur); };
    CHECK_THROWS_AS(constructBadEnergyDensityInput(node), UniformCommunicateError);    
    }
}

TEST_CASE("EnergyDensityInput::parseXML::axes", "[estimators]")
{
  using input = qmcplusplus::testing::ValidEnergyDensityInput;
  auto& input_xml = input::xml[input::valid::ION];
  Libxml2Document doc;
  bool okay       = doc.parseFromString(input_xml);
  xmlNodePtr node = doc.getRoot();
  EnergyDensityInput edi(node);
  auto sgis = edi.get_space_grid_inputs();
  auto rpts = edi.get_ref_points_input();
  CHECK(sgis.size() == 2);
}
 
}
