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

#include "ReferencePointsInput.h"

#include <stdio.h>
#include <sstream>

#include "ValidReferencePointsInput.h"
#include "OhmmsData/Libxml2Doc.h"

namespace qmcplusplus
{

TEST_CASE("ReferencePointsInput::parseXML", "[estimators]")
{
  for (auto input_xml : testing::valid_reference_points_input_sections)
  {
    Libxml2Document doc;
    bool okay = doc.parseFromString(input_xml);
    REQUIRE(okay);
    xmlNodePtr node = doc.getRoot();

    ReferencePointsInput rpi(node);
  }

  for (auto input_xml : testing::invalid_reference_points_input_sections)
  {
    Libxml2Document doc;
    bool okay = doc.parseFromString(input_xml);
    REQUIRE(okay);
    xmlNodePtr node = doc.getRoot();

    CHECK_THROWS(ReferencePointsInput(node));
  }

}

TEST_CASE("ReferencePointsInput::Construction", "[estimators]")
{
  ReferencePointsInput::Points points{{"r1",{0.0, 0.0, 1.0}}, {"r2",{0.0, 1.0, 0.0}}, {"r3",{1.0, 0.0, 0.0}}};
  ReferencePointsInput{points, ReferencePointsInput::CoordForm::CELL};
}
  
}


