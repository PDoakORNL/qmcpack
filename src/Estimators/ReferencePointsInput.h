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

#ifndef QMCPLUSPLUS_REFERENCE_POINTS_INPUT_H
#define QMCPLUSPLUS_REFERENCE_POINTS_INPUT_H

#include <any>
#include <type_traits>

#include "Configuration.h"
#include "InputSection.h"
#include "OhmmsData/Libxml2Doc.h"

namespace qmcplusplus
{

class ReferencePoints;

class ReferencePointsInput
{
public:
  using Real     = QMCTraits::RealType;
  using Consumer = ReferencePoints;
  using Point    = TinyVector<Real, OHMMS_DIM>;
  using Points   = std::map<std::string, Point>;

  enum class CoordForm
  {
    CELL = 0,
    CARTESIAN
  };

  inline static const std::unordered_map<std::string, std::any> lookup_input_enum_value{{"coord-cell", CoordForm::CELL},
                                                                                        {"coord-cartesian",
                                                                                         CoordForm::CARTESIAN}};

  class ReferencePointsInputSection : public InputSection
  {
  public:
    ReferencePointsInputSection()
    {
      section_name = "ReferencePoints";
      attributes   = {"coord"};
      enums        = {"coord"};
      required     = {"coord"};
    }

    std::any assignAnyEnum(const std::string& name) const override;
  };

  ReferencePointsInput(xmlNodePtr cur);
  ReferencePointsInput(ReferencePointsInput& other)  = default;
  ReferencePointsInput(ReferencePointsInput&& other) = default;
  ReferencePointsInput(const Points& points, const CoordForm coord_form = CoordForm::CARTESIAN);

  CoordForm get_coord_form() const { return coord_form_; }
  const Points& get_points() const { return points_; }

private:
  bool readXML(xmlNodePtr);

  ReferencePointsInputSection input_section_;
  Points points_;
  CoordForm coord_form_;
};

} // namespace qmcplusplus
#endif
