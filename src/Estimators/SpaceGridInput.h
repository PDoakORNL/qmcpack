//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// Some code refactored from: SpaceGrid.h & SpaceGrid.cpp
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_SPACEGRID_INPUT_H
#define QMCPLUSPLUS_SPACEGRID_INPUT_H

#include <Configuration.h>
#include "InputSection.h"

namespace qmcplusplus
{
namespace testing
{
template<typename T>
class EnergyDensityTests;
}

class SpaceGrid;

class SpaceGridInput
{
public:
  using Consumer = SpaceGrid;
  using Real     = double;

  enum class CoordForm
  {
    CARTESIAN = 0,
    CYLINDRICAL,
    SPHERICAL,
    VORONOI
  };

  inline static const std::unordered_map<std::string, std::any>
      lookup_input_enum_value{{"coord-cartesian", CoordForm::CARTESIAN},
                              {"coord-cyclindrical", CoordForm::CYLINDRICAL},
                              {"coord-spherical", CoordForm::SPHERICAL},
                              {"coord-voronoi", CoordForm::VORONOI}};


  // legal labels for each coordinate type.  These are actually effectively enums
  inline static const std::array<const std::string, 3> ax_cartesian{"x", "y", "z"};
  inline static const std::array<const std::string, 3> ax_cylindrical{"r", "phi", "z"};
  inline static const std::array<const std::string, 3> ax_spherical{"r", "phi", "theta"};

  class SpaceGridAxisInput
  {
    class SpaceGridAxisInputSection : public InputSection
    {
    public:
      SpaceGridAxisInputSection()
      {
        section_name = "axis";
        attributes   = {"label", "grid", "p1", "p2", "scale"};
        strings      = {"label", "p1", "p2"};
	multi_reals = {"grid"},
        reals        = {"scale"};
        required     = {"label","p1"};
      }
    };

  public:
    SpaceGridAxisInput(xmlNodePtr cur) { input_section_.readXML(cur); }
    static std::any makeAxis(xmlNodePtr cur, std::string& value_label)
    {
      SpaceGridAxisInput space_grid_axis{cur};
      value_label= "axis";
      return space_grid_axis;
    }
    const SpaceGridAxisInputSection& get_input() { return input_section_; }
  private:
    SpaceGridAxisInputSection input_section_;
  };

  class SpaceGridOriginInput
  {
    class SpaceGridOriginInputSection : public InputSection
    {
    public:
      SpaceGridOriginInputSection()
      {
        section_name = "origin";
        attributes   = {"p1", "p2", "fraction"};
	required = {"p1"};
        strings      = {"p1", "p2"};
        reals        = {"fraction"};
      }
    };

  public:
    SpaceGridOriginInput(xmlNodePtr cur) { input_section_.readXML(cur); }

    static std::any makeOrigin(xmlNodePtr cur, std::string& value_label)
    {
      SpaceGridOriginInput space_grid_origin{cur};
      value_label = "origin";
      return space_grid_origin;
    }

  private:
    SpaceGridOriginInputSection input_section_;
  };

  class SpaceGridInputSection : public InputSection
  {
  public:
    SpaceGridInputSection()
    {
      section_name = "SpaceGrid";
      attributes   = {"coord"};
      enums        = {"coord"};
      delegates    = {"origin", "axis"};
      multiple     = {"axis"};
      registerDelegate("origin", SpaceGridOriginInput::makeOrigin);
      registerDelegate("axis", SpaceGridAxisInput::makeAxis);
    }
    std::any assignAnyEnum(const std::string& name) const override;
    void checkParticularValidity() override;
  };

  SpaceGridInput(xmlNodePtr cur);

  void checkAxes(std::vector<std::any>& axes);
  
  CoordForm get_coord_form() const { return coord_form_; }
  const std::array<std::string, OHMMS_DIM>& get_axis_p1() const { return axis_p1_; }
  const std::array<Real, OHMMS_DIM>& get_axis_scale() const { return axis_scale_; }
  const std::array<std::string, OHMMS_DIM>& get_axis_label() const { return axis_label_; }
  const std::array<std::vector<double>, OHMMS_DIM>& get_axis_grid() const { return axis_grid_; }
  const std::string& get_origin_p1() const { return origin_p1_; }
  const std::string& get_origin_p2() const { return origin_p2_; }
  Real get_origin_fraction() const { return origin_fraction_; }
private:
  SpaceGridInputSection input_section_;
  CoordForm coord_form_;
  std::string origin_p1_{""};
  std::string origin_p2_{"zero"};
  Real origin_fraction_ = 0.0;
  std::array<std::string, OHMMS_DIM> axis_p1_;
  std::array<Real, OHMMS_DIM> axis_scale_;
  std::array<std::string, OHMMS_DIM> axis_label_;
  std::array<std::vector<double>, OHMMS_DIM> axis_grid_;
};

std::any makeSpaceGridInput(xmlNodePtr, std::string& value_label);

} // namespace qmcplusplus
#endif
