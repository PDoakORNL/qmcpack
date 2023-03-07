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

#include <set>
#include <unordered_map>

#include <Configuration.h>
#include "InputSection.h"
#include "EstimatorInput.h"
#include "ParseGridInput.hpp"

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
    SPHERICAL
  };

  inline static const std::unordered_map<std::string, std::any> lookup_input_enum_value{{"coord-cartesian",
                                                                                         CoordForm::CARTESIAN},
                                                                                        {"coord-cyclindrical",
                                                                                         CoordForm::CYLINDRICAL},
                                                                                        {"coord-spherical",
                                                                                         CoordForm::SPHERICAL}};

  using label_set = std::set<std::string_view>;
  // legal labels for each coordinate type.  These are actually effectively enums
  inline static const label_set ax_cartesian{"x", "y", "z"};
  inline static const label_set ax_cylindrical{"r", "phi", "z"};
  inline static const label_set ax_spherical{"r", "phi", "theta"};

  inline static const std::unordered_map<CoordForm, label_set> axes_label_sets{{CoordForm::CARTESIAN, ax_cartesian},
                                                                               {CoordForm::CYLINDRICAL, ax_cylindrical},
                                                                               {CoordForm::SPHERICAL, ax_spherical}};

  class SpaceGridAxisInput
  {
    class SpaceGridAxisInputSection : public InputSection
    {
    public:
      SpaceGridAxisInputSection()
      {
        section_name = "axis";
        attributes   = {"label", "p1", "p2", "scale"};
        strings      = {"label", "p1", "p2"};
        custom_attributes = {"grid"},
	reals = {"scale"};
        required = {"label", "p1"};
      }
      void setFromStreamCustom(const std::string& ename, const std::string& name, std::istringstream& svalue) override;
    };

  public:
    SpaceGridAxisInput(xmlNodePtr cur);

    static std::any makeAxis(xmlNodePtr cur, std::string& value_label)
    {
      SpaceGridAxisInput space_grid_axis{cur};

      value_label = "axis";
      return space_grid_axis;
    }

    const SpaceGridAxisInputSection& get_input() { return input_section_; }

    std::string get_label() const { return label_; }
    Real get_scale() const { return scale_; }
    std::string get_p1() const { return p1_; }
    std::string get_p2() const { return p2_; }
    AxisGrid<Real> get_grid() const { return grid_; }

  private:
    SpaceGridAxisInputSection input_section_;
    std::string label_ = "";
    Real scale_        = 1.0;
    std::string p1_    = "";
    std::string p2_    = "zero";
    AxisGrid<Real> grid_;
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
        required     = {"p1"};
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
  const std::array<std::string, OHMMS_DIM>& get_axis_p1s() const { return axis_p1s_; }
  const std::array<std::string, OHMMS_DIM>& get_axis_p2s() const { return axis_p2s_; }

  const std::array<Real, OHMMS_DIM>& get_axis_scales() const { return axis_scales_; }
  const std::array<std::string, OHMMS_DIM>& get_axis_labels() const { return axis_labels_; }
  const std::array<AxisGrid<Real>, OHMMS_DIM>& get_axis_grids() const { return axis_grids_; }
  const std::string& get_origin_p1() const { return origin_p1_; }
  const std::string& get_origin_p2() const { return origin_p2_; }
  Real get_origin_fraction() const { return origin_fraction_; }
  /** axes_label_set accessor, avoids a bunch of switch statements
   *  at must be used because std::unordered_map::operator[] can't return a const reference
   */
  const label_set& get_axes_label_set() const { return axes_label_sets.at(coord_form_); }

private:
  SpaceGridInputSection input_section_;
  CoordForm coord_form_;
  std::string origin_p1_{""};
  std::string origin_p2_{"zero"};
  Real origin_fraction_ = 0.0;
  std::array<std::string, OHMMS_DIM> axis_labels_;
  std::array<std::string, OHMMS_DIM> axis_p1s_;
  std::array<std::string, OHMMS_DIM> axis_p2s_;
  std::array<Real, OHMMS_DIM> axis_scales_;
  std::array<AxisGrid<Real>, OHMMS_DIM> axis_grids_;
};

std::any makeSpaceGridInput(xmlNodePtr, std::string& value_label);

} // namespace qmcplusplus
#endif
