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

namespace qmcplusplus
{
namespace testing
{
template<typename T>
class EnergyDensityTests;
}

class SpaceGridInput
{
public:
  using Consumer = SpaceGrid;
  using Real     = double;

  enum class SpaceGridCoord
  {
    CARTESIAN = 0,
    CYLINDRICAL,
    SPHERICAL,
    VORONOI
  };

  inline static const std::unordered_map<std::string, std::any>
      lookup_input_enum_value{{"spacegridcoord-cartesian", SpaceGridCoord::CARTESIAN},
                              {"spacegridcoord-cyclindrical", SpaceGridCoord::CYLINDRICAL},
                              {"spacegridcoord-spherical", SpaceGridCoord::SPHERICAL},
                              {"spacegridcoord-voronoi", SpaceGridCoord::VORONOI}};

  class SpaceGridAxisInput
  {
    class SpaceGridAxisInputSection : public InputSection
    {
    public:
      SpaceGridAxisInputSection()
      {
        section_name = "axis";
        attributes   = {"label", "grid", "p1", "p2", "scale"};
        strings      = {"label", "grid", "p1", "p2"};
        reals        = {"scale"};
      }
    };

  public:
    SpaceGridAxisInput(xmlNodePtr cur) { input_section_.readXML(cur); }
    static std::any makeAxis(xmlNodePtr cur, std::string& value_label)
    {
      SpaceGridAxisInput space_grid_axis{cur};
      value_label = "spacegridorigin";
      return space_grid_axis;
    }
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
        strings      = {"p1", "p2"};
        reals        = {"fraction"};
      }
    };

  public:
    SpaceGridOriginInput(xmlNodePtr cur) { input_section_.readXML(cur); }

    static std::any makeOrigin(xmlNodePtr cur, std::string& value_label)
    {
      SpaceGridOriginInput space_grid_origin{cur};
      value_label = "spacegridorigin";
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
      registerDelegate("origin", SpaceGridOriginInput::makeOrigin);
      registerDelegate("axis", SpaceGridAxisInput::makeAxis);
    }
    std::any assignAnyEnum(const std::string& name) const override;
  };

private:
  SpaceGridInputSection input_section_;
  SpaceGridCoord space_grid_coord;
  std::string origin_p1;
  std::string origin_p2;
  Real origin_fraction;
  std::array<std::string, OHMMS_DIM> axis_p1;
  std::array<Real, OHMMS_DIM> axis_scale;
  std::array<std::string, OHMMS_DIM> axis_label;
  std::array<std::string, OHMMS_DIM> axis_grid;
};

std::any makeSpaceGridInput(xmlNodePtr, std::string& value_label);

} // namespace qmcplusplus
#endif
