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

namespace qmcplusplus
{
namespace testing
{
template<typename T>
class EnergyDensityTests;

class SpaceGridInput
{
public:
  using Consumer = SpaceGrid;

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

  class SpaceGridAxisInputSection : public InputSection
  {
    SpaceGridAxisInputSection()
    {
      section_name = "axis";
      attributes   = {"label", "grid", "p1", "p2", "scale"};
      strings      = {"label", "grid", "p1", "p2"};
      reals        = {"scale"};
    }
  };

  std::any makeAxis(xmlNodePtr cur)
  {
    SpaceGridAxisInputSection space_grid_axis;
    space_grid_axis.readXML(cur);
    return space_grid_axis;
  }

  class SpaceGridOriginInputSection : public InputSection
  {
    SpaceGridOriginInputSection()
    {
      section_name = "origin";
      attributes   = {"p1", "p2", "fraction"};
      strings      = {"p1", "p2"};
      reals        = {"fraction"};
    }
  };

  std::any makeOrigin(xmlNodePtr cur)
  {
    SpaceGridOriginInputSection space_grid_origin;
    space_origin_origin.readXML(cur);
    return space_grid_origin;
  }

  class SpaceGridInputSection : public InputSection
  {
  public:
    SpaceGridInputSection()
    {
      section_name = "SpaceGrid";
      attributes   = {"coord"};
      enums        = {"coord"};
      delegates    = {"origin", "axis"};
      registerDelegate("origin", makeOrigin);
      registerDelegate("axis", makeAxis);
    }
    std::any assignAnyEnum(const std::string& name) const override;
  };

private:
  SpaceGridCoord space_grid_coord;
  std::string origin_p1;
  std::string origin_p2;
  Real origin_fraction;
  std::array<std::string, OHMMSDIM> axis_p1;
  std::array<Real, OHMMSDIM> axis_scale;
  std::array<std::string, OHMMSDIM> axis_label;
  std::array<std::string, OHMMSDIM> axis_grid;
}
}; // namespace testing
} // namespace qmcplusplus
#endif
