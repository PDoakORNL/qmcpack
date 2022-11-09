//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// Some code refactored from: DensityMatrices1b.h
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_ENERGY_DENSITY_INPUT_H
#define QMCPLUSPLUS_ENERGY_DENSITY_INPUT_H

#include "ReferencePointsInput.h"
#include "SpaceGrid.h"

namespace qmcplusplus
{
namespace testing
{
template<typename T>
class EnergyDensityTests;
}

class EnergyDensityEstimator;

/** EnergyDensity has two other XML input reading objects that it delegates to.
    I don't think these need to be handled with a variant vector because there
    is no expectation it be expanded and they in general occur together and not as options.
 */
// using EnergyDensityDelegate =
//     std::variant<std::monostate, ReferencePointsInput, SpaceGridInput >;
// using EnergyDensityDelegates = std::vector<EnergyDensityDelegate>;

/** Native representation for DensityMatrices1B Estimator's inputs
 */
class EnergyDensityInput
{
public:
  using Consumer = EnergyDensityEstimator;

  class EnergyDensityInputSection : public InputSection
  {
    public:
    EnergyDensityInputSection()
    {
      section_name = "EnergyDensity";
      attributes = {"name", "dynamic", "static", "ion_points", "type"};
      parameters = {"reference_points","spacegrid"};
      strings = {"name", "type", "dynamic", "static"};
      booleans = {"ion_points"};
      delegates = {"reference_points","spacegrid"};
      required = {"type", "spacegrid"};
      registerDelegate("reference_points", makeReferencePointsInput);
    }
    /** Here the delegate input object is registered */
    EnergyDensityInputSection();
    EnergyDensityInputSection(const EnergyDensityInputSection&) = default;
    bool setFromStreamCustom(const std::string& name, istringstream sstream) override;
  };

  using Real = QMCTraits::RealType;

  EnergyDensityInput(const EnergyDensityInput&) = default;
  EnergyDensityInput(xmlNodePtr cur);

  const std::string& get_name() const { return name_; }
  const std::string& get_dynamic() const { return dynamic_; }
  const std::string& get_static() const { return static_; }
  const ReferencePointsInput& get_ref_points_input() const { return ref_points_input_; }
  const std::vector<SpaceGridsInput>& get_space_grid_inputs() const { return space_grid_inputs_; }
private:
  std::string name_;
  std::string dynamic_;
  std::string static_;
  ReferencePointsInput ref_points_input_;
  std::vector<SpaceGridsInput> space_grid_inputs_;
}

}

#endif
