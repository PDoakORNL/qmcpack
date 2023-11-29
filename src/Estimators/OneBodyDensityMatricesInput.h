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

#ifndef QMCPLUSPLUS_ONE_BODY_DENSITY_MATRICES_INPUT_H
#define QMCPLUSPLUS_ONE_BODY_DENSITY_MATRICES_INPUT_H

#include "Configuration.h"
#include "InputSection.h"

namespace qmcplusplus
{
namespace testing
{
template<typename T>
class OneBodyDensityMatricesTests;
}

class OneBodyDensityMatrices;

/** Native representation for DensityMatrices1B Estimator's inputs
 */
class OneBodyDensityMatricesInput
{
public:
  using Consumer = OneBodyDensityMatrices;

  enum class Integrator
  {
    UNIFORM_GRID,
    UNIFORM,
    DENSITY
  };

  enum class Evaluator
  {
    LOOP,
    MATRIX
  };

  /** mapping for enumerated options of OneBodyDensityMatrices
   *  This data object is the basis of input enum string values
   *  translation to native C++ scoped enums. This boiler plate is
   *  generated by utils/code_tools/qmcpack-elisp-generators.el
   *  qmcp-add-enum-string-map
   *
   *  This plus the virtual assignAnyEnum method are needed by InputSection to
   *  validate and assign enum values from input.
   *
   *  In testing code we assume this map is bidirectional.
   */
  inline static const std::unordered_map<std::string, std::any>
      lookup_input_enum_value{{"integrator-uniform_grid", Integrator::UNIFORM_GRID},
                              {"integrator-uniform", Integrator::UNIFORM},
                              {"integrator-density", Integrator::DENSITY},
                              {"evaluator-loop", Evaluator::LOOP},
                              {"evaluator-matrix", Evaluator::MATRIX}};

  class OneBodyDensityMatricesInputSection : public InputSection
  {
  public:
    /** parse time definition of input parameters */
    OneBodyDensityMatricesInputSection()
    {
      // clang-format off
      section_name  = "OneBodyDensityMatrices";
      attributes    = {"name", "type"};
      parameters    = {"basis", "energy_matrix", "integrator", "evaluator", "scale",
                       "corner", "center", "points", "samples", "warmup", "timestep",
                       "use_drift", "check_overlap", "check_derivatives", "acceptance_ratio", "rstats",
                       "normalized", "volumed_normed"};
      bools         = {"energy_matrix", "use_drift", "normalized", "volume_normed",
                       "check_overlap", "check_derivatives", "rstats", "acceptance_ratio"};
      enums         = {"integrator", "evaluator"};
      strings       = {"name", "type"};
      multi_strings = {"basis"};
      integers      = {"points", "samples"};
      reals         = {"scale", "timestep"};
      positions     = {"center", "corner"};
      required      = {"name", "basis"};
      // I'd much rather see the default defined in simple native c++ as below
      // clang-format on
    }
    OneBodyDensityMatricesInputSection(const OneBodyDensityMatricesInputSection&) = default;
    /** do parse time checks of input */
    void checkParticularValidity() override;
    std::any assignAnyEnum(const std::string& name) const override;
  };

  using Real     = QMCTraits::FullPrecRealType;
  using Position = TinyVector<Real, OHMMS_DIM>;

  /** default copy constructor
   *  This is required due to OBDMI being part of a variant used as a vector element.
   */
  OneBodyDensityMatricesInput(const OneBodyDensityMatricesInput&) = default;
  OneBodyDensityMatricesInput(xmlNodePtr cur);

private:
  OneBodyDensityMatricesInputSection input_section_;

  // Default parameters for OneBodyDensityMatrices
  bool energy_matrix_          = false;
  bool use_drift_              = false;
  bool normalized_             = true;
  bool volume_normalized_      = true;
  bool check_overlap_          = false;
  bool check_derivatives_      = false;
  bool rstats_                 = false;
  bool write_acceptance_ratio_ = false;
  /// This flag is derived from input so if you construct an OBDMI directly with center it must be set.
  bool center_defined_   = false;
  bool corner_defined_   = false;
  Integrator integrator_ = Integrator::UNIFORM_GRID;
  Evaluator evaluator_   = Evaluator::LOOP;
  Real scale_            = 1.0;
  /// center_ does not have a default. The estimator sets if from input Lattice if it isn't set
  Position center_;
  Position corner_;
  Real timestep_      = 0.5;
  int points_         = 10;
  int samples_        = 10;
  int warmup_samples_ = 30;
  std::vector<std::string> basis_sets_;

public:
  bool get_energy_matrix() const { return energy_matrix_; }
  bool get_use_drift() const { return use_drift_; }
  bool get_normalized() const { return normalized_; }
  bool get_volume_normalized() const { return volume_normalized_; }
  bool get_check_overlap() const { return check_overlap_; }
  bool get_check_derivatives() const { return check_derivatives_; }
  bool get_rstats() const { return rstats_; }
  bool get_write_acceptance_ratio() const { return write_acceptance_ratio_; }
  Integrator get_integrator() const { return integrator_; }
  Evaluator get_evaluator() const { return evaluator_; }
  Real get_scale() const { return scale_; }
  Position get_center() const { return center_; }
  Position get_corner() const { return corner_; }
  bool get_center_defined() const { return center_defined_; }
  bool get_corner_defined() const { return corner_defined_; }
  Real get_timestep() const { return timestep_; }
  int get_points() const { return points_; }
  int get_samples() const { return samples_; }
  int get_warmup_samples() const { return warmup_samples_; }
  const std::vector<std::string>& get_basis_sets() const { return basis_sets_; }
  template<typename T>
  friend class testing::OneBodyDensityMatricesTests;
};

} // namespace qmcplusplus

#endif
