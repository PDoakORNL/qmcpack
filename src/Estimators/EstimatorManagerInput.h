//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File refactored from: EstimatorManagerNew.h
//////////////////////////////////////////////////////////////////////////////////////

/** \file
 *  declares the list  supported estimator input types and declares the input type for 
 *  EstimatorManagerNew.
 */
#ifndef QMCPLUSPLUS_ESIMATORMANAGERINPUT_H
#define QMCPLUSPLUS_ESIMATORMANAGERINPUT_H


#include "type_traits/template_types.hpp"
#include <functional>
#include <vector>
#include <variant>
#include <libxml/tree.h>
#include "io/InputNode.hpp"

namespace qmcplusplus
{

class SpinDensityInput;
class MomentumDistributionInput;
class OneBodyDensityMatricesInput;

using EstimatorInput = std::variant<RefW<SpinDensityInput>, RefW<MomentumDistributionInput>, RefW<OneBodyDensityMatricesInput>>;

/** These are the estimator input types EstimatorManagerInput delegates to.
 *  We of course know all the estimator types at compile time and it is useful to have type safety for their usage.
 */
using EstimatorInputs = std::vector<EstimatorInput>;

namespace testing
{
/** class to allow private method testing in unit tests
 */
class EstimatorManagerInputTests;
}

/** Input type for EstimatorManagerNew
 *  Parses Estimators level of input and and delegates child estimator nodes
 *  to appropriate estimator input class
 *  Will later provide access to estimator input instances for esimator construction.
 */
class EstimatorManagerInput
{
public:
  EstimatorManagerInput() = default;
  EstimatorManagerInput(EstimatorManagerInput&& emi) = default;
  EstimatorManagerInput(xmlNodePtr cur);
  EstimatorInputs& get_estimator_inputs() { return estimator_inputs; };
private:
  /** read <Estimators> node
   */
  void readXML(xmlNodePtr cur);

  /// this is a vector of variants for typesafe access to the estimator inputs
  EstimatorInputs estimator_inputs;
  /// this is a type erased ownership vector for the estimator inputs
  UPtrVector<InputNode> estimator_input_storage;

  template<typename T, typename... Args>
  void appendEstimatorInput(Args&&... args)
  {
    estimator_input_storage.push_back(std::make_unique<T>(std::forward<T>(args)...));
    estimator_inputs.push_back(static_cast<T&>(*estimator_input_storage.back()));
  }

  friend class testing::EstimatorManagerInputTests;
};

} // namespace qmcplusplus

#endif
