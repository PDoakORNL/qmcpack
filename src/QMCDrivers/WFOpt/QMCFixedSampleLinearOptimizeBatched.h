//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


/** @file QMCFixedSampleLinearOptimizeBatched.h
 * @brief Definition of QMCDriver which performs VMC and optimization.
 */
#ifndef QMCPLUSPLUS_QMCFSLINEAROPTIMIZATION_BATCHED_H
#define QMCPLUSPLUS_QMCFSLINEAROPTIMIZATION_BATCHED_H

#include "QMCDrivers/QMCDriverNew.h"
#include "QMCDrivers/QMCDriverInput.h"
#include "QMCDrivers/VMC/VMCDriverInput.h"
#include "QMCDrivers/VMC/VMCBatched.h"
#include "Optimize/NRCOptimization.h"
#include "Optimize/NRCOptimizationFunctionWrapper.h"
#ifdef HAVE_LMY_ENGINE
#include "formic/utils/matrix.h"
#include "formic/utils/lmyengine/engine.h"
#endif
#include "QMCDrivers/Optimizers/DescentEngine.h"
#include "QMCDrivers/Optimizers/HybridEngine.h"
#include "OutputMatrix.h"
#include "LinearMethod.h"

namespace qmcplusplus
{
/** @ingroup QMCDrivers
 * @brief Implements wave-function optimization
 *
 * Optimization by correlated sampling method with configurations
 * generated from VMC.
 */


///forward declaration of a cost function
class QMCCostFunctionBase;

class GradientTest;


class QMCFixedSampleLinearOptimizeBatched : public QMCDriverNew, LinearMethod
{
public:
  ///Constructor.
  QMCFixedSampleLinearOptimizeBatched(const ProjectData& project_data,
                                      QMCDriverInput&& qmcdriver_input,
                                      const std::optional<EstimatorManagerInput>& global_emi,
                                      VMCDriverInput&& vmcdriver_input,
                                      WalkerConfigurations& wc,
                                      MCPopulation&& population,
                                      SampleStack& samples,
                                      Communicate* comm);

  ///Destructor
  ~QMCFixedSampleLinearOptimizeBatched() override;

  QMCRunType getRunType() override { return QMCRunType::LINEAR_OPTIMIZE; }

  void setWaveFunctionNode(xmlNodePtr cur) { wfNode = cur; }

  ///Run the Optimization algorithm.
  bool run() override;
  ///preprocess xml node
  void process(xmlNodePtr cur) override;
  ///process xml node value (parameters for both VMC and OPT) for the actual optimization
  bool processOptXML(xmlNodePtr cur, const std::string& vmcMove, bool reportH5, bool useGPU);

  RealType costFunc(RealType dl);

  ///common operation to start optimization
  void start();

#ifdef HAVE_LMY_ENGINE
  using ValueType = QMCTraits::ValueType;
  void engine_start(cqmc::engine::LMYEngine<ValueType>* EngineObj,
                    DescentEngine& descentEngineObj,
                    std::string MinMethod);
#endif


  ///common operation to finish optimization, used by the derived classes
  void finish();

  void generateSamples();


private:
  NRCOptimizationFunctionWrapper<QMCFixedSampleLinearOptimizeBatched> objFuncWrapper_;

  inline bool ValidCostFunction(bool valid)
  {
    if (!valid)
      app_log() << " Cost Function is Invalid. If this frequently, try reducing the step size of the line minimization "
                   "or reduce the number of cycles. "
                << std::endl;
    return valid;
  }

  // check if the proposed new cost function value is the best available
  bool is_best_cost(const int ii,
                    const std::vector<RealType>& cv,
                    const std::vector<double>& sh,
                    const RealType ic) const;

  // perform the adaptive three-shift update
  bool adaptive_three_shift_run();

  // perform the single-shift update, no sample regeneration
  bool one_shift_run();

  // perform optimization using a gradient descent algorithm
  bool descent_run();

  // Previous linear optimizers ("quartic" and "rescale")
  bool previous_linear_methods_run();


#ifdef HAVE_LMY_ENGINE
  // use hybrid approach of descent and blocked linear method for optimization
  bool hybrid_run();
#endif

  // Perform test of gradients
  bool test_run();

  std::unique_ptr<GradientTest> testEngineObj;


  void solveShiftsWithoutLMYEngine(const std::vector<double>& shifts_i,
                                   const std::vector<double>& shiffts_s,
                                   std::vector<std::vector<RealType>>& parameterDirections);

#ifdef HAVE_LMY_ENGINE
  formic::VarDeps vdeps;
  cqmc::engine::LMYEngine<ValueType>* EngineObj;
#endif

  //engine for running various gradient descent based algorithms for optimization
  std::unique_ptr<DescentEngine> descentEngineObj;

  //engine for controlling a optimization using a hybrid combination of linear method and descent
  std::unique_ptr<HybridEngine> hybridEngineObj;

  // prepare a vector of shifts to try
  std::vector<double> prepare_shifts(const double central_shift) const;

  // previous update
#ifdef HAVE_LMY_ENGINE
  std::vector<formic::ColVec<double>> previous_update;
#endif

  void print_cost_summary_header();
  void print_cost_summary(const double si,
                          const double ss,
                          const RealType mc,
                          const RealType cv,
                          const int ind,
                          const int bi,
                          const bool gu);

  // ------------------------------------
  // Used by legacy linear method algos

  std::vector<RealType> optdir, optparam;

  ///Number of iterations maximum before generating new configurations.
  int Max_iterations;

  RealType param_tol;
  //-------------------------------------


  ///Number of iterations maximum before generating new configurations.
  int nstabilizers;
  RealType stabilizerScale, bigChange, exp0;
  /// the previous best identity shift
  RealType bestShift_i;
  /// the previous best overlap shift
  RealType bestShift_s;
  /// current shift_i, shift_s input values
  RealType shift_i_input, shift_s_input;
  /// accept history, remember the last 2 iterations, value 00, 01, 10, 11
  std::bitset<2> accept_history;
  /// Shift_s adjustment base
  RealType shift_s_base;

  // ------------------------------------
  // Parameters in this struct are used by one or more of the adaptive LM, descent, or hybrid optimizers

  //String inputs are listed outside the struct
  ///name of the current optimization method, updated by processOptXML before run
  std::string MinMethod;

  //LMY related input
  struct LMYOptions
  {
    /// number of shifts we will try
    int num_shifts = 3;
    /// the maximum relative change in the cost function for the adaptive three-shift scheme
    RealType max_relative_cost_change = 10.0;
    ///max amount a parameter may change relative to current wave function weight
    RealType max_param_change = 0.3;
    /// the tolerance to cost function increases when choosing the best shift in the adaptive shift method
    RealType cost_increase_tol = 0.0;
    /// the shift_i value that the adaptive shift method should aim for
    RealType target_shift_i = -1.0;
    ///whether we are targeting an excited state
    bool targetExcited = false;
    ///whether we are doing block algorithm
    bool block_lm = false;
    ///number of blocks used in block algorithm
    int nblocks = 1;
    ///number of old updates kept
    int nolds = 1;
    ///number of directions kept
    int nkept = 1;
    ///number of samples to do in correlated sampling part
    int nsamp_comp = 0;
    ///the shift to use when targeting an excited state
    RealType omega_shift = 0.0;
    ///whether to do the first part of block lm
    bool block_first = true;
    ///whether to do the second part of block lm
    bool block_second = false;
    ///whether to do the third part of block lm
    bool block_third = false;
    ///whether to filter parameters for the lm
    bool filter_param = false;
    ///whether to filter parameters for the lm
    bool filter_info = false;
    ///threshold for filtering parameters for the lm
    double ratio_threshold = 0.0;
    ///whether to store samples for the lm
    bool store_samples = false;
    ///type of the previous optimization method, updated by processOptXML before run
    OptimizerType previous_optimizer_type = OptimizerType::NONE;
    ///type of the current optimization method, updated by processOptXML before run
    OptimizerType current_optimizer_type = OptimizerType::NONE;
    ///whether to use hybrid method
    bool doHybrid = false;
  };

  /// LMY engine related options
  LMYOptions options_LMY_;

  // ------------------------------------

  // Test parameter gradients
  bool doGradientTest;

  // Output Hamiltonian and overlap matrices
  bool do_output_matrices_csv_;

  // Output Hamiltonian and overlap matrices in HDF format
  bool do_output_matrices_hdf_;

  // Flag to open the files on first pass and print header line
  bool output_matrices_initialized_;

  OutputMatrix output_hamiltonian_;
  OutputMatrix output_overlap_;

  // Freeze variational parameters.  Do not update them during each step.
  bool freeze_parameters_;

  NewTimer& generate_samples_timer_;
  NewTimer& initialize_timer_;
  NewTimer& eigenvalue_timer_;
  NewTimer& involvmat_timer_;
  NewTimer& line_min_timer_;
  NewTimer& cost_function_timer_;

  ///xml node to be dumped
  xmlNodePtr wfNode;

  ParameterSet m_param;

  ///target cost function to optimize
  std::unique_ptr<QMCCostFunctionBase> optTarget;

  ///vmc engine
  std::unique_ptr<VMCBatched> vmcEngine;

  VMCDriverInput vmcdriver_input_;
  SampleStack& samples_;

  /// This is retained in order to construct and reconstruct the vmcEngine.
  const std::optional<EstimatorManagerInput> global_emi_;
};
} // namespace qmcplusplus
#endif
