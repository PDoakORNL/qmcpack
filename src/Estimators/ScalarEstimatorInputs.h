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

#ifndef QMCPLUSPLUS_SCALAR_ESTIMATOR_INPUTS_H
#define QMCPLUSPLUS_SCALAR_ESTIMATOR_INPUTS_H

#include "InputSection.h"

namespace qmcplusplus
{

class LocalEnergyEstimator;
  
class LocalEnergyInput
{
  class LocalEnergyInputSection : public InputSection
  {
  public:
    LocalEnergyInputSection()
    {
      section_name = "LocalEnergy";
      attributes   = {"type", "hdf5"};
      bools        = {"hdf5"};
      strings      = {"type"};
    };
  };
public:
  using Consumer = LocalEnergyEstimator;
  LocalEnergyInput() = default;
  LocalEnergyInput(xmlNodePtr cur);
  bool get_use_hdf5() const { return use_hdf5_; }
private:
  LocalEnergyInputSection input_section_;
  bool use_hdf5_ = true;
};

struct CSEnergyEstimator;

class CSLocalEnergyInput
{
  class CSLocalEnergyInputSection : public InputSection
  {
  public:
    CSLocalEnergyInputSection()
    {
      section_name = "CSLocalEnergy";
      attributes   = {"npsi", "type"};
      integers     = {"npsi"};
      strings      = {"type"};
    }
  };
public:
  using Consumer = CSEnergyEstimator;
  CSLocalEnergyInput() = default;
  CSLocalEnergyInput(xmlNodePtr cur);
  int get_n_psi() const { return n_psi_; }

private:
  CSLocalEnergyInputSection input_section_;
  int n_psi_ = 1;
};

} // namespace qmcplusplus

#endif
