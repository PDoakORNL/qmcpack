//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_LISTENER_HPP
#define QMCPLUSPLUS_LISTENER_HPP
/** \file
 *  Listener and its supporting types serves to deliver "trace" values from
 *  QMCHamiltonian to estimators.
 *  It is designed to be lightweight and minimal but while still minimizing the
 *  changes in QMCHamiltonian at this time.
 */

#include <functional>

#include "type_traits/template_types.hpp"

#include "OhmmsPETE/OhmmsMatrix.h"
#include "OhmmsPETE/OhmmsVector.h"

namespace qmcplusplus
{

template<typename REAL>
class ListenerVar
{
public:
  ListenerVar(const std::string& name);
  std::function<void(const int walker_index, const std::string& name, const std::vector<REAL>&)> report;
  const std::string& get_name() const { return name_; }

private:
  const std::string name_;
};


/** This is the type registers a listener expecting a callback with a vector of values 1 per particle, called once per walker.
 *  The name is primarily for debugging purposes, if have decided against have the QMCHamiltonian use it to estable routeing.
 *  Instead the register functions are specfic for what the listener wants to listen to.
 */
template<typename REAL>
class ListenerVector
{
public:
  using Report = std::function<void(const int walker_index, const std::string& name, const Vector<REAL>&)>;
  ListenerVector(const std::string& name,
                 Report report_func)
      : report(report_func), name_(name)
  {}
  /** Report vector for a Hamiltonian to make a per particle report to a listener
   */
  Report report;
  const std::string& get_name() const { return name_; }

private:
  const std::string name_;
};

/** Its unclear there is any need for this.
 */
template<typename REAL>
class ListenerCombined
{
public:
  ListenerCombined(
      const std::string& name,
      std::function<void(const int walker_index, const RefVector<const std::string&> names, const Matrix<REAL>&)>
          report_func);
  std::function<void(const int walker_index, const RefVector<const std::string&> name, const Matrix<REAL>&)> report;
  const std::string& get_name() const { return name_; }

private:
  const std::string name_;
};

} // namespace qmcplusplus


#endif
