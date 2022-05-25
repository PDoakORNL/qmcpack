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

#include "OhmmsPETE/OhmmsMatrix.h"
#include "OhmmsPETE/OhmmsVector.h"

namespace qmcplusplus
{

template<typename REAL>
class ListenerVar
{
public:
  ListenerVar(const std::string& name);
  std::function<void(const int walker_index, const std::vector<REAL>&)> report;
  const std::string& get_name() const { return name_; }
private:
  const std::string name_;  
};

template<typename REAL>
class ListenerVector
{
public:
  ListenerVector(const std::string& name, std::function<void(const int walker_index, const Vector<REAL>&)>report_func) : report(report_func), name_(name) {}
  std::function<void(const int walker_index, const Vector<REAL>&)> report;
  const std::string& get_name() const { return name_; }
private:
  const std::string name_;
};


template<typename REAL>
class ListenerCombined
{
public:
  ListenerCombined(const std::string& name, std::function<void(const int walker_index, const Matrix<REAL>&)>report_func);
  std::function<void(const int walker_index, const std::string& name, const Matrix<REAL>&)> report;
  const std::string& get_name() const { return name_; }
private:
  const std::string name_;
};


template<typename REAL>
class Listener
{
public:
  std::vector<std::string> vars;
  std::vector<std::string> multi_vars;
  std::function<void(const std::vector<REAL>&)> report_vars;
  std::function<void(const std::string& name, Vector<REAL>&)> report_component;
  std::function<void(const Matrix<REAL>&)> report_multi_vars;
};


} // namespace qmcplusplus


#endif
