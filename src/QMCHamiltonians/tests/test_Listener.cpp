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

#include "catch.hpp"

#include "Listener.hpp"
#include <string>
#include <vector>
#include "Configuration.h"

namespace qmcplusplus
{

using QMCT = QMCTraits;
using Real = QMCT::RealType;

namespace testing
{

class Talker
{
private:
  std::vector<Listener<Real>> listeners;
  std::vector<ListenerVector<Real>> listener_vectors;
public:
  void speakTo(Listener<Real> listener) { listeners.push_back(listener); }
  void registerVector(ListenerVector<Real>& listener_vector) { listener_vectors.push_back(listener_vector); }
  void report()
  {
    int listener_count = 0;
    for (auto& listener : listeners) {
      std::vector<Real> vals(listener.vars.size());
      std::iota(vals.begin(), vals.end(), listener_count); 
      listener.report_vars(vals);
      Matrix<Real> matrix(3, listener_count + 1);
      std::iota(matrix.begin(), matrix.end(), listener_count);
      listener.report_multi_vars(matrix);
      ++listener_count;      
    }
    
  }
  void reportVector() {
    Vector<Real> vec_part(4);
    std::iota(vec_part.begin(), vec_part.end(), 0);
    listener_vectors[0].report(0,vec_part);
  }
    
};

class TestReceiver
{
public:
  void doRegister(Talker& talker)
  {
    std::vector<std::string> variable_tags{"test","test2"};
    std::vector<std::string> multi_variable_tags{"mtest"};
    values_.resize(variable_tags.size());
    auto takeReport = [this](const std::vector<Real>& values) { for(int i=0; i < values.size(); ++i)
	values_[i] = values[i];};
    auto takeComponentReport = [this](const std::string& name, Vector<Real>& values) { component_values_ = values; };
    auto takeMultiVarReport = [this](const Matrix<Real>& values) { multi_vars_ = values; };
    Listener<Real> my_listener{variable_tags, multi_variable_tags, takeReport, takeComponentReport, takeMultiVarReport};
    talker.speakTo(my_listener);
  }

  auto getParticularListener(Vector<Real>& part_val) {
    return [&part_val] (const int walker_index, const Vector<Real>& values) ->void { part_val = values; };
  }
  
  Vector<Real> component_values_;
  std::vector<Real> values_;
  Matrix<Real> multi_vars_;

  Vector<Real> particular_values_;
};

TEST_CASE("Listener::basic", "[hamiltonian]")
{
  Talker talker;
  TestReceiver test_listener;
  test_listener.doRegister(talker);
  talker.report();
  CHECK(test_listener.values_[0] == 0);
  CHECK(test_listener.values_[1] == 1);

  TestReceiver test_listener2;
  test_listener2.doRegister(talker);
  talker.report();

  CHECK(test_listener.values_[0] == Approx(0.0));
  CHECK(test_listener.values_[1] == Approx(1.0));
  CHECK(test_listener2.values_[0] ==Approx(1));
  CHECK(test_listener2.values_[1] == Approx(2));

  CHECK(test_listener.multi_vars_(0,0) == Approx(0));
  CHECK(test_listener.multi_vars_(1,0) == Approx(1));
  CHECK(test_listener.multi_vars_(2,0) == Approx(2));

  CHECK(test_listener2.multi_vars_(0,0) == Approx(1));
  CHECK(test_listener2.multi_vars_(1,0) == Approx(3));
  CHECK(test_listener2.multi_vars_(2,0) == Approx(5));
}

TEST_CASE("ListenerVector", "[hamiltonian]")
{
  Talker talker;
  TestReceiver test_receiver;
  ListenerVector<Real> listen_vector("weight", test_receiver.getParticularListener(test_receiver.particular_values_));
  
  talker.registerVector(listen_vector);

  talker.reportVector();
  CHECK(test_receiver.particular_values_[0] == 0);
  CHECK(test_receiver.particular_values_[3] == 3);
}

} // namespace testing
} // namespace qmcplusplus
