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
  std::vector<ListenerVector<Real>> listener_vectors_;
  const std::string name_{"Talker"};
public:
  void registerVector(ListenerVector<Real>& listener_vector) { listener_vectors_.push_back(listener_vector); }
  void reportVector() {
    Vector<Real> vec_part(4);
    std::iota(vec_part.begin(), vec_part.end(), 0);
    for (auto& listener : listener_vectors_)
      listener.report(0,name_,vec_part);
  }
    
};

class TestReceiver
{
public:
  /** Listener frunction that has captured an object data member.
   *  This leaves access to the listening object quick controlled but allows a great deal of flexibility in dealing with
   *  the vector report.  In this case we just copy it, but it could be a reduce or something more involved.
   */  
  auto getParticularListener(Vector<Real>& part_val) {
    return [&part_val] (const int walker_index, const std::string& name, const Vector<Real>& values) ->void { part_val = values; };
  }
  
  Vector<Real> particular_values_;
};

TEST_CASE("ListenerVector", "[hamiltonian]")
{
  Talker talker;
  TestReceiver test_receiver;
  ListenerVector<Real> listen_vector("kinetic", test_receiver.getParticularListener(test_receiver.particular_values_));
  
  talker.registerVector(listen_vector);

  talker.reportVector();
  CHECK(test_receiver.particular_values_[0] == 0);
  CHECK(test_receiver.particular_values_[3] == 3);
}

} // namespace testing
} // namespace qmcplusplus
