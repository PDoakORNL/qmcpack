//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
//////////////////////////////////////////////////////////////////////////////////////

/** @file
 *  A minimally functional wrapper for the since c++11 <random>
 */
#ifndef QMCPLUSPLUS_STDRANDNOBASE_HPP
#define QMCPLUSPLUS_STDRANDNOBASE_HPP

#include <vector>
#include <random>
#include <sstream>
#include <iterator>
#include "StdRandom.h"
namespace qmcplusplus
{

/** generating real type random numbers [min, max) in a uniform distribution.
 * This is intended to match the behavior of random::uniform_real_distribution in the boost libraries
 * Its behavor is different from std::uniform_real_distribution which uses std::generate_canonical
 * to generate enough entropy according to the precision of T.
 */
// template<typename T>
// class uniform_real_distribution_as_boost
// {
// public:
//   using result_type = T;
//   static_assert(std::is_floating_point<T>::value);

//   uniform_real_distribution_as_boost(T min = T(0.0), T max = T(1.0)) : min_(min), max_(max) {}
//   ///Generating functions.
//   template<typename RNG>
//   result_type operator()(RNG& eng)
//   {
//     return static_cast<result_type>(eng() - eng.min()) / (static_cast<result_type>(eng.max() - eng.min()) + 1) *
//         (max_ - min_) + min_;
//   }

// private:
//   T min_;
//   T max_;
// };

template<typename T>
class StdRandomNoBase
{
public:
  using Engine      = std::mt19937;
  using result_type = T;
  using uint_type   = uint_fast32_t;

  StdRandomNoBase(uint_type iseed = 911);

  void init(int iseed_in)
  {
    uint_type baseSeed = iseed_in;
    engine.seed(baseSeed);
  }

  void seed(uint_type aseed)  { engine.seed(aseed); }
  result_type operator()() ;
  void write(std::ostream& rout) const  { rout << engine; }
  void read(std::istream& rin)  { rin >> engine; }
  size_t state_size() const  { return stream_state_size; }

  void load(const std::vector<uint_type>& newstate) ;
  void save(std::vector<uint_type>& curstate) const ;
  std::unique_ptr<StdRandomNoBase<T>> clone() const  { return std::make_unique<StdRandomNoBase<T>>(*this); }

  // Non const allows use of default copy constructor
  std::string ClassName{"StdRand"};
  std::string EngineName{"std::mt19937"};

private:
  ///random number generator [0,1)
  uniform_real_distribution_as_boost<T> distribution;
  Engine engine;
  /// the number count of streaming states. Must match read/write/load/save
  std::size_t stream_state_size;
};

} // namespace qmcplusplus

#endif
