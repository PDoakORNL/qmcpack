//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2023 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


/** \file
 *  This implements micro benchmarking on RandomGenerator
 */

#include "catch.hpp"

#include <algorithm>

#include <OhmmsPETE/TinyVector.h>
#include <MCCoords.hpp>
#include "StdRandom.h"
#include <ParticleBase/RandomSeqGenerator.h>

namespace qmcplusplus
{

// Mechanism to pretty print benchmark names.
struct RandomBenchmarkParameters;
std::ostream& operator<<(std::ostream& out, const RandomBenchmarkParameters& dcbmp);
struct RandomBenchmarkParameters
{
  std::string name;
  int n_walkers;
  int batch_size;
  int steps;
  int particles;

  std::string str()
  {
    std::stringstream stream;
    stream << *this;
    return stream.str();
  }
};

std::ostream& operator<<(std::ostream& out, const RandomBenchmarkParameters& rbp)
{
  out << rbp.name << " n=" << rbp.n_walkers << " batch_size=" << rbp.batch_size << " steps=" << rbp.steps << " particles=" << rbp.particles;
  return out;
}
  
/** This benchmark runs by default
 */
TEST_CASE("Benchmark of makeGaussRandomWithEngine single thread", "[utilities][random][benchmark]")
{
  RandomBenchmarkParameters params;
  params.name       = "Batched Paticle Move";
  params.particles          = 32;
  params.n_walkers          = 256;
  params.steps      = 128;
  params.batch_size    = 128;
  
  StdRandom<double> std_rng;

  MCCoords<CoordsType::POS> deltas(params.particles * params.batch_size);
  int threads = 1;
  
  int batches = params.n_walkers / params.batch_size;

  BENCHMARK_ADVANCED(params.str())(Catch::Benchmark::Chronometer meter)
  {
    meter.measure([&] {
  for (int batch = 0; batch < batches / threads; ++batch)
    for (int step = 0; step < params.steps; ++step)
      {
	makeGaussRandomWithEngine(deltas, std_rng);
      }
    });
  };
}


} // namespace qmcplusplus
