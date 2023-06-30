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
#include "StdRandomNoBase.hpp"
#include "FakeRandom.h"
#include "RandomBase.h"
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
  out << rbp.name << " nw:" << rbp.n_walkers << " stp:" << rbp.steps << " ps:" << rbp.particles;
  return out;
}

/** This benchmark runs by default
 */
TEST_CASE("makeGaussRandomWithEngine TinyVec bsize 16 1thread StdRandomNoBase", "[utilities][random][benchmark]")
{
  RandomBenchmarkParameters params;
  params.name       = "Std Batched";
  params.particles  = 64;
  params.n_walkers  = 64;
  params.steps      = 128;
  params.batch_size = 16;

  StdRandomNoBase<double> std_rng;

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

  params.particles = 128;
  MCCoords<CoordsType::POS> deltas2(params.particles * params.batch_size);

  BENCHMARK_ADVANCED(params.str())(Catch::Benchmark::Chronometer meter)
  {
    meter.measure([&] {
      for (int batch = 0; batch < batches / threads; ++batch)
        for (int step = 0; step < params.steps; ++step)
        {
          makeGaussRandomWithEngine(deltas2, std_rng);
        }
    });
  };

  params.particles = 256;
  MCCoords<CoordsType::POS> deltas3(params.particles * params.batch_size);
  BENCHMARK_ADVANCED(params.str())(Catch::Benchmark::Chronometer meter)
  {
    meter.measure([&] {
      for (int batch = 0; batch < batches / threads; ++batch)
        for (int step = 0; step < params.steps; ++step)
        {
          makeGaussRandomWithEngine(deltas3, std_rng);
        }
    });
  };
}

TEST_CASE("makeGaussRandomWithEngine TinyVec bsize 16 1thread 2 RandomBase child", "[utilities][random][benchmark]")
{
  RandomBenchmarkParameters params;
  params.name       = "Std Batched";
  params.particles  = 64;
  params.n_walkers  = 64;
  params.steps      = 128;
  params.batch_size = 16;

  StdRandom<double> std_rng;
  FakeRandom<double> fake_rng;
  RandomBase<double>& base_rng = std_rng;

  MCCoords<CoordsType::POS> deltas(params.particles * params.batch_size);
  int threads = 1;

  int batches = params.n_walkers / params.batch_size;

  BENCHMARK_ADVANCED(params.str())(Catch::Benchmark::Chronometer meter)
  {
    meter.measure([&] {
      for (int batch = 0; batch < batches / threads; ++batch)
        for (int step = 0; step < params.steps; ++step)
        {
          makeGaussRandomWithEngine(deltas, base_rng);
        }
    });
  };

  params.particles = 128;
  MCCoords<CoordsType::POS> deltas2(params.particles * params.batch_size);

  BENCHMARK_ADVANCED(params.str())(Catch::Benchmark::Chronometer meter)
  {
    meter.measure([&] {
      for (int batch = 0; batch < batches / threads; ++batch)
        for (int step = 0; step < params.steps; ++step)
        {
          makeGaussRandomWithEngine(deltas2, base_rng);
        }
    });
  };

  params.particles = 256;
  MCCoords<CoordsType::POS> deltas3(params.particles * params.batch_size);
  BENCHMARK_ADVANCED(params.str())(Catch::Benchmark::Chronometer meter)
  {
    meter.measure([&] {
      for (int batch = 0; batch < batches / threads; ++batch)
        for (int step = 0; step < params.steps; ++step)
        {
          makeGaussRandomWithEngine(deltas3, base_rng);
        }
    });
  };

  
  RandomBase<double>& base_rng2 = fake_rng;
  params.name                   = "Fake Batched";

  BENCHMARK_ADVANCED(params.str())(Catch::Benchmark::Chronometer meter)
  {
    meter.measure([&] {
      for (int batch = 0; batch < batches / threads; ++batch)
        for (int step = 0; step < params.steps; ++step)
        {
          makeGaussRandomWithEngine(deltas, base_rng2);
        }
    });
  };

  params.particles = 128;
  BENCHMARK_ADVANCED(params.str())(Catch::Benchmark::Chronometer meter)
  {
    meter.measure([&] {
      for (int batch = 0; batch < batches / threads; ++batch)
        for (int step = 0; step < params.steps; ++step)
        {
          makeGaussRandomWithEngine(deltas2, base_rng2);
        }
    });
  };

  params.particles = 256;
  BENCHMARK_ADVANCED(params.str())(Catch::Benchmark::Chronometer meter)
  {
    meter.measure([&] {
      for (int batch = 0; batch < batches / threads; ++batch)
        for (int step = 0; step < params.steps; ++step)
        {
          makeGaussRandomWithEngine(deltas3, base_rng2);
        }
    });
  };

}

template<class RNG>
bool dummyRandomFunc(RNG& rng, int electrons)
{
  using Real = typename RNG::result_type;
  Real compare_to = 0.5;
  int accum{0};
  for (int i = 0; i < electrons; ++i)
    if (rng() > compare_to)
      accum += 1;
    else
      accum -= 1;
  return static_cast<bool>(accum % 2);
}

template<class RNG>
auto dummyRandomRot(RNG& rng, int ions)
{
  using Real = typename RNG::result_type;
  Real x{0};
  Real y{0};
  Real z{0};
  for (int i = 0; i < ions; ++i)
  {
    x += rng();
    y += rng();
    z += rng();
  }
  return x + y + z;
}

TEST_CASE("Nonlocal single StdRandomNoBase", "[utilities][random][benchmark]")
{
  RandomBenchmarkParameters params;
  params.name       = "Std Nonlocal";
  params.particles  = 256;
  params.n_walkers  = 256;
  params.steps      = 128;
  params.batch_size = 128;

  StdRandomNoBase<double> std_rng;

  MCCoords<CoordsType::POS> deltas(params.particles * params.batch_size);
  int threads = 1;

  int batches = params.n_walkers / params.batch_size;

  BENCHMARK_ADVANCED(params.str())(Catch::Benchmark::Chronometer meter)
  {
    meter.measure([&] {
      for (int batch = 0; batch < batches / threads; ++batch)
        for (int step = 0; step < params.steps; ++step)
        {
          bool res = dummyRandomFunc(std_rng, params.particles);
          auto sum = dummyRandomRot(std_rng, params.particles / 4);
        }
    });
  };
}

TEST_CASE("Nonlocal single thread 2 child", "[utilities][random][benchmark]")
{
  RandomBenchmarkParameters params;
  params.name       = "Std Nonlocal";
  params.particles  = 256;
  params.n_walkers  = 256;
  params.steps      = 128;
  params.batch_size = 128;

  StdRandom<double> std_rng;
  FakeRandom<double> fake_rng;
  RandomBase<double>& base_rng = std_rng;

  MCCoords<CoordsType::POS> deltas(params.particles * params.batch_size);
  int threads = 1;

  int batches = params.n_walkers / params.batch_size;

  BENCHMARK_ADVANCED(params.str())(Catch::Benchmark::Chronometer meter)
  {
    meter.measure([&] {
      for (int batch = 0; batch < batches / threads; ++batch)
        for (int step = 0; step < params.steps; ++step)
        {
          bool res = dummyRandomFunc(base_rng, params.particles);
          auto sum = dummyRandomRot(base_rng, params.particles / 4);
        }
    });
  };

  RandomBase<double>& base_rng2 = fake_rng;
  params.name                   = "Fake Nonlocal";

  BENCHMARK_ADVANCED(params.str())(Catch::Benchmark::Chronometer meter)
  {
    meter.measure([&] {
      for (int batch = 0; batch < batches / threads; ++batch)
        for (int step = 0; step < params.steps; ++step)
        {
          bool res = dummyRandomFunc(base_rng, params.particles);
          auto sum = dummyRandomRot(base_rng, params.particles / 4);
        }
    });
  };
}

} // namespace qmcplusplus
