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

#include <catch.hpp>
#include "CUDA/CUDALinearAlgebraHandles.h"

namespace qmcplusplus
{
TEST_CASE("CUDALinearAlgebraHandles_copy_destroy", "[CUDA]")
{
  CUDALinearAlgebraHandles clah;
  {
    CUDALinearAlgebraHandles clah2(clah);
    CHECK(clah.getStream() == clah2.getStream());
  }
  cudaStream_t cublas_stream;
  cublasErrorCheck(cublasGetStream(clah.getCuBLAS(), &cublas_stream), "cublasGetStream failed!");
  CHECK(clah.getStream() == cublas_stream);
  // won't crash if there is no access violating second destroy.
}
} // namespace qmcplusplus
