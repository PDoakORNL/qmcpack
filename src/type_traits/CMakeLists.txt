#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2021 QMCPACK developers.
#//
#// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
#//
#// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
#//////////////////////////////////////////////////////////////////////////////////////

add_library(qmc_type_traits INTERFACE)
target_include_directories(qmc_type_traits INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})

if(BUILD_UNIT_TESTS)
  add_subdirectory(tests)
endif()
