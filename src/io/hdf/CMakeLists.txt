#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2020 QMCPACK developers.
#//
#// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
#//
#// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
#//////////////////////////////////////////////////////////////////////////////////////

set(IO_HDF5_SRC hdf_wrapper_functions.cpp hdf_archive.cpp hdf_path.cpp hdf_pete_hyperslab.cpp)

add_library(qmcio_hdf ${IO_HDF5_SRC})

#TARGET_INCLUDE_DIRECTORIES(qmcio_hdf PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(qmcio_hdf PUBLIC containers IO::HDF5 message)

if(BUILD_UNIT_TESTS)
  add_subdirectory(tests)
endif()
