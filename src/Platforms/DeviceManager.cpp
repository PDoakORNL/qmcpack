//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
//////////////////////////////////////////////////////////////////////////////////////


#include "DeviceManager.h"
#include <memory>
#include "Host/OutputManager.h"

namespace qmcplusplus
{

DeviceManager::DeviceManager(int local_rank, int local_size)
    : default_device_num(-1),
      num_devices(0)
#if defined(ENABLE_CUDA)
      ,
      cuda_dm_(default_device_num, num_devices, local_rank, local_size)
#endif
#if defined(ENABLE_OFFLOAD)
      ,
      omptarget_dm_(default_device_num, num_devices, local_rank, local_size)
#endif
{
  if (num_devices > 0)
  {
    if (local_size % num_devices != 0)
      app_warning() << "The number of MPI ranks on the node is not divisible by the number of accelerators. "
                    << "Imbalanced load may cause performance loss.\n";
  }
}

DeviceManager::~DeviceManager() = default;

std::unique_ptr<DeviceManager> DeviceManager::global;

} // namespace qmcplusplus
