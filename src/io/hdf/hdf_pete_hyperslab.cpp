//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2025 QMCPACK developers
//
// File developed by: Peter W. Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter W. Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include "hdf_pete_hyperslab.h"

/**@file
 * @brief actual explicit instantiation of hyperslab_proxies for some OhmmsTypes
 *
 * I found this necessary to get tracable debug information for these
 * wrappers. It has the desirable side effect of reducing needless rebuilding
 * of io/hdf code.
 */

namespace qmcplusplus
{
#define HYPERSLAB_PROXY_INST(CONTAINER, RANK, INTE)                                                               \
  template struct hyperslab_proxy<CONTAINER, RANK>;                                                               \
  template hyperslab_proxy<CONTAINER, RANK>::hyperslab_proxy(CONTAINER& a, const std::array<INTE, RANK>& dims_in, \
                                                             const std::array<INTE, RANK>& selected_in,           \
                                                             const std::array<INTE, RANK>& offsets_in);           \
  template void hyperslab_proxy<CONTAINER, RANK>::adaptShape<INTE>(const std::vector<INTE>& sizes_file);          \
  template struct h5data_proxy<hyperslab_proxy<CONTAINER, RANK>>;

HYPERSLAB_PROXY_INST(Vector<double>, 1, std::size_t);
HYPERSLAB_PROXY_INST(Vector<double>, 2, std::size_t);
HYPERSLAB_PROXY_INST(Vector<double>, 3, std::size_t);
HYPERSLAB_PROXY_INST(Vector<double>, 4, std::size_t);
HYPERSLAB_PROXY_INST(Vector<std::complex<double>>, 1, std::size_t);
HYPERSLAB_PROXY_INST(Vector<std::complex<double>>, 2, std::size_t);
HYPERSLAB_PROXY_INST(Vector<std::complex<double>>, 3, std::size_t);
HYPERSLAB_PROXY_INST(Vector<std::complex<double>>, 4, std::size_t);

} // namespace qmcplusplus
