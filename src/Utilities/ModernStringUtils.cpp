//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//////////////////////////////////////////////////////////////////////////////////////

#include "ModernStringUtils.hpp"
#include <algorithm>
#include <cctype>

namespace qmcplusplus
{
using std::string_view;
using std::size_t;


std::string lowerCase(const std::string_view s)
{
  std::string lower_str{s};
  std::transform(lower_str.begin(), lower_str.end(), lower_str.begin(),
                 [](unsigned char c) { return std::tolower(c); });
  return lower_str;
}

namespace modernstrutil {
std::vector<std::string_view> split(const string_view s, const string_view delimiters)
{
  std::vector<std::string_view> tokens;
  size_t right = 0;
  size_t left = 0;
  while(true)
  {
    left = s.find_first_not_of(delimiters, right);
    right = s.find_first_of(delimiters, left);
    if(left != s.npos)
      tokens.push_back(s.substr(left, right - left));
    else
      break;
  }
  return tokens;
}
}


} // namespace qmcplusplus
