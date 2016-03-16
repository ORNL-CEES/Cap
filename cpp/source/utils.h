/* Copyright (c) 2016, the Cap authors.
 *
 * This file is subject to the Modified BSD License and may not be distributed
 * without copyright and license information. Please refer to the file LICENSE
 * for the text and further information on this license.
 */

#ifndef CAP_UTILS_H
#define CAP_UTILS_H

#include <boost/property_tree/ptree.hpp>
#include <string>
#include <vector>
#include <iterator>

namespace cap
{
template <typename T>
std::vector<T> to_vector(std::string const &s);

template <typename T>
std::string to_string(std::vector<T> const &v);
} // end namespace cap

#endif // CAP_UTILS_H
