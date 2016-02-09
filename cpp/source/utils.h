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

template <class InputIterator1, class InputIterator2, class OutputIterator,
          class T>
void approximate_integral_with_trapezoidal_rule(InputIterator1 first1,
                                                InputIterator1 last1,
                                                InputIterator2 first2,
                                                OutputIterator result, T init)
{
  *result = init;
  ++first1;
  ++first2;
  ++result;
  while (first1 != last1)
  {
    *result =
        *std::prev(result) +
        0.5 * (*first1 - *std::prev(first1)) * (*first2 + *std::prev(first2));
    ++first1;
    ++first2;
    ++result;
  }
}

template <typename T>
std::vector<T> to_vector(std::string const &s);

template <typename T>
std::string to_string(std::vector<T> const &v);

// template <typename T>
// void traverse_recursive
//    ( boost::property_tree::ptree const &            parent
//    , boost::property_tree::ptree::path_type const & child_path
//    , boost::property_tree::ptree const &            child
//    , T & method
//    )
//{
//    method(parent, child_path, child);
//    for (boost::property_tree::ptree::const_iterator it = child.begin(); it !=
//    child.end(); ++it) {
//        boost::property_tree::ptree::path_type current_path = child_path /
//        boost::property_tree::ptree::path_type(it->first);
//        traverse_recursive(parent, current_path, it->second, method);
//    }
//}
//
// template <typename T>
// void traverse(const boost::property_tree::ptree & parent, T & method)
//{
//  traverse_recursive(parent, "", parent, method);
//}
//
// void merge
//    ( boost::property_tree::ptree /* const */ &            parent
//    , boost::property_tree::ptree::path_type const & child_path
//    , boost::property_tree::ptree const &            child
//    )
//{
//    parent.put(child_path, child.data());
//}

} // end namespace cap

#endif // CAP_UTILS_H
