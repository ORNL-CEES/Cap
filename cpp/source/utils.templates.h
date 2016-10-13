/* Copyright (c) 2016, the Cap authors.
 *
 * This file is subject to the Modified BSD License and may not be distributed
 * without copyright and license information. Please refer to the file LICENSE
 * for the text and further information on this license.
 */

#include <cap/utils.h>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/foreach.hpp>
#include <type_traits>

namespace cap
{

namespace internal
{
void filter_boolean(std::string &b)
{
  boost::algorithm::to_lower(b);
  if (b.compare("true") == 0)
    b = "1";
  else if (b.compare("false") == 0)
    b = "0";
}
}

template <typename T>
std::vector<T> to_vector(std::string const &s)
{
  std::vector<T> v;
  std::stringstream ss(s);
  std::string item;
  while (std::getline(ss, item, ','))
  {
    boost::algorithm::trim(item);
    if (std::is_same<bool, T>::value)
      internal::filter_boolean(item);
    v.push_back(boost::lexical_cast<T>(item));
  }
  return v;
}

template <typename T>
std::string to_string(std::vector<T> const &v)
{
  std::string s;
  BOOST_FOREACH (T const &item, v)
  {
    s.append(boost::lexical_cast<std::string>(item) + ",");
  }
  return s;
}

template <typename T>
std::map<std::string, T> to_map(std::string const &s)
{
  std::map<std::string, T> m;
  std::vector<std::string> pairs;
  boost::algorithm::split(pairs, s, boost::algorithm::is_any_of(","));
  for (auto const &p : pairs)
  {
    std::vector<std::string> kv;
    boost::algorithm::split(kv, p, boost::algorithm::is_any_of("="));
    if (kv.size() != 2)
      throw std::runtime_error("invalid key-value pair " + p);
    for (auto &x : kv)
      boost::algorithm::trim(x);
    if (std::is_same<bool, T>::value)
      internal::filter_boolean(kv[1]);
    m.emplace(kv[0], boost::lexical_cast<T>(kv[1]));
  }
  return m;
}

} // end namespace cap
