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
    {
      boost::algorithm::to_lower(item);
      if (item.compare("true") == 0)
        item = "1";
      else if (item.compare("false") == 0)
        item = "0";
    }
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

} // end namespace cap
