#include <cap/utils.h>
#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>

namespace cap {

template <typename T>
std::vector<T> to_vector(std::string const & s)
{
    std::vector<T> v;
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, ',')) { v.push_back(boost::lexical_cast<T>(item)); }
    return v;
}

template <typename T>
std::string to_string(std::vector<T> const & v)
{
    std::string s;
    BOOST_FOREACH(T const & item, v) { s.append(boost::lexical_cast<std::string>(item)+","); }
    return s;
}

} // end namespace cap

