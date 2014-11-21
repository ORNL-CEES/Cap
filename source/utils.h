#ifndef CACHE_UTILS_H
#define CACHE_UTILS_H

#include <string>
#include <vector>

namespace cache {

template <typename T>
std::vector<T> to_vector(std::string const & s);

template <typename T>
std::string to_string(std::vector<T> const & v);

} // end namespace cache

#endif // CACHE_UTILS_H
