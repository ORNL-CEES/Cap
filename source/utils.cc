#include <cache/utils.templates.h>

namespace cache {

template std::vector<int>         to_vector(std::string const & s);
template std::vector<float>       to_vector(std::string const & s);
template std::vector<double>      to_vector(std::string const & s);
template std::vector<std::string> to_vector(std::string const & s);

template std::string to_string(std::vector<int> const         & v);
template std::string to_string(std::vector<float> const       & v);
template std::string to_string(std::vector<double> const      & v);
template std::string to_string(std::vector<std::string> const & v);

} // end namespace cache

