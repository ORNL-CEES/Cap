#include <boost/property_tree/ptree.hpp>
#include <memory>

namespace cap {

void compute_equivalent_circuit(std::shared_ptr<boost::property_tree::ptree const> input_database,
                                std::shared_ptr<boost::property_tree::ptree      > output_database);

} // end namespace cap
