/* Copyright (c) 2016, the Cap authors.
 *
 * This file is subject to the Modified BSD License and may not be distributed
 * without copyright and license information. Please refer to the file LICENSE
 * for the text and further information on this license.
 */

#include <boost/property_tree/ptree.hpp>
#include <memory>

namespace cap
{

void compute_equivalent_circuit(
    boost::property_tree::ptree const &input_database,
    boost::property_tree::ptree &output_database);

} // end namespace cap
