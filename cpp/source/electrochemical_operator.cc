/* Copyright (c) 2016, the Cap authors.
 *
 * This file is subject to the Modified BSD License and may not be distributed
 * without copyright and license information. Please refer to the file LICENSE
 * for the text and further information on this license.
 */

#include <cap/electrochemical_operator.templates.h>

namespace cap
{

template class ElectrochemicalOperatorParameters<2>;
template class ElectrochemicalOperator<2>;
template class ElectrochemicalOperatorParameters<3>;
template class ElectrochemicalOperator<3>;

} // end namespace cap
