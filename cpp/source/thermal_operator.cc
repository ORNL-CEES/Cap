/* Copyright (c) 2016, the Cap authors.
 *
 * This file is subject to the Modified BSD License and may not be distributed
 * without copyright and license information. Please refer to the file LICENSE
 * for the text and further information on this license.
 */

#include <cap/thermal_operator.templates.h>

namespace cap
{

template class ThermalOperatorParameters<2>;
template class ThermalOperator<2>;
template class ThermalOperatorParameters<3>;
template class ThermalOperator<3>;

} // end namespace cap
