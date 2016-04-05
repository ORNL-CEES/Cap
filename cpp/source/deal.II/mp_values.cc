/* Copyright (c) 2016, the Cap authors.
 *
 * This file is subject to the Modified BSD License and may not be distributed
 * without copyright and license information. Please refer to the file LICENSE
 * for the text and further information on this license.
 */

#include <cap/mp_values.templates.h>

namespace cap
{

template class MPValuesParameters<2>;
template class MPValues<2>;
template class MPValuesParameters<3>;
template class MPValues<3>;
template class PorousElectrodeMPValues<2>;
template class PorousElectrodeMPValues<3>;
template class MetalFoilMPValues<2>;
template class MetalFoilMPValues<3>;

} // end samespace cap
