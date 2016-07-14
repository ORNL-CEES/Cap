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
template class MPValuesParameters<3>;
template class SuperCapacitorMPValues<2>;
template class SuperCapacitorMPValues<3>;
template class InhomogeneousSuperCapacitorMPValues<2>;
template class InhomogeneousSuperCapacitorMPValues<3>;
template class SuperCapacitorMPValuesFactory<2>;
template class SuperCapacitorMPValuesFactory<3>;
template class PorousElectrodeMPValues<2>;
template class PorousElectrodeMPValues<3>;
template class MetalFoilMPValues<2>;
template class MetalFoilMPValues<3>;

template class CompositeMat<2>;
template class CompositeMat<3>;
template class CompositePro<2>;
template class CompositePro<3>;
template class UniformConstantMPValues<2>;
template class UniformConstantMPValues<3>;

} // end samespace cap
