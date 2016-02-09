/* Copyright (c) 2016, the Cap authors.
 *
 * This file is subject to the Modified BSD License and may not be distributed
 * without copyright and license information. Please refer to the file LICENSE
 * for the text and further information on this license. 
 */

#include <cap/boundary_values.templates.h>

namespace cap
{

template class BoundaryValuesParameters<2>;
template class BoundaryValues<2>;
template class SuperCapacitorBoundaryValuesParameters<2>;
template class SuperCapacitorBoundaryValues<2>;
template class BoundaryValuesParameters<3>;
template class BoundaryValues<3>;
template class SuperCapacitorBoundaryValuesParameters<3>;
template class SuperCapacitorBoundaryValues<3>;

} // end namespace cap
