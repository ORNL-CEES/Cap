/* Copyright (c) 2016, the Cap authors.
 *
 * This file is subject to the Modified BSD License and may not be distributed
 * without copyright and license information. Please refer to the file LICENSE
 * for the text and further information on this license. 
 */

#include <cap/geometry.templates.h>

namespace cap
{

template class Geometry<2>;
template class Geometry<3>;
template class SuperCapacitorGeometry<2>;
template class SuperCapacitorGeometry<3>;
template class DummyGeometry<2>;
template class DummyGeometry<3>;

} // end namespace cap
