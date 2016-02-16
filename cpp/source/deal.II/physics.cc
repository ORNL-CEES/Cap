/* Copyright (c) 2016, the Cap authors.
 *
 * This file is subject to the Modified BSD License and may not be distributed
 * without copyright and license information. Please refer to the file LICENSE
 * for the text and further information on this license.
 */

#include <cap/deal.II/physics.templates.h>

namespace cap
{
template class PhysicsParameters<2>;
template class PhysicsParameters<3>;
template class Physics<2>;
template class Physics<3>;
}
