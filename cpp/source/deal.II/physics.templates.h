/* Copyright (c) 2016, the Cap authors.
 *
 * This file is subject to the Modified BSD License and may not be distributed
 * without copyright and license information. Please refer to the file LICENSE
 * for the text and further information on this license.
 */

#ifndef CAP_PHYSICS_TEMPLATES_H
#define CAP_PHYSICS_TEMPLATES_H

#include <cap/physics.h>
#include <deal.II/base/function.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/vector_tools.h>

namespace cap
{
template <int dim>
Physics<dim>::Physics(std::shared_ptr<PhysicsParameters<dim> const> parameters)
    : dof_handler(parameters->dof_handler), mp_values(parameters->mp_values)
{
}
}

#endif
