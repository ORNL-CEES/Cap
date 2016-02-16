/* Copyright (c) 2016, the Cap authors.
 *
 * This file is subject to the Modified BSD License and may not be distributed
 * without copyright and license information. Please refer to the file LICENSE
 * for the text and further information on this license.
 */

#ifndef CAP_PHYSICS_TEMPLATES_H
#define CAP_PHYSICS_TEMPLATES_H

#include <cap/deal.II/physics.h>
#include <deal.II/base/function.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/vector_tools.h>

namespace cap
{
template <int dim>
Physics<dim>::Physics(std::shared_ptr<PhysicsParameters<dim> const> parameters)
    : dof_handler(parameters->dof_handler), 
      n_dofs(parameters->n_dofs),
      mp_values(parameters->mp_values)
{
  dealii::DoFTools::make_hanging_node_constraints(*dof_handler, constraint_matrix);
  // TODO temporary Dirichlet boundary condition.
  dealii::VectorTools::interpolate_boundary_values(*dof_handler, 0, 
      dealii::ZeroFunction<dim>(), constraint_matrix);
  constraint_matrix.close();

  // Create sparsity pattern
  unsigned int const max_couplings = dof_handler->max_couplings_between_dofs();
  sparsity_pattern.reinit(n_dofs, n_dofs,
                          max_couplings);
  dealii::DoFTools::make_sparsity_pattern(
      *dof_handler, sparsity_pattern, constraint_matrix);
  sparsity_pattern.compress();

  // initialize matrices and vectors
  system_matrix.reinit(sparsity_pattern);
  mass_matrix.reinit(sparsity_pattern);
  system_rhs.reinit(n_dofs);
}
}

#endif
