/* Copyright (c) 2016, the Cap authors.
 *
 * This file is subject to the Modified BSD License and may not be distributed
 * without copyright and license information. Please refer to the file LICENSE
 * for the text and further information on this license.
 */

#include <cap/operator.h>

#include <deal.II/dofs/dof_tools.h>

namespace cap
{

//////////////////////// OPERATOR ////////////////////////////
template <int dim>
Operator<dim>::Operator(
    std::shared_ptr<OperatorParameters<dim> const> parameters)
    : dof_handler(parameters->dof_handler),
      constraint_matrix(parameters->constraint_matrix),
      mp_values(parameters->mp_values), b_values(parameters->boundary_values)
{
  this->stiffness_matrix.reinit(*(parameters->sparsity_pattern));
  this->mass_matrix.reinit(*(parameters->sparsity_pattern));
  this->load_vector.reinit(*(parameters->some_vector));
}

template <int dim>
void Operator<dim>::set_null_space(unsigned int const component,
                                   dealii::types::material_id const material_id)
{
  unsigned int const n_components =
      dealii::DoFTools::n_components(*dof_handler);
  std::vector<bool> mask(n_components, false);
  mask[component] = true;
  dealii::ComponentMask component_mask(mask);
  std::vector<bool> selected_dofs((*dof_handler).n_dofs());
  dealii::DoFTools::extract_dofs(*dof_handler, component_mask, selected_dofs);
  dealii::FiniteElement<dim> const &fe = (*dof_handler).get_fe();
  unsigned int const dofs_per_cell = fe.dofs_per_cell;
  std::vector<dealii::types::global_dof_index> local_dof_indices(dofs_per_cell);
  for (auto cell : dof_handler->active_cell_iterators())
  {
    cell->get_dof_indices(local_dof_indices);
    if (cell->material_id() != material_id)
    {
      for (unsigned int dof = 0; dof < dofs_per_cell; ++dof)
      {
        selected_dofs[local_dof_indices[dof]] = false;
      } // end for dof
    }   // end if
  }     // end for cell
  for (dealii::types::global_dof_index dof = 0; dof < (*dof_handler).n_dofs();
       ++dof)
  {
    if (selected_dofs[dof])
    {
      null_space_dof_indices.push_back(dof - this->dof_shift);
    } // end if
  }   // end for dof
}

} // end namespace cap
