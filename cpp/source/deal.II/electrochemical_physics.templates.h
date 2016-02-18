/* Copyright (c) 2016, the Cap authors.
 *
 * This file is subject to the Modified BSD License and may not be distributed
 * without copyright and license information. Please refer to the file LICENSE
 * for the text and further information on this license.
 */

#ifndef CAP_DEAL_II_ELECTROCHEMICAL_PHYSICS_TEMPLATES_H
#define CAP_DEAL_II_ELECTROCHEMICAL_PHYSICS_TEMPLATES_H

#include <cap/deal.II/electrochemical_physics.h>
#include <boost/assert.hpp>
#include <deal.II/base/function.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/numerics/vector_tools.h>

namespace cap
{
template <int dim>
ElectrochemicalPhysics<dim>::ElectrochemicalPhysics(
    std::shared_ptr<PhysicsParameters<dim> const> parameters)
    : Physics<dim>(parameters)
{
  std::shared_ptr<boost::property_tree::ptree const> database =
      parameters->database;

  // clang-format off
  this->solid_potential_component  = database->get<unsigned int>("solid_potential_component");
  this->liquid_potential_component = database->get<unsigned int>("liquid_potential_component");
  this->temperature_component      = database->get<unsigned int>("temperature_component");
  // clang-format on

  alpha             = database->get<double>("material_properties.alpha", 0.0);
  anode_boundary_id = database->get<dealii::types::boundary_id>(
      "boundary_values.anode_boundary_id");
  cathode_boundary_id = database->get<dealii::types::boundary_id>(
      "boundary_values.cathode_boundary_id");

  std::shared_ptr<
      ElectrochemicalPhysicsParameters<dim> const> electrochemical_parameters =
      std::dynamic_pointer_cast<ElectrochemicalPhysicsParameters<dim> const>(
          parameters);
  BOOST_ASSERT_MSG(electrochemical_parameters == nullptr,
                   "Problem during dowcasting the pointer");

  // Take care of hanging nodes
  dealii::DoFTools::make_hanging_node_constraints(*(this->dof_handler),
                                                  this->constraint_matrix);

  // Take care of Dirichlet boundary condition.
  // The anode is always set in Earth (Dirichlet value of 0).
  // If we impose a the voltage, the cathode is also a Dirichlet condition.
  unsigned int const n_components =
      dealii::DoFTools::n_components(*(this->dof_handler));
  std::vector<bool> mask(n_components, false);
  mask[this->solid_potential_component] = true;
  dealii::ComponentMask component_mask(mask);
  typename dealii::FunctionMap<dim>::type dirichlet_boundary_condition;
  dealii::ZeroFunction<dim> homogeneous_bc(n_components);
  dirichlet_boundary_condition[anode_boundary_id] = &homogeneous_bc;
  dealii::Function<dim> *anode_dirichlet_bc = nullptr;
  if (electrochemical_parameters->charge_type == ConstantVoltage)
  {
    anode_dirichlet_bc = new dealii::ConstantFunction<dim>(
        electrochemical_parameters->constant_voltage, n_components);
    dirichlet_boundary_condition[cathode_boundary_id] = anode_dirichlet_bc;
  }

  dealii::VectorTools::interpolate_boundary_values(
      *(this->dof_handler), dirichlet_boundary_condition,
      this->constraint_matrix, component_mask);

  // Finally close the ConstraintMatrix.
  this->constraint_matrix.close();

  // Free anode_dirichlet_bc if necessary.
  if (anode_dirichlet_bc != nullptr)
  {
    delete anode_dirichlet_bc;
    anode_dirichlet_bc = nullptr;
  }

  // Create sparsity pattern
  unsigned int const max_couplings =
      this->dof_handler->max_couplings_between_dofs();
  this->sparsity_pattern.reinit(this->n_dofs, this->n_dofs, max_couplings);
  dealii::DoFTools::make_sparsity_pattern(
      *(this->dof_handler), this->sparsity_pattern, this->constraint_matrix);
  this->sparsity_pattern.compress();

  // Initialize matrices and vectors
  this->system_matrix.reinit(this->sparsity_pattern);
  this->mass_matrix.reinit(this->sparsity_pattern);
  this->system_rhs.reinit(this->n_dofs);

  assemble_system(parameters);
}

template <int dim>
void ElectrochemicalPhysics<dim>::assemble_system(
    std::shared_ptr<PhysicsParameters<dim> const> parameters)
{
  std::shared_ptr<
      ElectrochemicalPhysicsParameters<dim> const> electrochemical_parameters =
      std::dynamic_pointer_cast<ElectrochemicalPhysicsParameters<dim> const>(
          parameters);
  BOOST_ASSERT_MSG(electrochemical_parameters == nullptr,
                   "Problem during dowcasting the pointer");

  dealii::DoFHandler<dim> const &dof_handler = *(this->dof_handler);
  dealii::FiniteElement<dim> const &fe       = dof_handler.get_fe();

  dealii::FEValuesExtractors::Scalar const solid_potential(
      this->solid_potential_component);
  dealii::FEValuesExtractors::Scalar const liquid_potential(
      this->liquid_potential_component);
  dealii::QGauss<dim> quadrature_rule(fe.degree + 1);
  dealii::FEValues<dim> fe_values(
      fe, quadrature_rule, dealii::update_values | dealii::update_gradients |
                               dealii::update_JxW_values);

  unsigned int const dofs_per_cell = fe.dofs_per_cell;
  unsigned int const n_q_points    = quadrature_rule.size();
  unsigned int const time_step = electrochemical_parameters->time_step;
  dealii::Vector<double> cell_rhs(dofs_per_cell);
  dealii::FullMatrix<double> cell_system_matrix(dofs_per_cell, dofs_per_cell);
  dealii::FullMatrix<double> cell_mass_matrix(dofs_per_cell, dofs_per_cell);
  std::vector<double> solid_phase_diffusion_coefficient_values(n_q_points);
  std::vector<double> liquid_phase_diffusion_coefficient_values(n_q_points);
  std::vector<double> specific_capacitance_values(n_q_points);
  std::vector<double> faradaic_reaction_coefficient_values(n_q_points);
  std::vector<dealii::types::global_dof_index> local_dof_indices(dofs_per_cell);

  this->system_matrix = 0.0;
  this->mass_matrix   = 0.0;
  this->system_rhs    = 0.0;

  for (auto cell : dof_handler.active_cell_iterators())
  {
    cell_system_matrix = 0.0;
    cell_mass_matrix   = 0.0;
    cell_rhs           = 0.0;

    // clang-format off
      (this->mp_values)->get_values("specific_capacitance", cell, specific_capacitance_values);
      (this->mp_values)->get_values("solid_electrical_conductivity", cell, solid_phase_diffusion_coefficient_values);
      (this->mp_values)->get_values("liquid_electrical_conductivity", cell, liquid_phase_diffusion_coefficient_values);
      (this->mp_values)->get_values("faradaic_reaction_coefficient", cell, faradaic_reaction_coefficient_values);
    // clang-format on

    // The coefficients are zeros when the physics does not make sense.
    for (unsigned int q = 0; q < n_q_points; ++q)
      for (unsigned int i = 0; i < dofs_per_cell; ++i)
        for (unsigned int j = 0; j < dofs_per_cell; ++j)
        {
          // Mass matrix terms
          cell_mass_matrix(i, j) +=
              (specific_capacitance_values[q] *
                   (fe_values[solid_potential].value(i, q) *
                    fe_values[solid_potential].value(j, q)) -
               specific_capacitance_values[q] *
                   (fe_values[solid_potential].value(i, q) *
                    fe_values[liquid_potential].value(j, q)) -
               specific_capacitance_values[q] *
                   (fe_values[liquid_potential].value(i, q) *
                    fe_values[solid_potential].value(j, q)) +
               specific_capacitance_values[q] *
                   (fe_values[liquid_potential].value(i, q) *
                    fe_values[liquid_potential].value(j, q))) *
              fe_values.JxW(q);
          cell_system_matrix(i, j) +=
              cell_mass_matrix(i, j) +
              time_step *
                  // Stiffness matrix terms
                  (solid_phase_diffusion_coefficient_values[q] *
                       (fe_values[solid_potential].gradient(i, q) *
                        fe_values[solid_potential].gradient(j, q)) +
                   liquid_phase_diffusion_coefficient_values[q] *
                       (fe_values[liquid_potential].gradient(i, q) *
                        fe_values[liquid_potential].gradient(j, q)) +
                   faradaic_reaction_coefficient_values[q] *
                       (fe_values[solid_potential].value(i, q) *
                        fe_values[solid_potential].value(j, q)) -
                   faradaic_reaction_coefficient_values[q] *
                       (fe_values[liquid_potential].value(i, q) *
                        fe_values[solid_potential].value(j, q)) -
                   faradaic_reaction_coefficient_values[q] *
                       (fe_values[solid_potential].value(i, q) *
                        fe_values[liquid_potential].value(j, q)) +
                   faradaic_reaction_coefficient_values[q] *
                       (fe_values[liquid_potential].value(i, q) *
                        fe_values[liquid_potential].value(j, q))) *
                  fe_values.JxW(q);
        }

    cell->get_dof_indices(local_dof_indices);
    this->constraint_matrix.distribute_local_to_global(
        cell_system_matrix, cell_rhs, local_dof_indices, this->system_matrix,
        this->system_rhs);
    this->constraint_matrix.distribute_local_to_global(
        cell_mass_matrix, local_dof_indices, this->mass_matrix);
  }
}
}

#endif
