/* Copyright (c) 2016, the Cap authors.
 *
 * This file is subject to the Modified BSD License and may not be distributed
 * without copyright and license information. Please refer to the file LICENSE
 * for the text and further information on this license.
 */

#include <cap/post_processor.h>
#include <cap/utils.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_values.h>
#include <algorithm>
#include <functional>
#include <limits>

namespace cap
{

//////////////////////// POSTPROCESSOR ////////////////////////////
template <int dim>
Postprocessor<dim>::Postprocessor(
    std::shared_ptr<PostprocessorParameters<dim> const> parameters,
    boost::mpi::communicator mpi_communicator)
    : _communicator(mpi_communicator), dof_handler(parameters->dof_handler),
      solution(parameters->solution), mp_values(parameters->mp_values)
{
  BOOST_ASSERT_MSG(dof_handler != nullptr, "Invalid DoFHandler");
  BOOST_ASSERT_MSG(solution != nullptr, "Invalid solution vector.");
  BOOST_ASSERT_MSG(mp_values != nullptr, "Invalid MPValues.");
}

template <int dim>
dealii::Vector<double> const &
Postprocessor<dim>::get(std::string const &key) const
{
  std::unordered_map<std::string, dealii::Vector<double>>::const_iterator it =
      this->vectors.find(key);
  AssertThrow(it != this->vectors.end(), dealii::StandardExceptions::ExcMessage(
                                             "Key " + key + " doesn't exist"));
  return it->second;
}

template <int dim>
void Postprocessor<dim>::get(std::string const &key, double &value) const
{
  std::unordered_map<std::string, double>::const_iterator it =
      this->values.find(key);
  AssertThrow(it != this->values.end(), dealii::StandardExceptions::ExcMessage(
                                            "Key " + key + " doesn't exist"));
  value = it->second;
}

template <int dim>
std::vector<std::string> Postprocessor<dim>::get_vector_keys() const
{
  std::vector<std::string> keys;
  std::unordered_map<std::string, dealii::Vector<double>>::const_iterator it =
      this->vectors.begin();
  std::unordered_map<std::string, dealii::Vector<double>>::const_iterator
      end_it = this->vectors.end();
  for (; it != end_it; ++it)
    keys.push_back(it->first);
  return keys;
}

//////////////////////// SUPERCAPACITOR POSTPROCESSOR PARAMETERS //////
template <int dim>
SuperCapacitorPostprocessorParameters<dim>::
    SuperCapacitorPostprocessorParameters(
        std::shared_ptr<boost::property_tree::ptree const> d,
        std::shared_ptr<dealii::DoFHandler<dim>> dof_handler)
    : PostprocessorParameters<dim>(d, dof_handler)
{
}

//////////////////////// SUPERCAPACITOR POSTPROCESSOR /////////////////
template <int dim>
SuperCapacitorPostprocessor<dim>::SuperCapacitorPostprocessor(
    std::shared_ptr<PostprocessorParameters<dim> const> parameters,
    boost::mpi::communicator mpi_communicator)
    : Postprocessor<dim>(parameters, mpi_communicator)
{
  dealii::DoFHandler<dim> const &dof_handler = *(this->dof_handler);
  this->values["voltage"] = 0.0;
  this->values["current"] = 0.0;
  this->values["joule_heating"] = 0.0;
  this->values["surface_area"] = 0.0;
  this->values["volume"] = 0.0;
  this->values["mass"] = 0.0;
  this->values["anode_electrode_interfacial_surface_area"] = 0.0;
  this->values["anode_electrode_mass_of_active_material"] = 0.0;
  this->values["cathode_electrode_interfacial_surface_area"] = 0.0;
  this->values["cathode_electrode_mass_of_active_material"] = 0.0;

  std::shared_ptr<boost::property_tree::ptree const> database =
      parameters->database;

  this->debug_material_properties = cap::to_vector<std::string>(
      database->get("debug.material_properties", ""));
  this->debug_solution_fields =
      cap::to_vector<std::string>(database->get("debug.solution_fields", ""));
  this->debug_solution_fluxes =
      cap::to_vector<std::string>(database->get("debug.solution_fluxes", ""));
  this->debug_boundary_ids = database->get("debug.boundary_ids", false);
  this->debug_material_ids = database->get("debug.material_ids", false);

  // n_active_cells is the total number of locally owned cells. This includes
  // the ghost cells and the artificial cells. They are filtered out by DataOut.
  unsigned int const n_active_cells =
      dof_handler.get_triangulation().n_active_cells();
  if (this->debug_material_ids)
    this->vectors["material_id"] = dealii::Vector<double>(n_active_cells);
  if (this->debug_boundary_ids)
    throw dealii::StandardExceptions::ExcMessage("not implemented yet");
  for (std::vector<std::string>::const_iterator it =
           this->debug_material_properties.begin();
       it != this->debug_material_properties.end(); ++it)
    this->vectors[*it] = dealii::Vector<double>(n_active_cells);
  for (std::vector<std::string>::const_iterator it =
           this->debug_solution_fields.begin();
       it != this->debug_solution_fields.end(); ++it)
    this->vectors[*it] = dealii::Vector<double>(n_active_cells);
  for (std::vector<std::string>::const_iterator it =
           this->debug_solution_fluxes.begin();
       it != this->debug_solution_fluxes.end(); ++it)
    for (int d = 0; d < dim; ++d)
      this->vectors[(*it) + "_" + std::to_string(d)] =
          dealii::Vector<double>(n_active_cells);
}

template <int dim>
void SuperCapacitorPostprocessor<dim>::reset(
    std::shared_ptr<PostprocessorParameters<dim> const> parameters)
{
  dealii::DoFHandler<dim> const &dof_handler = *(this->dof_handler);
  dealii::Trilinos::MPI::BlockVector const &solution = *(this->solution);

  std::for_each(this->values.begin(), this->values.end(),
                [](std::unordered_map<std::string, double>::value_type &p)
                {
                  p.second = 0.0;
                });

  std::shared_ptr<boost::property_tree::ptree const> database =
      parameters->database;

  // clang-format off
  dealii::types::material_id const anode_electrode_material_id   = database->get<dealii::types::material_id>("geometry.anode_electrode_material_id");
  dealii::types::material_id const cathode_electrode_material_id = database->get<dealii::types::material_id>("geometry.cathode_electrode_material_id");
  dealii::types::boundary_id const cathode_boundary_id           = database->get<dealii::types::boundary_id>("boundary_values.cathode_boundary_id");
  dealii::FEValuesExtractors::Scalar const solid_potential (database->get<unsigned int>("solid_potential_component"));
  dealii::FEValuesExtractors::Scalar const liquid_potential(database->get<unsigned int>("liquid_potential_component"));
  // clang-format on

  dealii::Vector<double> joule_heating(
      dof_handler.get_triangulation().n_active_cells());

  dealii::FiniteElement<dim> const &fe =
      dof_handler.get_fe(); // TODO: don't want to use directly fe because we
                            // might create postprocessor that will only know
                            // about dof_handler
  dealii::QGauss<dim> quadrature_rule(fe.degree + 1);
  dealii::QGauss<dim - 1> face_quadrature_rule(fe.degree + 1);
  dealii::FEValues<dim> fe_values(
      fe, quadrature_rule, dealii::update_values | dealii::update_gradients |
                               dealii::update_JxW_values);
  dealii::FEFaceValues<dim> fe_face_values(
      fe, face_quadrature_rule,
      dealii::update_values | dealii::update_gradients |
          dealii::update_JxW_values | dealii::update_normal_vectors);
  //    unsigned int const dofs_per_cell   = fe.dofs_per_cell;
  unsigned int const n_q_points = quadrature_rule.size();
  unsigned int const n_face_q_points = face_quadrature_rule.size();
  std::vector<double> solid_electrical_conductivity_values(n_q_points);
  std::vector<double> liquid_electrical_conductivity_values(n_q_points);
  std::vector<double> density_values(n_q_points);
  std::vector<double> density_of_active_material_values(n_q_points);
  std::vector<double> specific_surface_area_values(n_q_points);
  std::vector<dealii::Tensor<1, dim>> solid_potential_gradients(n_q_points);
  std::vector<dealii::Tensor<1, dim>> liquid_potential_gradients(n_q_points);
  std::vector<double> solid_potential_values(n_q_points);
  std::vector<double> liquid_potential_values(n_q_points);
  double anode_electrode_potential = 0.0;
  double cathode_electrode_potential = 0.0;
  double anode_electrode_volume = 0.0;
  double cathode_electrode_volume = 0.0;

  std::vector<double> face_solid_electrical_conductivity_values(
      n_face_q_points);
  std::vector<double> face_solid_potential_values(n_face_q_points);
  std::vector<dealii::Tensor<1, dim>> face_solid_potential_gradients(
      n_face_q_points);
  std::vector<dealii::Tensor<1, dim>> normal_vectors(n_face_q_points);

  dealii::IndexSet locally_relevant_dofs;
  dealii::DoFTools::extract_locally_relevant_dofs(dof_handler,
                                                  locally_relevant_dofs);
  std::vector<dealii::IndexSet> index_sets(1, locally_relevant_dofs);
  dealii::Trilinos::MPI::BlockVector relevant_solution(index_sets);
  relevant_solution = solution;

  for (auto cell : dof_handler.active_cell_iterators())
  {
    if (cell->is_locally_owned())
    {
      fe_values.reinit(cell);
      this->mp_values->get_values("solid_electrical_conductivity", cell,
                                  solid_electrical_conductivity_values);
      this->mp_values->get_values("liquid_electrical_conductivity", cell,
                                  liquid_electrical_conductivity_values);
      this->mp_values->get_values("density", cell, density_values);
      this->mp_values->get_values("density_of_active_material", cell,
                                  density_of_active_material_values);
      this->mp_values->get_values("specific_surface_area", cell,
                                  specific_surface_area_values);
      if (*std::max_element(solid_electrical_conductivity_values.begin(),
                            solid_electrical_conductivity_values.end()) >
          1e-300)
      {
        fe_values[solid_potential].get_function_gradients(
            relevant_solution, solid_potential_gradients);
        fe_values[solid_potential].get_function_values(relevant_solution,
                                                       solid_potential_values);
      }
      if (*std::max_element(liquid_electrical_conductivity_values.begin(),
                            liquid_electrical_conductivity_values.end()) >
          1e-300)
      {
        fe_values[liquid_potential].get_function_gradients(
            relevant_solution, liquid_potential_gradients);
        fe_values[liquid_potential].get_function_values(
            relevant_solution, liquid_potential_values);
      }
      for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
      {
        this->values["joule_heating"] +=
            (solid_electrical_conductivity_values[q_point] *
                 solid_potential_gradients[q_point] *
                 solid_potential_gradients[q_point] +
             liquid_electrical_conductivity_values[q_point] *
                 liquid_potential_gradients[q_point] *
                 liquid_potential_gradients[q_point]) *
            fe_values.JxW(q_point);
        this->values["volume"] += fe_values.JxW(q_point);
        this->values["mass"] +=
            density_values[q_point] * fe_values.JxW(q_point);
        if (cell->material_id() == anode_electrode_material_id)
        {
          anode_electrode_potential += (solid_potential_values[q_point] -
                                        liquid_potential_values[q_point]) *
                                       fe_values.JxW(q_point);
          anode_electrode_volume += fe_values.JxW(q_point);
          this->values["anode_electrode_interfacial_surface_area"] +=
              specific_surface_area_values[q_point] * fe_values.JxW(q_point);
          this->values["anode_electrode_mass_of_active_material"] +=
              density_of_active_material_values[q_point] *
              fe_values.JxW(q_point);
        }
        else if (cell->material_id() == cathode_electrode_material_id)
        {
          cathode_electrode_potential += (solid_potential_values[q_point] -
                                          liquid_potential_values[q_point]) *
                                         fe_values.JxW(q_point);
          cathode_electrode_volume += fe_values.JxW(q_point);
          this->values["cathode_electrode_interfacial_surface_area"] +=
              specific_surface_area_values[q_point] * fe_values.JxW(q_point);
          this->values["cathode_electrode_mass_of_active_material"] +=
              density_of_active_material_values[q_point] *
              fe_values.JxW(q_point);
        }
        else
        {
          // do nothing
        }
      } // end for quadrature point
      if (this->debug_material_ids)
        this->vectors["material_id"][cell->active_cell_index()] =
            static_cast<double>(cell->material_id());
      for (std::vector<std::string>::const_iterator it =
               this->debug_material_properties.begin();
           it != this->debug_material_properties.end(); ++it)
      {
        std::vector<double> values(n_q_points);
        this->mp_values->get_values(*it, cell, values);
        double cell_averaged_value = 0.0;
        for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
        {
          cell_averaged_value += values[q_point] * fe_values.JxW(q_point);
        }
        cell_averaged_value /= cell->measure();
        this->vectors[*it][cell->active_cell_index()] = cell_averaged_value;
      }
      for (std::vector<std::string>::const_iterator it =
               this->debug_solution_fields.begin();
           it != this->debug_solution_fields.end(); ++it)
      {
        std::vector<double> values(n_q_points);
        if (it->compare("solid_potential") == 0)
        {
          values = solid_potential_values;
        }
        else if (it->compare("liquid_potential") == 0)
        {
          values = liquid_potential_values;
        }
        else if (it->compare("overpotential") == 0)
        {
          std::transform(solid_potential_values.begin(),
                         solid_potential_values.end(),
                         liquid_potential_values.begin(), values.begin(),
                         std::minus<double>());
        }
        else if (it->compare("joule_heating") == 0)
        {
          for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
          {
            values[q_point] =
                liquid_electrical_conductivity_values[q_point] *
                    liquid_potential_gradients[q_point].norm_square() +
                solid_electrical_conductivity_values[q_point] *
                    solid_potential_gradients[q_point].norm_square();
          }
        }
        else
        {
          throw dealii::StandardExceptions::ExcMessage(
              "Solution field '" + (*it) + "' is not recognized");
        }
        double cell_averaged_value = 0.0;
        for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
        {
          cell_averaged_value += values[q_point] * fe_values.JxW(q_point);
        }
        cell_averaged_value /= cell->measure();
        this->vectors[*it][cell->active_cell_index()] = cell_averaged_value;
      }
      for (std::vector<std::string>::const_iterator it =
               this->debug_solution_fluxes.begin();
           it != this->debug_solution_fluxes.end(); ++it)
      {
        std::vector<dealii::Tensor<1, dim>> values(n_q_points);
        if (it->compare("solid_current_density") == 0)
        {
          std::transform(solid_electrical_conductivity_values.begin(),
                         solid_electrical_conductivity_values.end(),
                         solid_potential_gradients.begin(), values.begin(),
                         [](double const x, dealii::Tensor<1, dim> const &y)
                         {
                           return x * y;
                         });
        }
        else if (it->compare("liquid_current_density") == 0)
        {
          std::transform(liquid_electrical_conductivity_values.begin(),
                         liquid_electrical_conductivity_values.end(),
                         liquid_potential_gradients.begin(), values.begin(),
                         [](double const x, dealii::Tensor<1, dim> const &y)
                         {
                           return x * y;
                         });
        }
        else
        {
          throw dealii::StandardExceptions::ExcMessage(
              "Solution flux '" + (*it) + "' is not recognized");
        }
        dealii::Tensor<1, dim> cell_averaged_value;
        cell_averaged_value = 0.0;
        for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
        {
          cell_averaged_value += values[q_point] * fe_values.JxW(q_point);
        }
        cell_averaged_value /= cell->measure();
        for (int d = 0; d < dim; ++d)
          this->vectors[(*it) + "_" +
                        std::to_string(d)][cell->active_cell_index()] =
              cell_averaged_value[d];
      }

      if (cell->at_boundary())
      {
        for (unsigned int face = 0;
             face < dealii::GeometryInfo<dim>::faces_per_cell; ++face)
        {
          if (cell->face(face)->at_boundary())
          {
            if (cell->face(face)->boundary_id() == cathode_boundary_id)
            {
              fe_face_values.reinit(cell, face);
              this->mp_values->get_values(
                  "solid_electrical_conductivity", cell,
                  face_solid_electrical_conductivity_values); // TODO: should
                                                              // take
                                                              // face as an
                                                              // argument...
              fe_face_values[solid_potential].get_function_gradients(
                  relevant_solution, face_solid_potential_gradients);
              fe_face_values[solid_potential].get_function_values(
                  relevant_solution, face_solid_potential_values);
              normal_vectors = fe_face_values.get_all_normal_vectors();
              for (unsigned int face_q_point = 0;
                   face_q_point < n_face_q_points; ++face_q_point)
              {
                this->values["current"] +=
                    (face_solid_electrical_conductivity_values[face_q_point] *
                     face_solid_potential_gradients[face_q_point] *
                     normal_vectors[face_q_point]) *
                    fe_face_values.JxW(face_q_point);
                this->values["voltage"] +=
                    face_solid_potential_values[face_q_point] *
                    fe_face_values.JxW(face_q_point);
                this->values["surface_area"] +=
                    fe_face_values.JxW(face_q_point);
              } // end for face quadrature point
            }   // end if cathode
          }     // end if face at boundary
        }       // end for face
      }         // end if cell at boundary
    }
  } // end for cell
  // AllReduce to get the scalar quantities
  this->values["current"] =
      dealii::Utilities::MPI::sum(this->values["current"], this->_communicator);
  this->values["surface_area"] = dealii::Utilities::MPI::sum(
      this->values["surface_area"], this->_communicator);
  this->values["voltage"] =
      dealii::Utilities::MPI::sum(this->values["voltage"], this->_communicator);
  anode_electrode_potential = dealii::Utilities::MPI::sum(
      anode_electrode_potential, this->_communicator);
  anode_electrode_volume =
      dealii::Utilities::MPI::sum(anode_electrode_volume, this->_communicator);
  cathode_electrode_potential = dealii::Utilities::MPI::sum(
      cathode_electrode_potential, this->_communicator);
  cathode_electrode_volume = dealii::Utilities::MPI::sum(
      cathode_electrode_volume, this->_communicator);

  this->values["voltage"] /= this->values["surface_area"];
  anode_electrode_potential /= anode_electrode_volume;
  cathode_electrode_potential /= cathode_electrode_volume;
  this->values["anode_potential"] = anode_electrode_potential;
  this->values["cathode_potential"] = cathode_electrode_potential;
}

} // end namespace cap
