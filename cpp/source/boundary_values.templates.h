/* Copyright (c) 2016, the Cap authors.
 *
 * This file is subject to the Modified BSD License and may not be distributed
 * without copyright and license information. Please refer to the file LICENSE
 * for the text and further information on this license. 
 */

#include <cap/boundary_values.h>
#include <stdexcept>

namespace cap
{

template <int dim, int spacedim>
void BoundaryValues<dim, spacedim>::get_values(
    std::string const &key, active_cell_iterator const &cell,
    unsigned int const face, std::vector<double> &values) const

{
  std::ignore = key;
  std::ignore = cell;
  std::ignore = face;
  std::ignore = values;
  throw std::runtime_error("Should have been overloaded...");
}

template <int dim, int spacedim>
void BoundaryValues<dim, spacedim>::get_values(
    std::string const &key, active_cell_iterator const &cell,
    unsigned int const face,
    std::vector<dealii::Tensor<1, spacedim>> &values) const

{
  std::ignore = key;
  std::ignore = cell;
  std::ignore = face;
  std::ignore = values;
  throw std::runtime_error("Should have been overloaded...");
}

//////////////////////// SUPERCAPACITORS ////////////////////////////
template <int dim, int spacedim>
SuperCapacitorBoundaryValues<dim, spacedim>::SuperCapacitorBoundaryValues(
    BoundaryValuesParameters<dim, spacedim> const &parameters)
    : BoundaryValues<dim, spacedim>(parameters)
{
  std::shared_ptr<boost::property_tree::ptree> database = parameters.database;

  //    SuperCapacitorBoundaryValuesParameters<dim, spacedim> const *
  //    super_capacitor_parameters =
  //        dynamic_cast<SuperCapacitorBoundaryValuesParameters<dim, spacedim>
  //        const *>(&parameters);

  // clang-format off
  this->separator_material_id         = database->get<dealii::types::material_id>("separator_material_id");
  this->anode_electrode_material_id   = database->get<dealii::types::material_id>("anode_electrode_material_id");
  this->anode_collector_material_id   = database->get<dealii::types::material_id>("anode_collector_material_id");
  this->cathode_electrode_material_id = database->get<dealii::types::material_id>("cathode_electrode_material_id");
  this->cathode_collector_material_id = database->get<dealii::types::material_id>("cathode_collector_material_id");
  // clang-format on

  // clang-format off
  this->cathode_boundary_id = database->get<dealii::types::boundary_id>("cathode_boundary_id");
  this->anode_boundary_id   = database->get<dealii::types::boundary_id>("anode_boundary_id");
  this->upper_boundary_id   = database->get<dealii::types::boundary_id>("upper_boundary_id");
  this->lower_boundary_id   = database->get<dealii::types::boundary_id>("lower_boundary_id");
  this->other_boundary_id   = database->get<dealii::types::boundary_id>("other_boundary_id");
  //clang-format on

  // clang-format off
  this->upper_ambient_temperature       = database->get<double>("ambient_temperature");
  this->lower_ambient_temperature       = database->get<double>("ambient_temperature");
  this->upper_heat_transfer_coefficient = database->get<double>("heat_transfer_coefficient");
  this->lower_heat_transfer_coefficient = 0.0;
  // clang-format on
}

template <int dim, int spacedim>
void SuperCapacitorBoundaryValues<dim, spacedim>::get_values(
    std::string const &key, active_cell_iterator const &cell,
    unsigned int const face, std::vector<double> &values) const
{
  if (key.compare("ambient_temperature") == 0)
  {
    if ((cell->face(face)->boundary_id() == this->anode_boundary_id) ||
        (cell->face(face)->boundary_id() == this->cathode_boundary_id) ||
        (cell->face(face)->boundary_id() == this->other_boundary_id))
    {
      std::fill(values.begin(), values.end(), 0.0);
    }
    else if (cell->face(face)->boundary_id() == this->upper_boundary_id)
    {
      std::fill(values.begin(), values.end(), this->upper_ambient_temperature);
    }
    else if (cell->face(face)->boundary_id() == this->lower_boundary_id)
    {
      std::fill(values.begin(), values.end(), this->lower_ambient_temperature);
    }
    else
    {
      throw std::runtime_error("Invalid boundary id");
    } // end if boundary id
  }
  else if (key.compare("heat_transfer_coefficient") == 0)
  {
    if ((cell->face(face)->boundary_id() == this->anode_boundary_id) ||
        (cell->face(face)->boundary_id() == this->cathode_boundary_id) ||
        (cell->face(face)->boundary_id() == this->other_boundary_id))
    {
      std::fill(values.begin(), values.end(), 0.0);
    }
    else if (cell->face(face)->boundary_id() == this->upper_boundary_id)
    {
      std::fill(values.begin(), values.end(),
                this->upper_heat_transfer_coefficient);
    }
    else if (cell->face(face)->boundary_id() == this->lower_boundary_id)
    {
      std::fill(values.begin(), values.end(),
                this->lower_heat_transfer_coefficient);
    }
    else
    {
      throw std::runtime_error("Invalid boundary id");
    } // end if boundary id
  }
  else
  {
    throw std::runtime_error("Invalid key");
  } // end if key
}

} // end namespace cap
