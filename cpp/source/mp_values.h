/* Copyright (c) 2016, the Cap authors.
 *
 * This file is subject to the Modified BSD License and may not be distributed
 * without copyright and license information. Please refer to the file LICENSE
 * for the text and further information on this license.
 */

#ifndef CAP_MP_VALUES_H
#define CAP_MP_VALUES_H

#include <cap/geometry.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>
#include <boost/property_tree/ptree.hpp>

namespace cap
{

//////////////////////// MP VALUES PARAMETERS ////////////////////////////
template <int dim, int spacedim = dim>
class MPValuesParameters
{
public:
  MPValuesParameters(std::shared_ptr<boost::property_tree::ptree const> d)
      : database(d)
  {
  }
  virtual ~MPValuesParameters() = default;
  // keep public for now
  std::shared_ptr<boost::property_tree::ptree const> database;
  std::shared_ptr<Geometry<dim> const> geometry;
};

//////////////////////// MP VALUES ////////////////////////////
template <int dim, int spacedim = dim>
class MPValues
{
public:
  typedef typename dealii::DoFHandler<dim, spacedim>::active_cell_iterator
      active_cell_iterator;

  MPValues(MPValuesParameters<dim, spacedim> const &params);

  virtual ~MPValues() = default;

  virtual void get_values(std::string const &key,
                          active_cell_iterator const &cell,
                          std::vector<double> &values) const;

  virtual void
  get_values(std::string const &key, active_cell_iterator const &cell,
             std::vector<dealii::Tensor<1, spacedim>> &values) const;

protected:
  std::unordered_map<dealii::types::material_id, std::shared_ptr<MPValues<dim>>>
      materials;
};

//////////////////////// NEW STUFF ////////////////////////////
template <int dim, int spacedim = dim>
class NewStuffMPValues : public MPValues<dim, spacedim>
{
public:
  typedef typename dealii::DoFHandler<dim, spacedim>::active_cell_iterator
      active_cell_iterator;

  NewStuffMPValues(MPValuesParameters<dim, spacedim> const &parameters);

  // Needed to fix hidding of get_values.
  using MPValues<dim, spacedim>::get_values;

  void get_values(std::string const &key, active_cell_iterator const &cell,
                  std::vector<double> &values) const override;

protected:
  std::unordered_map<std::string,
                     std::function<void(active_cell_iterator const &,
                                        std::vector<double> &)>> properties;
};

template <int dim, int spacedim = dim>
class PorousElectrodeMPValues : public NewStuffMPValues<dim, spacedim>
{
public:
  typedef typename dealii::DoFHandler<dim, spacedim>::active_cell_iterator
      active_cell_iterator;
  PorousElectrodeMPValues(MPValuesParameters<dim, spacedim> const &parameters);
};

template <int dim, int spacedim = dim>
class MetalFoilMPValues : public NewStuffMPValues<dim, spacedim>
{
public:
  typedef typename dealii::DoFHandler<dim, spacedim>::active_cell_iterator
      active_cell_iterator;
  MetalFoilMPValues(MPValuesParameters<dim, spacedim> const &parameters);
};

template <int dim, int spacedim = dim>
std::shared_ptr<MPValues<dim>>
buildMaterial(std::string const &material_name,
              std::shared_ptr<boost::property_tree::ptree const> database)
{
  std::shared_ptr<boost::property_tree::ptree> material_database =
      std::make_shared<boost::property_tree::ptree>(
          database->get_child(material_name));
  std::string const type = material_database->get<std::string>("type");
  if (type.compare("porous_electrode") == 0)
  {
    std::shared_ptr<boost::property_tree::ptree> dummy_database =
        std::make_shared<boost::property_tree::ptree>(*database);
    dummy_database->put("ugly_hack", material_name);
    return std::make_shared<PorousElectrodeMPValues<dim>>(
        MPValuesParameters<dim>(dummy_database));
  }
  else if (type.compare("permeable_membrane") == 0)
  {
    std::shared_ptr<boost::property_tree::ptree> dummy_database =
        std::make_shared<boost::property_tree::ptree>(*database);
    dummy_database->put("ugly_hack", material_name);
    std::string const matrix_phase =
        dummy_database->get<std::string>(material_name + "." + "matrix_phase");
    dummy_database->put(matrix_phase + "." + "differential_capacitance", 0.0);
    dummy_database->put(matrix_phase + "." + "exchange_current_density", 0.0);
    dummy_database->put(matrix_phase + "." + "electrical_resistivity",
                        std::numeric_limits<double>::max());
    return std::make_shared<PorousElectrodeMPValues<dim>>(
        MPValuesParameters<dim>(dummy_database));
  }
  else if (type.compare("current_collector") == 0)
  {
    std::shared_ptr<boost::property_tree::ptree> dummy_database =
        std::make_shared<boost::property_tree::ptree>(*database);
    dummy_database->put("ugly_hack", material_name);
    return std::make_shared<MetalFoilMPValues<dim>>(
        MPValuesParameters<dim>(dummy_database));
  }
  else
  {
    throw std::runtime_error("Invalid material type " + type);
  }
}

} // end namespace cap

#endif // CAP_MP_VALUES_H
