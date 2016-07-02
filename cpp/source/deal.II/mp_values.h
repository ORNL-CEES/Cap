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
      : database(d), geometry(nullptr)
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

  MPValues() = default;

  virtual ~MPValues() = default;

  virtual void get_values(std::string const &key,
                          active_cell_iterator const &cell,
                          std::vector<double> &values) const = 0;
};

template <int dim, int spacedim = dim>
class CompositeMat : public MPValues<dim, spacedim>
{
public:
  using active_cell_iterator =
      typename dealii::DoFHandler<dim, spacedim>::active_cell_iterator;

  CompositeMat() = default;

  void get_values(std::string const &key, active_cell_iterator const &cell,
                  std::vector<double> &values) const override;

protected:
  std::unordered_map<dealii::types::material_id,
                     std::shared_ptr<MPValues<dim, spacedim>>> _materials = {};
};

template <int dim, int spacedim = dim>
class CompositePro : public MPValues<dim, spacedim>
{
public:
  using active_cell_iterator =
      typename dealii::DoFHandler<dim, spacedim>::active_cell_iterator;

  CompositePro() = default;

  void get_values(std::string const &key, active_cell_iterator const &cell,
                  std::vector<double> &values) const override;

protected:
  std::unordered_map<std::string, std::shared_ptr<MPValues<dim, spacedim>>>
      _properties = {};
};

template <int dim, int spacedim = dim>
class UniformConstantMPValues : public MPValues<dim, spacedim>
{
public:
  using active_cell_iterator =
      typename dealii::DoFHandler<dim, spacedim>::active_cell_iterator;

  UniformConstantMPValues(double const &val);

  void get_values(std::string const &key, active_cell_iterator const &cell,
                  std::vector<double> &values) const override;

protected:
  // get_values(...) will assign _val to all elements in the vector values.
  double _val;
};

template <int dim, int spacedim = dim>
class SuperCapacitorMPValues : public CompositeMat<dim, spacedim>
{
public:
  SuperCapacitorMPValues(MPValuesParameters<dim, spacedim> const &params);
  // Use build(...) to create either an homogeneous
  // (SuperCapacitorMPValues) or an inhomogeneous model
  // (InhomogeneousSuperCapacitorMPValues)
  static std::unique_ptr<MPValues<dim, spacedim>>
  build(MPValuesParameters<dim, spacedim> const &params);
};

template <int dim, int spacedim = dim>
class InhomogeneousSuperCapacitorMPValues : public MPValues<dim, spacedim>
{
public:
  using active_cell_iterator =
      typename dealii::DoFHandler<dim, spacedim>::active_cell_iterator;

  InhomogeneousSuperCapacitorMPValues(
      MPValuesParameters<dim, spacedim> const &params);

  void get_values(std::string const &key, active_cell_iterator const &cell,
                  std::vector<double> &values) const override;

protected:
  std::map<dealii::CellId, std::shared_ptr<MPValues<dim, spacedim>>> _map = {};
};

template <int dim, int spacedim = dim>
class PorousElectrodeMPValues : public CompositePro<dim, spacedim>
{
public:
  PorousElectrodeMPValues(std::string const &material_name,
                          MPValuesParameters<dim, spacedim> const &params);
};
template <int dim, int spacedim = dim>
class MetalFoilMPValues : public CompositePro<dim, spacedim>
{
public:
  MetalFoilMPValues(std::string const &material_name,
                    MPValuesParameters<dim, spacedim> const &params);
};

} // end namespace cap

#endif // CAP_MP_VALUES_H
