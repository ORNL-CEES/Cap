/* Copyright (c) 2016, the Cap authors.
 *
 * This file is subject to the Modified BSD License and may not be distributed
 * without copyright and license information. Please refer to the file LICENSE
 * for the text and further information on this license.
 */

#ifndef CAP_MP_VALUES_H
#define CAP_MP_VALUES_H

#include <cap/geometry.h>
#include <deal.II/grid/cell_id.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/function.h>
#include <boost/property_tree/ptree.hpp>

namespace cap
{

//////////////////////// MP VALUES PARAMETERS ////////////////////////////
template <int dim>
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
  // cannot be const because get_triangulation(...) is not const...
  std::shared_ptr<Geometry<dim>> geometry;
};

//////////////////////// MP VALUES ////////////////////////////
template <int dim>
class MPValues
{
public:
  MPValues() = default;

  virtual ~MPValues() = default;

  virtual void get_values(std::string const &key,
                          dealii::FEValues<dim> const &fe_values,
                          std::vector<double> &values) const = 0;
};

template <int dim>
class CompositeMat : public MPValues<dim>
{
public:
  CompositeMat() = default;

  void get_values(std::string const &key,
                  dealii::FEValues<dim> const &fe_values,
                  std::vector<double> &values) const override;

protected:
  std::unordered_map<dealii::types::material_id, std::shared_ptr<MPValues<dim>>>
      _materials = {};
};

template <int dim>
class CompositePro : public MPValues<dim>
{
public:
  CompositePro() = default;

  void get_values(std::string const &key,
                  dealii::FEValues<dim> const &fe_values,
                  std::vector<double> &values) const override;

protected:
  std::unordered_map<std::string, std::shared_ptr<MPValues<dim>>> _properties =
      {};
};

template <int dim>
class UniformConstantMPValues : public MPValues<dim>
{
public:
  UniformConstantMPValues(double const &val);

  void get_values(std::string const &key,
                  dealii::FEValues<dim> const &fe_values,
                  std::vector<double> &values) const override;

protected:
  // get_values(...) will assign _val to all elements in the vector values.
  double _val;
};

template <int dim>
class FunctionSpaceMPValues : public MPValues<dim>
{
public:
  FunctionSpaceMPValues(std::shared_ptr<dealii::Function<dim>> const &function);

  /*!
    \note
    Here \p fe_values must be constructed with the flag \c
    dealii::update_quadrature_points.
  */
  void get_values(std::string const &key,
                  dealii::FEValues<dim> const &fe_values,
                  std::vector<double> &values) const override;

private:
  std::shared_ptr<dealii::Function<dim>> const _function;
};

template <int dim>
class SuperCapacitorMPValuesFactory
{
public:
  static std::unique_ptr<MPValues<dim>>
  build(MPValuesParameters<dim> const &params);
};

template <int dim>
class SuperCapacitorMPValues : public CompositeMat<dim>
{
public:
  SuperCapacitorMPValues(MPValuesParameters<dim> const &params);
};

template <int dim>
class InhomogeneousSuperCapacitorMPValues : public MPValues<dim>
{
public:
  InhomogeneousSuperCapacitorMPValues(MPValuesParameters<dim> const &params);

  void get_values(std::string const &key,
                  dealii::FEValues<dim> const &fe_values,
                  std::vector<double> &values) const override;

protected:
  std::map<dealii::CellId, std::shared_ptr<MPValues<dim>>> _map = {};
};

template <int dim>
class PorousElectrodeMPValues : public CompositePro<dim>
{
public:
  PorousElectrodeMPValues(std::string const &material_name,
                          MPValuesParameters<dim> const &params);
};

template <int dim>
class MetalFoilMPValues : public CompositePro<dim>
{
public:
  MetalFoilMPValues(std::string const &material_name,
                    MPValuesParameters<dim> const &params);
};

} // end namespace cap

#endif // CAP_MP_VALUES_H
