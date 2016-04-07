/* Copyright (c) 2016, the Cap authors.
 *
 * This file is subject to the Modified BSD License and may not be distributed
 * without copyright and license information. Please refer to the file LICENSE
 * for the text and further information on this license.
 */

#ifndef CAP_GEOMETRY_H
#define CAP_GEOMETRY_H

#include <cap/types.h>
#include <deal.II/base/types.h>
#include <deal.II/distributed/tria.h>
#include <boost/mpi.hpp>
#include <boost/property_tree/ptree.hpp>
#include <memory>
#include <unordered_map>

namespace cap
{
/**
 * This class allows access to the underlying
 * triangulation and the material id map. It also creates the triangulation by
 * reading the mesh from a ucd file.
 */
template <int dim>
class Geometry
{
public:
  /**
   * This contructor uses a mesh in ucd format. The database is used to get the
   * name of the mesh file and to create the materials map.
   */
  Geometry(std::shared_ptr<boost::property_tree::ptree const> database,
           boost::mpi::communicator mpi_communicator);

  /**
   * This contructor uses the dealii::distributed::Triangulation @p
   * triangulation. If the database is necessary to create the materials map.
   */
  Geometry(
      std::shared_ptr<boost::property_tree::ptree const> database,
      std::shared_ptr<dealii::distributed::Triangulation<dim>> triangulation);

  virtual ~Geometry() = default;

  inline std::shared_ptr<dealii::distributed::Triangulation<dim> const>
  get_triangulation() const
  {
    return this->triangulation;
  }

  inline boost::mpi::communicator get_mpi_communicator() const
  {
    return mpi_communicator;
  }

  inline std::shared_ptr<std::unordered_map<
      std::string, std::vector<dealii::types::material_id>> const>
  get_materials() const
  {
    return this->materials;
  }

private:
  /**
   * Helper function for the constructor.
   */
  void fill_materials_map(
      std::shared_ptr<boost::property_tree::ptree const> database);

  boost::mpi::communicator mpi_communicator;
  std::shared_ptr<dealii::distributed::Triangulation<dim>> triangulation;
  std::shared_ptr<std::unordered_map<
      std::string, std::vector<dealii::types::material_id>>> materials;
};
} // end namespace cap

#endif // CAP_GEOMETRY_H
