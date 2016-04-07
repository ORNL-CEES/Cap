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
   * triangulation. The database is necessary to create the materials map.
   */
  Geometry(
      std::shared_ptr<boost::property_tree::ptree const> database,
      std::shared_ptr<dealii::distributed::Triangulation<dim>> triangulation);

  virtual ~Geometry() = default;

  dealii::types::boundary_id get_anode_boundary_id() const
  {
    return _anode_boundary_id;
  }

  dealii::types::boundary_id get_cathode_boundary_id() const
  {
    return _cathode_boundary_id;
  }

  std::shared_ptr<dealii::distributed::Triangulation<dim> const>
  get_triangulation() const
  {
    return _triangulation;
  }

  boost::mpi::communicator get_mpi_communicator() const
  {
    return _communicator;
  }

  std::shared_ptr<
      std::unordered_map<std::string, std::vector<dealii::types::material_id>>>
  get_materials() const
  {
    return _materials;
  }

private:
  /**
   * Helper function for the constructor, when the mesh is loaded from a mesh.
   */
  void fill_materials_map(
      std::shared_ptr<boost::property_tree::ptree const> database);

  /**
   * Helper function for the constructor, when the mesh is generated from a
   * database.
   */
  void fill_materials_map();

  /**
   * Create a mesh from a property tree.
   */
  void mesh_generator(boost::property_tree::ptree const &database);

  /**
   * Set the boundary IDs on the cathode and the anode.
   */
  void set_boundary_ids(double const collector_top,
                        double const collector_bottom);

  boost::mpi::communicator _communicator;
  dealii::types::boundary_id _anode_boundary_id;
  dealii::types::boundary_id _cathode_boundary_id;
  std::shared_ptr<dealii::distributed::Triangulation<dim>> _triangulation;
  std::shared_ptr<std::unordered_map<
      std::string, std::vector<dealii::types::material_id>>> _materials;
};
} // end namespace cap

#endif // CAP_GEOMETRY_H
