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
#include <set>

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
  Geometry(std::shared_ptr<boost::property_tree::ptree> database,
           boost::mpi::communicator mpi_communicator);

  /**
   * This contructor uses the dealii::distributed::Triangulation @p
   * triangulation. The database is necessary to create the materials map.
   */
  Geometry(
      std::shared_ptr<dealii::distributed::Triangulation<dim>> triangulation,
      std::shared_ptr<std::unordered_map<
          std::string, std::set<dealii::types::material_id>>> materials,
      std::shared_ptr<std::unordered_map<
          std::string, std::set<dealii::types::boundary_id>>> boundaries);

  virtual ~Geometry() = default;

  /**
   * Read the weights of the cells and repartion the Triangulation.
   */
  void repartition();

  std::shared_ptr<dealii::distributed::Triangulation<dim>> get_triangulation()
  {
    return _triangulation;
  }

  boost::mpi::communicator get_mpi_communicator() const
  {
    return _communicator;
  }

  std::shared_ptr<
      std::unordered_map<std::string, std::set<dealii::types::material_id>>>
  get_materials() const
  {
    return _materials;
  }

  std::shared_ptr<
      std::unordered_map<std::string, std::set<dealii::types::boundary_id>>>
  get_boundaries() const
  {
    return _boundaries;
  }

  void set_materials(
      std::shared_ptr<std::unordered_map<
          std::string, std::set<dealii::types::material_id>>> materials)
  {
    _materials = materials;
  }

  void set_boundaries(
      std::shared_ptr<std::unordered_map<
          std::string, std::set<dealii::types::boundary_id>>> boundaries)
  {
    _boundaries = boundaries;
  }

private:
  /**
   * Compute the weight used to do load balancing. This is necessary because the
   * physics solved in the collectors, the electrodes, and the separator are
   * different.
   */
  unsigned int compute_cell_weight(
      typename dealii::Triangulation<dim, dim>::cell_iterator const &cell,
      std::array<unsigned int, 4> const &weights) const;

  /**
   * Helper function for the constructor, when the mesh is loaded from a mesh.
   */
  void fill_material_and_boundary_maps(
      std::shared_ptr<boost::property_tree::ptree const> database);

  /**
   * Helper function for the constructor, when the mesh is generated from a
   * database.
   */
  void fill_material_and_boundary_maps();

  /**
   * Convert the geometry database to one that can be used to generate a mesh.
   */
  void convert_geometry_database(
      std::shared_ptr<boost::property_tree::ptree> database);

  /**
   * Create a mesh from a property tree.
   */
  void mesh_generator(boost::property_tree::ptree const &database);

  /**
   * Set the boundary IDs on the cathode and the anode.
   */
  void set_boundary_ids(double const collector_top,
                        double const collector_bottom);

  /**
   * Helper function that serialize and save the geometry in @p filename.
   */
  void output_coarse_mesh(std::string const &filename);

  boost::mpi::communicator _communicator;
  std::shared_ptr<dealii::distributed::Triangulation<dim>> _triangulation;
  std::shared_ptr<std::unordered_map<
      std::string, std::set<dealii::types::material_id>>> _materials;
  std::shared_ptr<std::unordered_map<
      std::string, std::set<dealii::types::boundary_id>>> _boundaries;
  std::unordered_map<std::string, unsigned int> _weights = {};
};
} // end namespace cap

#endif // CAP_GEOMETRY_H
