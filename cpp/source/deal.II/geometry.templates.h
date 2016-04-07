/* Copyright (c) 2016, the Cap authors.
 *
 * This file is subject to the Modified BSD License and may not be distributed
 * without copyright and license information. Please refer to the file LICENSE
 * for the text and further information on this license.
 */

#include <cap/geometry.h>
#include <cap/utils.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/base/geometry_info.h>
#include <fstream>
#include <tuple>

namespace cap
{
template <int dim>
void Geometry<dim>::fill_materials_map(
    std::shared_ptr<boost::property_tree::ptree const> database)
{
  this->materials = std::make_shared<std::unordered_map<
      std::string, std::vector<dealii::types::material_id>>>();
  int const n_materials = database->get<int>("materials");
  for (int m = 0; m < n_materials; ++m)
  {
    boost::property_tree::ptree material_database =
        database->get_child("material_" + std::to_string(m));
    std::vector<dealii::types::material_id> const material_ids =
        to_vector<dealii::types::material_id>(
            material_database.get<std::string>("material_id"));
    std::string const material_name =
        material_database.get<std::string>("name");
    (*this->materials).emplace(material_name, material_ids);
  }
}

template <int dim>
Geometry<dim>::Geometry(
    std::shared_ptr<boost::property_tree::ptree const> database,
    boost::mpi::communicator mpi_communicator)
    : mpi_communicator(mpi_communicator)
{
  this->triangulation =
      std::make_shared<dealii::distributed::Triangulation<dim>>(
          mpi_communicator);
  std::string const mesh_file = database->get<std::string>("mesh_file");
  dealii::GridIn<dim> mesh_reader;
  mesh_reader.attach_triangulation(*(this->triangulation));
  std::fstream fin;
  fin.open(mesh_file.c_str(), std::fstream::in);
  std::string const file_extension =
      mesh_file.substr(mesh_file.find_last_of(".") + 1);
  if (file_extension.compare("ucd") == 0)
  {
    mesh_reader.read_ucd(fin);
  }
  else
  {
    throw std::runtime_error("Bad mesh file extension ." + file_extension +
                             " in mesh file " + mesh_file);
  }
  fin.close();

  fill_materials_map(database);
}

template <int dim>
Geometry<dim>::Geometry(
    std::shared_ptr<boost::property_tree::ptree const> database,
    std::shared_ptr<dealii::distributed::Triangulation<dim>> triangulation)
    : mpi_communicator(boost::mpi::communicator(
          triangulation->get_communicator(), boost::mpi::comm_duplicate)),
      triangulation(triangulation)
{
  fill_materials_map(database);
}
} // end namespace cap
