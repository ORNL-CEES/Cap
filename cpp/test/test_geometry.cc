/* Copyright (c) 2016, the Cap authors.
 *
 * This file is subject to the Modified BSD License and may not be distributed
 * without copyright and license information. Please refer to the file LICENSE
 * for the text and further information on this license.
 */

#define BOOST_TEST_MODULE TestGeometry

#include "main.cc"

#include <cap/energy_storage_device.h>
#include <cap/supercapacitor.h>
#include <cap/utils.h>
#include <cap/geometry.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/info_parser.hpp>
#include <boost/test/unit_test.hpp>
#include <fstream>
#include <unordered_map>

// - Check that a mesh can be loaded, that the areas are computed correctly, and
// check that a mesh can be written.
// - Check that we can build a 3D geometry.

template <int dim>
void write_mesh(
    std::string const &mesh_file,
    std::shared_ptr<dealii::distributed::Triangulation<dim>> triangulation)
{
  dealii::GridOut mesh_writer;
  std::fstream fout;
  fout.open(mesh_file.c_str(), std::fstream::out);
  mesh_writer.write_vtu(*triangulation, fout);
  fout.close();
}

BOOST_AUTO_TEST_CASE(test_reset_geometry)
{
  std::shared_ptr<boost::property_tree::ptree> params(
      new boost::property_tree::ptree);
  params->put("anode_collector_material_id", 4);
  params->put("anode_electrode_material_id", 1);
  params->put("separator_material_id", 2);
  params->put("cathode_electrode_material_id", 3);
  params->put("cathode_collector_material_id", 5);
  params->put("type", "file");
  params->put("mesh_file", "mesh_2d.ucd");
  params->put("materials", 1);
  params->put("material_0.name", "all");
  params->put("material_0.material_id", "1,2,3,4,5");
  params->put("boundaries", 2);
  params->put("boundary_0.name", "anode");
  params->put("boundary_0.boundary_id", "1");
  params->put("boundary_1.name", "cathode");
  params->put("boundary_1.boundary_id", "2");

  cap::Geometry<2> geo(params, boost::mpi::communicator());
  write_mesh("output_test_geometry_0.vtu", geo.get_triangulation());

  dealii::distributed::Triangulation<2> const &tria = *geo.get_triangulation();
  std::cout << "cells=" << tria.n_active_cells() << "  "
            << "faces=" << tria.n_active_faces() << "  "
            << "vertices=" << tria.n_used_vertices() << "\n";

  std::unordered_map<dealii::types::material_id, double> measure;
  auto cell = tria.begin_active();
  auto end_cell = tria.end();
  double m(0.);
  for (; cell != end_cell; ++cell)
    if (cell->is_locally_owned())
    {
      measure[cell->material_id()] += cell->measure();
      m += cell->measure();
    }

  for (auto x : measure)
    std::cout << std::to_string(x.first) << "  " << x.second << "\n";

  double const percent_tolerance = 1.0e-12;
  std::array<double, 5> reference_measure = {6.25e-10, 1.25e-9, 1.25e-9,
                                             1.5e-10, 1.5e-10};
  unsigned int pos(0);
  for (std::string const layer :
       {"separator", "anode_electrode", "cathode_electrode", "anode_collector",
        "cathode_collector"})
  {
    BOOST_CHECK_CLOSE(measure[params->get<dealii::types::material_id>(
                          layer + "_material_id")],
                      reference_measure[pos], percent_tolerance);
    ++pos;
  }
}

BOOST_AUTO_TEST_CASE(test_throw_geometry)
{
  std::shared_ptr<boost::property_tree::ptree> params(
      new boost::property_tree::ptree);
  params->put("type", "supercapacitor");
  params->put("anode_collector_thickness", 5.0e-4);
  params->put("anode_electrode_thickness", 50.0e-4);
  params->put("separator_thickness", 25.0e-4);
  params->put("cathode_electrode_thickness", 50.0e-4);
  params->put("cathode_collector_thickness", 50.0e-4);
  params->put("geometric_area", 25.0e-2);
  params->put("tab_height", 5.0e-4);

  // For now, both collectors must have the same dimensions.
  BOOST_CHECK_THROW(cap::Geometry<2> geo(params, boost::mpi::communicator()),
                    std::runtime_error);
}

BOOST_AUTO_TEST_CASE(test_3d_geometry)
{
  boost::property_tree::ptree device_database;
  boost::property_tree::info_parser::read_info("super_capacitor.info",
                                               device_database);
  boost::property_tree::ptree geometry_database;
  boost::property_tree::info_parser::read_info("generate_mesh.info",
                                               geometry_database);
  std::vector<unsigned int> divisions(3);
  // Collector
  divisions[0] = 3;
  divisions[1] = 3;
  divisions[2] = 3;
  geometry_database.put("collector.divisions", cap::to_string(divisions));
  // Anode
  divisions[0] = 3;
  divisions[1] = 3;
  divisions[2] = 2;
  geometry_database.put("anode.divisions", cap::to_string(divisions));
  // Separator
  divisions[0] = 3;
  divisions[1] = 3;
  divisions[2] = 2;
  geometry_database.put("separator.divisions", cap::to_string(divisions));
  // Cathode
  divisions[0] = 3;
  divisions[1] = 3;
  divisions[2] = 2;
  geometry_database.put("cathode.divisions", cap::to_string(divisions));

  device_database.put_child("geometry", geometry_database);
  device_database.put("dim", 3);

  std::shared_ptr<cap::EnergyStorageDevice> device =
      cap::EnergyStorageDevice::build(device_database,
                                      boost::mpi::communicator());
  std::shared_ptr<cap::SuperCapacitor<3>> supercapacitor =
      std::static_pointer_cast<cap::SuperCapacitor<3>>(device);
  std::shared_ptr<cap::Geometry<3>> geometry = supercapacitor->get_geometry();
  std::shared_ptr<dealii::distributed::Triangulation<3> const> triangulation =
      geometry->get_triangulation();
  const unsigned int n_cells = 6912;
  BOOST_CHECK(n_cells == triangulation->n_active_cells());
}

BOOST_AUTO_TEST_CASE(n_refinements)
{
  boost::property_tree::ptree ptree;
  boost::property_tree::read_info("super_capacitor.info", ptree);
  auto build_geometry = [](boost::property_tree::ptree &ptree)
  {
    return cap::Geometry<2>(std::make_shared<boost::property_tree::ptree>(
                                ptree.get_child("geometry")),
                            boost::mpi::communicator());
  };
  auto geometry_coarse = build_geometry(ptree);
  ptree.put("geometry.n_refinements", 1);
  auto geometry_fine = build_geometry(ptree);
  ptree.put("geometry.n_refinements", 2);
  auto geometry_finer = build_geometry(ptree);

  auto n_cells = [](cap::Geometry<2> &geometry)
  {
    return geometry.get_triangulation()->n_active_cells();
  };
  BOOST_TEST(n_cells(geometry_fine) == 4 * n_cells(geometry_coarse));
  BOOST_TEST(n_cells(geometry_finer) == 16 * n_cells(geometry_coarse));
}

BOOST_AUTO_TEST_CASE(check_no_overlap)
{
  auto mat = std::make_shared<
      std::unordered_map<std::string, std::set<dealii::types::material_id>>>(
      std::initializer_list<
          std::pair<std::string const, std::set<dealii::types::material_id>>>{
          {"foo", std::set<dealii::types::material_id>{1, 2}},
          {"bar", std::set<dealii::types::material_id>{3}}});

  auto bnd = std::make_shared<
      std::unordered_map<std::string, std::set<dealii::types::boundary_id>>>(
      std::initializer_list<
          std::pair<std::string const, std::set<dealii::types::boundary_id>>>{
          {"aaa", std::set<dealii::types::boundary_id>{1, 2}}});

  boost::mpi::communicator world;
  auto tria = std::make_shared<dealii::distributed::Triangulation<2>>(world);

  // No overlap
  BOOST_CHECK_NO_THROW(cap::Geometry<2>(tria, mat, bnd));

  // Material id "3" listed twice
  mat->emplace("already listed", std::set<dealii::types::material_id>{3, 4});
  BOOST_CHECK_THROW(cap::Geometry<2>(tria, mat, bnd), std::runtime_error);
  // Cleanup
  mat->erase("already listed");
  BOOST_CHECK_NO_THROW(cap::Geometry<2>(tria, mat, bnd));

  // Boundary id "1" listed twice
  bnd->emplace("bbb", std::set<dealii::types::boundary_id>{1});
  BOOST_CHECK_THROW(cap::Geometry<2>(tria, mat, bnd), std::runtime_error);
}

// Check that repetitions create a connected triangulation
BOOST_AUTO_TEST_CASE(check_repetitions)
{
  boost::property_tree::ptree geometry_database;
  boost::property_tree::info_parser::read_info("generate_mesh.info",
                                               geometry_database);
  geometry_database.put("n_repetitions", 5);
  cap::Geometry<2> geo(
      std::make_shared<boost::property_tree::ptree>(geometry_database),
      boost::mpi::communicator());

  // Get a sparsity pattern in which nonzero entries indicate that two cells are
  // connected via a common face. Then, check that we can traverse the entire
  // triangulation using cells that share a face.
  std::shared_ptr<dealii::distributed::Triangulation<2>> tria =
      geo.get_triangulation();
  dealii::DynamicSparsityPattern sparsity_pattern;
  dealii::GridTools::get_face_connectivity_of_cells(*tria, sparsity_pattern);

  std::unordered_set<unsigned int> cells_done;
  std::unordered_set<unsigned int> cells_to_do;
  cells_to_do.insert(0);
  unsigned int const n_cells = tria->n_active_cells();
  while (cells_to_do.size() != 0)
  {
    unsigned int const current_cell = *cells_to_do.begin();
    for (unsigned int i = 0; i < n_cells; ++i)
      if (sparsity_pattern.exists(current_cell, i))
        if (cells_done.count(i) == 0)
          cells_to_do.insert(i);

    // The current cell may not be the first element anymore
    std::unordered_set<unsigned int>::iterator cell_it =
        cells_to_do.find(current_cell);
    cells_to_do.erase(cell_it);
    cells_done.insert(current_cell);
  }

  BOOST_CHECK(cells_done.size() == n_cells);
}

BOOST_AUTO_TEST_CASE(test_hyper_trapezoid_2d_geometry)
{
  boost::property_tree::ptree device_database;
  boost::property_tree::info_parser::read_info("super_capacitor.info",
                                               device_database);
  boost::property_tree::ptree geometry_database;
  boost::property_tree::info_parser::read_info("generate_mesh.info",
                                               geometry_database);
  std::vector<unsigned int> divisions(2);
  // Collector
  divisions[0] = 3;
  divisions[1] = 2;
  geometry_database.put("collector.divisions", cap::to_string(divisions));
  // Anode
  divisions[0] = 1;
  divisions[1] = 1;
  geometry_database.put("anode.shape", "hyper_trapezoid_1");
  geometry_database.put("anode.divisions", cap::to_string(divisions));
  // Separator
  geometry_database.put("separator.shape", "hyper_trapezoid_2");
  geometry_database.put("separator.divisions", cap::to_string(divisions));
  // Cathode
  geometry_database.put("cathode.shape", "hyper_trapezoid_3");
  geometry_database.put("cathode.divisions", cap::to_string(divisions));

  device_database.put_child("geometry", geometry_database);

  std::shared_ptr<cap::EnergyStorageDevice> device =
      cap::EnergyStorageDevice::build(device_database,
                                      boost::mpi::communicator());
  std::shared_ptr<cap::SuperCapacitor<2>> supercapacitor =
      std::static_pointer_cast<cap::SuperCapacitor<2>>(device);
  std::shared_ptr<cap::Geometry<2>> geometry = supercapacitor->get_geometry();
  std::shared_ptr<dealii::distributed::Triangulation<2> const> triangulation =
      geometry->get_triangulation();
  const unsigned int n_cells = 240;
  write_mesh("output_test_geometry_1.vtu", geometry->get_triangulation());
  BOOST_CHECK(n_cells == triangulation->n_active_cells());
}
