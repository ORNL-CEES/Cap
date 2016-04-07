/* Copyright (c) 2016, the Cap authors.
 *
 * This file is subject to the Modified BSD License and may not be distributed
 * without copyright and license information. Please refer to the file LICENSE
 * for the text and further information on this license.
 */

#define BOOST_TEST_MODULE TestGeometry

#include "main.cc"

#include <cap/utils.h>
#include <cap/geometry.h>
#include <deal.II/grid/grid_out.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/test/unit_test.hpp>
#include <fstream>
#include <unordered_map>

// Check that a mesh can be loaded, that the areas are computed correctly, and
// check that a mesh can be written.

template <int dim>
void write_mesh(std::string const &mesh_file,
                std::shared_ptr<dealii::distributed::Triangulation<dim> const>
                    triangulation)
{
  dealii::GridOut mesh_writer;
  std::fstream fout;
  fout.open(mesh_file.c_str(), std::fstream::out);
  std::string const file_extension =
      mesh_file.substr(mesh_file.find_last_of(".") + 1);
  if (file_extension.compare("vtu") == 0)
  {
    mesh_writer.write_vtu(*triangulation, fout);
  }
  else
  {
    throw std::runtime_error("Bad output format ." + file_extension +
                             " in mesh file " + mesh_file);
  }
  fout.close();
}

BOOST_AUTO_TEST_CASE(test_reset_geometry)
{
  std::shared_ptr<boost::property_tree::ptree> params(
      new boost::property_tree::ptree);
  params->put("anode_collector_thickness", 45.0e-4);
  params->put("anode_electrode_thickness", 50.0e-4);
  params->put("separator_thickness", 5.0e-4);
  params->put("cathode_electrode_thickness", 250.0e-4);
  params->put("cathode_collector_thickness", 15.0e-4);
  params->put("geometric_area", 25.0e-2);
  params->put("tab_height", 20.0e-4);
  params->put("anode_collector_material_id", 4);
  params->put("anode_electrode_material_id", 1);
  params->put("separator_material_id", 2);
  params->put("cathode_electrode_material_id", 3);
  params->put("cathode_collector_material_id", 5);
  params->put("mesh_file", "mesh_2d.ucd");
  params->put("materials", 1);
  params->put("material_0.name", "all");
  params->put("material_0.material_id", "1,2,3,4,5");

  cap::Geometry<2> geo(params, boost::mpi::communicator());
  write_mesh("output_test_geometry_0.vtu", geo.get_triangulation());

  dealii::Triangulation<2> const &tria = *geo.get_triangulation();
  std::cout << "cells=" << tria.n_active_cells() << "  "
            << "faces=" << tria.n_active_faces() << "  "
            << "vertices=" << tria.n_used_vertices() << "\n";

  std::unordered_map<dealii::types::material_id, double> measure;
  auto cell = tria.begin_active();
  auto end_cell = tria.end();
  for (; cell != end_cell; ++cell)
  {
    measure[cell->material_id()] += cell->measure();
  }

  for (auto x : measure)
    std::cout << std::to_string(x.first) << "  " << x.second << "\n";

  double const percent_tolerance = 1.0e-12;
  for (std::string const layer :
       {"separator", "anode_electrode", "cathode_electrode"})
    BOOST_CHECK_CLOSE(measure[params->get<dealii::types::material_id>(
                          layer + "_material_id")],
                      0.01 * params->get<double>(layer + "_thickness") *
                          0.0001 * params->get<double>("geometric_area"),
                      percent_tolerance);
  for (std::string const layer : {"anode_collector", "cathode_collector"})
    BOOST_CHECK_CLOSE(measure[params->get<dealii::types::material_id>(
                          layer + "_material_id")],
                      0.01 * params->get<double>(layer + "_thickness") *
                          (0.0001 * params->get<double>("geometric_area") +
                           0.01 * params->get<double>("tab_height")),
                      percent_tolerance);
}
