/* Copyright (c) 2016, the Cap authors.
 *
 * This file is subject to the Modified BSD License and may not be distributed
 * without copyright and license information. Please refer to the file LICENSE
 * for the text and further information on this license.
 */

#define BOOST_TEST_MODULE MaterialPropertyValues

#include "main.cc"

#include <cap/geometry.h>
#include <cap/mp_values.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/info_parser.hpp>
#include <boost/test/unit_test.hpp>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/dofs/dof_handler.h>

// demonstrate how unifom constant and composite properties are used
BOOST_AUTO_TEST_CASE(mp_values)
{
  int constexpr dim = 2;
  std::shared_ptr<cap::MPValues<dim>> mp_values;

  // make a dummy triangulation to get an active cell iterator
  dealii::Triangulation<dim> triangulation;
  dealii::GridGenerator::hyper_cube(triangulation);
  dealii::DoFHandler<dim> dof_handler(triangulation);
  dealii::DoFHandler<dim>::active_cell_iterator cell =
      dof_handler.begin_active();

  // create a vector
  std::vector<double> values;
  // note that `MPValues::get_values(...)` does not allocate memory

  // build a constant uniform property
  double const value = 3.14;
  mp_values = std::make_shared<cap::UniformConstantMPValues<dim>>(value);

  // key property name is ignored for constant uniform property
  mp_values->get_values("doesnotmatter", cell, values);
  // vector is still empty
  BOOST_TEST(values.empty());
  // now resize it and try again
  values.resize(3, 0.0);
  mp_values->get_values("", cell, values);
  // all constant uniform property do is a `std::fill(...)`
  for (auto const &v : values)
    BOOST_TEST(v == value);

  // define a custom composite property
  class MyCompositeProMPValues : public cap::CompositePro<2>
  {
  public:
    MyCompositeProMPValues()
    {
      for (auto const &x : _map)
        (this->_properties)
            .emplace(x.first, std::make_shared<cap::UniformConstantMPValues<2>>(
                                  x.second));
    }
    std::map<std::string, double> _map = {{"foo", 1.0}, {"bar", 2.0}};
  };

  // build it
  mp_values = std::make_shared<MyCompositeProMPValues>();

  // if the property is not registered in the `this->_properties` map an
  // exception will be thrown
  BOOST_CHECK_THROW(mp_values->get_values("notregistered", cell, values),
                    std::runtime_error);

  // the composite property just redirects the call
  for (auto const &x :
       std::dynamic_pointer_cast<MyCompositeProMPValues>(mp_values)->_map)
  {
    mp_values->get_values(x.first, cell, values);
    for (auto const &v : values)
      BOOST_TEST(v == x.second);
  }

  // define another custom composite property
  // this time the map is based on `material_id`s
  class MyCompositeMatMPValues : public cap::CompositeMat<2>
  {
  public:
    MyCompositeMatMPValues()
    {
      for (auto const &x : _map)
        (this->_materials)
            .emplace(x.first, std::make_shared<cap::UniformConstantMPValues<2>>(
                                  x.second));
    }
    std::map<dealii::types::material_id, double> _map = {
        {3, 3.0}, {4, 4.0}, {5, 5.0}};
  };

  // build it
  mp_values = std::make_shared<MyCompositeMatMPValues>();

  // if the material is not registered in the `this->_materials` map an
  // exception will be thrown
  cell->set_material_id(1);
  BOOST_CHECK_THROW(mp_values->get_values("willbeignored", cell, values),
                    std::runtime_error);

  // the composite property just redirects the call
  for (auto const &x :
       std::dynamic_pointer_cast<MyCompositeMatMPValues>(mp_values)->_map)
  {
    cell->set_material_id(x.first);
    mp_values->get_values("", cell, values);
    for (auto const &v : values)
      BOOST_TEST(v == x.second);
  }
}

BOOST_AUTO_TEST_CASE(test_mp_values_throw)
{
  // ensure missing database will not segfault
  BOOST_CHECK_THROW(cap::SuperCapacitorMPValuesFactory<2>::build(
                        cap::MPValuesParameters<2>(nullptr)),
                    std::runtime_error);

  // same if geometry is missing
  auto empty_database = std::make_shared<boost::property_tree::ptree>();
  BOOST_CHECK_THROW(cap::SuperCapacitorMPValuesFactory<2>::build(
                        cap::MPValuesParameters<2>(empty_database)),
                    std::runtime_error);
}

BOOST_AUTO_TEST_CASE(test_inhomogenous_mp_values)
{
  boost::mpi::communicator world;

  // Use the standard input for build a SuperCapacitor
  boost::property_tree::ptree ptree;
  boost::property_tree::read_info("super_capacitor.info", ptree);

  // Extract the section for the material properties
  auto database = std::make_shared<boost::property_tree::ptree>(
      ptree.get_child("material_properties"));

  // Create an standard homogeneous SuperCapacitorMPValues using the build(...)
  // method which is static.
  cap::MPValuesParameters<2> params(database);
  params.geometry = std::make_shared<cap::Geometry<2>>(
      std::make_shared<boost::property_tree::ptree>(
          ptree.get_child("geometry")),
      world);
  std::shared_ptr<cap::MPValues<2>> mp_values =
      cap::SuperCapacitorMPValuesFactory<2>::build(params);
  BOOST_TEST(
      std::dynamic_pointer_cast<cap::SuperCapacitorMPValues<2>>(mp_values));

  // Now modify the material properties database to create the inhomogeneous
  // version of the MPValues.
  database->put("inhomogeneous", true);
  // NOTE: if parameters == 0, SuperCapacitorMPValues::build(...) still
  // instancies the map cell_id -> MPValues but the computation should be
  // identical to the homogeneous case that has a map material_id -> MPValues
  database->put("parameters", 2);
  database->put("parameter_0.path", "separator_material.void_volume_fraction");
  database->put("parameter_0.distribution_type", "uniform");
  database->put("parameter_0.range", "0.55, 0.65");
  database->put("parameter_1.path", "electrode_material.void_volume_fraction");
  database->put("parameter_1.distribution_type", "normal");
  database->put("parameter_1.mean", "0.67");
  database->put("parameter_1.standard_deviation", "0.04");
  mp_values = cap::SuperCapacitorMPValuesFactory<2>::build(params);
  BOOST_TEST(
      std::dynamic_pointer_cast<cap::InhomogeneousSuperCapacitorMPValues<2>>(
          mp_values));

  // Check that an exception is thrown if the same path is registered for
  // multiple parameters.
  database->put("parameter_1.path", "separator_material.void_volume_fraction");
  BOOST_CHECK_THROW(cap::SuperCapacitorMPValuesFactory<2>::build(params),
                    std::runtime_error);
  // Undo it.
  database->put("parameter_1.path", "electrode_material.void_volume_fraction");
  BOOST_CHECK_NO_THROW(cap::SuperCapacitorMPValuesFactory<2>::build(params));

  // Check that an exception is thrown the parameter path does not already exist
  // in the database. It most likely means that it is typo and we want to catch
  // that.
  database->put("parameter_0.path", "electrode_material.does_not_exist");
  BOOST_CHECK_THROW(cap::SuperCapacitorMPValuesFactory<2>::build(params),
                    std::runtime_error);
}

BOOST_AUTO_TEST_CASE(test_mp_values)
{
  // Fill in the geometry database
  std::shared_ptr<boost::property_tree::ptree> geometry_database(
      new boost::property_tree::ptree());

  boost::property_tree::ptree material_0_database;
  material_0_database.put("name", "anode");
  material_0_database.put("material_id", 1);

  boost::property_tree::ptree material_1_database;
  material_1_database.put("name", "separator");
  material_1_database.put("material_id", 2);

  boost::property_tree::ptree material_2_database;
  material_2_database.put("name", "cathode");
  material_2_database.put("material_id", 3);

  boost::property_tree::ptree material_3_database;
  material_3_database.put("name", "collector");
  material_3_database.put("material_id", 4);

  geometry_database->put("type", "file");
  geometry_database->put("mesh_file", "mesh_2d.ucd");
  geometry_database->put("anode_collector_thickness", 5.0e-4);
  geometry_database->put("anode_electrode_thickness", 50.0e-4);
  geometry_database->put("separator_thickness", 25.0e-4);
  geometry_database->put("cathode_electrode_thickness", 50.0e-4);
  geometry_database->put("cathode_collector_thickness", 5.0e-4);
  geometry_database->put("geometric_area", 25.0e-2);
  geometry_database->put("tab_height", 5.0e-4);
  geometry_database->put("materials", 4);
  geometry_database->put("anode_collector_material_id", 4);
  geometry_database->put("anode_electrode_material_id", 1);
  geometry_database->put("separator_material_id", 2);
  geometry_database->put("cathode_electrode_material_id", 3);
  geometry_database->put("cathode_collector_material_id", 5);
  geometry_database->put_child("material_0", material_0_database);
  geometry_database->put_child("material_1", material_1_database);
  geometry_database->put_child("material_2", material_2_database);
  geometry_database->put_child("material_3", material_3_database);
  geometry_database->put("boundary_values.anode_boundary_id", "1");
  geometry_database->put("boundary_values.cathode_boundary_id", "2");

  std::shared_ptr<cap::Geometry<2>> geometry =
      std::make_shared<cap::Geometry<2>>(geometry_database,
                                         boost::mpi::communicator());

  // Fill in the material properties database
  std::shared_ptr<boost::property_tree::ptree> material_properties_database(
      new boost::property_tree::ptree());

  boost::property_tree::ptree anode_database;
  anode_database.put("type", "porous_electrode");
  anode_database.put("matrix_phase", "electrode_material");
  anode_database.put("solution_phase", "electrolyte");

  boost::property_tree::ptree cathode_database;
  cathode_database.put("type", "porous_electrode");
  cathode_database.put("matrix_phase", "electrode_material");
  cathode_database.put("solution_phase", "electrolyte");

  boost::property_tree::ptree separator_database;
  separator_database.put("type", "permeable_membrane");
  separator_database.put("matrix_phase", "separator_material");
  separator_database.put("solution_phase", "electrolyte");

  boost::property_tree::ptree collector_database;
  collector_database.put("type", "current_collector");
  collector_database.put("metal_foil", "collector_material");

  boost::property_tree::ptree separator_material_database;
  separator_material_database.put("void_volume_fraction", 0.6);
  separator_material_database.put("tortuosity_factor", 1.29);
  separator_material_database.put("pores_characteristic_dimension", 1.5e-7);
  separator_material_database.put("pores_geometry_factor", 2.0);
  separator_material_database.put("mass_density", 3.2);
  separator_material_database.put("heat_capacity", 1.2528e3);
  separator_material_database.put("thermal_conductivity", 0.0019e2);

  boost::property_tree::ptree electrode_material_database;
  electrode_material_database.put("differential_capacitance", 3.134);
  electrode_material_database.put("exchange_current_density", 7.463e-10);
  electrode_material_database.put("void_volume_fraction", 0.67);
  electrode_material_database.put("tortuosity_factor", 2.3);
  electrode_material_database.put("pores_characteristic_dimension", 1.5e-7);
  electrode_material_database.put("pores_geometry_factor", 2.0);
  electrode_material_database.put("mass_density", 2.3);
  electrode_material_database.put("electrical_resistivity", 1.92);
  electrode_material_database.put("heat_capacity", 0.93e3);
  electrode_material_database.put("thermal_conductivity", 0.0011e2);

  boost::property_tree::ptree collector_material_database;
  collector_material_database.put("mass_density", 2.7);
  collector_material_database.put("electrical_resistivity", 28.2e-7);
  collector_material_database.put("heat_capacity", 2.7e3);
  collector_material_database.put("thermal_conductivity", 237.0);

  boost::property_tree::ptree electrolyte_database;
  electrolyte_database.put("mass_density", 1.2);
  electrolyte_database.put("electrical_resistivity", 1.49e3);
  electrolyte_database.put("heat_capacity", 0.0);
  electrolyte_database.put("thermal_conductivity", 0.0);

  material_properties_database->put_child("anode", anode_database);
  material_properties_database->put_child("cathode", cathode_database);
  material_properties_database->put_child("separator", separator_database);
  material_properties_database->put_child("collector", collector_database);
  material_properties_database->put_child("separator_material",
                                          separator_material_database);
  material_properties_database->put_child("electrode_material",
                                          electrode_material_database);
  material_properties_database->put_child("collector_material",
                                          collector_material_database);
  material_properties_database->put_child("electrolyte", electrolyte_database);

  cap::MPValuesParameters<2> params(material_properties_database);
  params.geometry = geometry;
  std::shared_ptr<cap::MPValues<2>> mp_values =
      std::make_shared<cap::SuperCapacitorMPValues<2>>(params);

  dealii::DoFHandler<2> dof_handler(*(geometry->get_triangulation()));
  dealii::DoFHandler<2>::active_cell_iterator cell = dof_handler.begin_active();

  std::vector<double> values(1);
  double const tolerance = 1e-2;
  mp_values->get_values("density", cell, values);
  BOOST_TEST(values[0] == 1563.);
  mp_values->get_values("solid_electrical_conductivity", cell, values);
  BOOST_TEST(std::abs(values[0] - 17.1785) < tolerance);

  // Move to another material
  for (unsigned int i = 0; i < 300; ++i)
    ++cell;

  mp_values->get_values("density", cell, values);
  BOOST_TEST(values[0] == 2000.);
  mp_values->get_values("solid_electrical_conductivity", cell, values);
  BOOST_TEST(std::abs(values[0]) == 0.);
}
