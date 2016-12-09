/* Copyright (c) 2016, the Cap authors.
 *
 * This file is subject to the Modified BSD License and may not be distributed
 * without copyright and license information. Please refer to the file LICENSE
 * for the text and further information on this license.
 */

#define BOOST_TEST_MODULE TestPostProcessor

#include "main.cc"

#include <cap/post_processor.h>
#include <cap/supercapacitor.h>
#include <cap/utils.h>
#include <deal.II/base/exceptions.h>
#include <boost/format.hpp>
#include <boost/foreach.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/info_parser.hpp>
#include <boost/test/unit_test.hpp>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdio>

BOOST_AUTO_TEST_CASE(test_compute_energy)
{
  // We need to use a large tolerance because the error increases at each time
  // step.
  double const tolerance = 1e-6;

  double const pi = M_PI;
  double const initial_time = 0.0;
  double const final_time = 2.0 * pi;
  std::size_t const n = 10001;
  std::vector<double> time(n);
  std::vector<double> power(n);

  // INITIALIZE SOLUTION
  for (std::size_t i = 0; i < n; ++i)
  {
    time[i] = initial_time +
              static_cast<double>(i) / (n - 1) * (final_time - initial_time);
    power[i] = std::cos(time[i]);
  }

  // COMPUTE ENERGY
  std::vector<double> energy(n);
  cap::compute_energy(time, power, energy);

  // CHECK THE ANSWER
  std::vector<double> exact(n);
  std::vector<double> error(n);
  exact[0] = 0.0;
  for (std::size_t i = 0; i < n; ++i)
  {
    exact[i] = std::sin(time[i]);
    error[i] = energy[i] - exact[i];
    BOOST_REQUIRE(std::abs(error[i]) < tolerance);
  }
}

// Check that the postprocessor output the quantities that were asked
BOOST_AUTO_TEST_CASE(test_supercapacitor_inspector)
{
  // Build an energy storage device
  boost::property_tree::ptree ptree;
  boost::property_tree::info_parser::read_info("super_capacitor.info", ptree);
  ptree.put("debug.material_ids", true);
  ptree.put("debug.material_properties",
            "solid_electrical_conductivity,liquid_electrical_conductivity");
  ptree.put("debug.solution_fields",
            "solid_potential,liquid_potential,overpotential,joule_heating");
  ptree.put("debug.solution_fluxes",
            "solid_current_density,liquid_current_density");
  boost::mpi::communicator world;
  std::shared_ptr<cap::EnergyStorageDevice> device =
      cap::EnergyStorageDevice::build(ptree, world);

  // Output the vtu file
  cap::SuperCapacitorInspector<2> supercap_inspector;
  supercap_inspector.inspect(device.get());

  if (world.rank() == 0)
  {
    // Check for that the fields are present in the output file
    unsigned int fields_found = 0;
    std::vector<std::string> fields;
    fields.push_back("material_id");
    fields.push_back("liquid_electrical_conductivity");
    fields.push_back("solid_electrical_conductivity");
    fields.push_back("solid_potential");
    fields.push_back("liquid_potential");
    fields.push_back("overpotential");
    fields.push_back("solid_current_density_0");
    fields.push_back("solid_current_density_1");
    fields.push_back("liquid_current_density_0");
    fields.push_back("liquid_current_density_1");
    fields.push_back("joule_heating");
    std::ifstream file;
    std::string line;
    std::string filename("solution-0000.0000.vtu");
    file.open(filename);
    if (file.is_open())
    {
      while (!file.eof())
      {
        std::getline(file, line);
        for (auto const &field : fields)
          if (line.find(field) != std::string::npos)
            ++fields_found;
      }
      file.close();
    }
    BOOST_CHECK_EQUAL(fields_found, fields.size());

    // Delete file
    std::remove(filename.c_str());
  }
}

// Check that the postprocessor throws when the input is invalid
BOOST_AUTO_TEST_CASE(test_throw)
{
  // Build an energy storage device
  boost::property_tree::ptree ptree;
  boost::property_tree::info_parser::read_info("super_capacitor.info", ptree);
  ptree.put("debug.material_properties", "invalid");
  boost::mpi::communicator world;
  std::string const prefix = "debug.";

  {
    // Check that the preprocessor throws with an invalid key
    ptree.put(prefix + "material_properties", "invalid");
    BOOST_CHECK_THROW(cap::EnergyStorageDevice::build(ptree, world),
                      std::runtime_error);
    // Clear and make sure it does not throw any more
    ptree.put(prefix + "material_properties", "");
    BOOST_CHECK_NO_THROW(cap::EnergyStorageDevice::build(ptree, world));
  }

  for (std::string const &s : {"solution_fields", "solution_fluxes"})
  {
    // Check that the preprocessor throws with an invalid key
    ptree.put(prefix + s, "invalid");
    BOOST_CHECK_THROW(cap::EnergyStorageDevice::build(ptree, world),
                      dealii::ExceptionBase);
    // Clear and make sure it does not throw any more
    ptree.put(prefix + s, "");
    BOOST_CHECK_NO_THROW(cap::EnergyStorageDevice::build(ptree, world));
  }
}
