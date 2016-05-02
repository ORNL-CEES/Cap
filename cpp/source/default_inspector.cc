/* Copyright (c) 2016, the Cap authors.
 *
 * This file is subject to the Modified BSD License and may not be distributed
 * without copyright and license information. Please refer to the file LICENSE
 * for the text and further information on this license.
 */

#include <cap/default_inspector.h>
#include <cap/supercapacitor.h>

namespace cap
{

DefaultInspector::DefaultInspector() : _data(std::map<std::string, double>()) {}

template <int dim>
std::map<std::string, double>
extract_data_from_super_capacitor(EnergyStorageDevice *device)
{
  std::map<std::string, double> data;
  auto super_capacitor = dynamic_cast<SuperCapacitor<dim> *>(device);
  if (super_capacitor)
  {
    // get some values from the post processor
    auto post_processor = super_capacitor->get_post_processor();
    BOOST_ASSERT_MSG(post_processor != nullptr,
                     "The Postprocessor does not exist.");
    double value;
    for (std::string const &key : {
             "anode_electrode_interfacial_surface_area",
             "anode_electrode_mass_of_active_material",
             "cathode_electrode_interfacial_surface_area",
             "cathode_electrode_mass_of_active_material",
         })
    {
      post_processor->get(key, value);
      data[key] = value;
    }

    // get other values from the property tree
    boost::property_tree::ptree const *ptree =
        super_capacitor->get_property_tree();
    data["geometric_area"] = ptree->get<double>("geometry.geometric_area");
    data["anode_electrode_thickness"] =
        ptree->get<double>("geometry.anode_electrode_thickness");
    data["cathode_electrode_thickness"] =
        ptree->get<double>("geometry.cathode_electrode_thickness");
    for (std::string const &electrode : {"anode", "cathode"})
    {
      std::string tmp = ptree->get<std::string>("material_properties." +
                                                electrode + ".matrix_phase");
      data[electrode + "_electrode_double_layer_capacitance"] =
          ptree->get<double>("material_properties." + tmp +
                             ".differential_capacitance");
    }
  }
  else
  {
    throw std::runtime_error("Downcasting failed");
  }
  return data;
}

void DefaultInspector::inspect(EnergyStorageDevice *device)
{
  if (dynamic_cast<SuperCapacitor<2> *>(device))
  {
    _data = extract_data_from_super_capacitor<2>(device);
  }
  else if (dynamic_cast<SuperCapacitor<3> *>(device))
  {
    _data = extract_data_from_super_capacitor<3>(device);
  }
  else
  {
    // do nothing
  }
}

std::map<std::string, double> DefaultInspector::get_data() { return _data; }

} // end namespace cap
