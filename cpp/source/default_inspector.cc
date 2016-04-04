#include <cap/default_inspector.h>
#include <cap/supercapacitor.h>

namespace cap {

void DefaultInspector::inspect(EnergyStorageDevice *device)
{
    auto super_capacitor = dynamic_cast<SuperCapacitor<2>*>(device);
    if (super_capacitor)
    {
        // get some values from the post processor
        auto post_processor = super_capacitor->get_post_processor();
        double value;
        for (std::string const & key : {
                 "anode_electrode_interfacial_surface_area",
                 "anode_electrode_mass_of_active_material",
                 "cathode_electrode_interfacial_surface_area",
                 "cathode_electrode_mass_of_active_material",
             })
        {
            post_processor->get(key, value);
            _data[key] = value;
        }

        // get other values from the property tree
        auto ptree = super_capacitor->get_property_tree();
        _data["geometric_area"] =
            ptree->get<double>("geometry.geometric_area");
        _data["anode_electrode_thickness"] =
            ptree->get<double>("geometry.anode_electrode_thickness");
        _data["cathode_electrode_thickness"] =
            ptree->get<double>("geometry.cathode_electrode_thickness");
        // TODO
        for (std::string const & electrode : {"anode", "cathode"})
        {
            std::string tmp =
                ptree->get<std::string>("material_properties."+electrode+".matrix_phase");
            _data[electrode+"_electrode_double_layer_capacitance"] =
                ptree->get<double>("material_properties."+tmp+".differential_capacitance");
        }
    }

}

std::map<std::string, double> DefaultInspector::get_data()
{
    return _data;
}

} // end namespace cap

