#include <cap/default_inspector.h>
#include <cap/supercapacitor.h>

namespace cap {

void DefaultInspector::inspect(EnergyStorageDevice *device)
{
    auto super_capacitor = dynamic_cast<SuperCapacitor<2>*>(device);
    if (super_capacitor)
    {
        auto post_processor = super_capacitor->get_post_processor();
        double value;
        post_processor->get("interfacial_surface_area", value);
        _data["surface_area"] = value;
        post_processor->get("mass", value);
        _data["mass"] = value;

        auto ptree = super_capacitor->get_property_tree();
        _data["geometric_area"] =
            ptree->get<double>("geometry.geometric_area");
        _data["electrode_thickness"] =
            ptree->get<double>("geometry.anode_electrode_thickness");
        // TODO
        std::string tmp =
            ptree->get<std::string>("material_properties.anode.matrix_phase");
        _data["capacitance"] =
            ptree->get<double>("material_properties."+tmp+".differential_capacitance");
    }

}

std::map<std::string, double> DefaultInspector::get_data()
{
    return _data;
}

} // end namespace cap

