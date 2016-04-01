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
        _data["interfacial_surface_area"] = value;
        post_processor->get("mass", value);
        _data["mass"] = value;
    }

}

std::map<std::string, double> DefaultInspector::get_data()
{
    return _data;
}

} // end namespace cap

