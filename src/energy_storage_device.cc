#include <cap/resistor_capacitor.h>

namespace cap {

std::shared_ptr<cap::EnergyStorageDevice>
buildEnergyStorageDevice(std::shared_ptr<cap::Parameters const> params)
{
    std::shared_ptr<boost::property_tree::ptree const> database = params->database;
    std::string const type = database->get<std::string>("type", "unknown_type");
    if (type.compare("SeriesRC") == 0)
        return std::make_shared<cap::SeriesRC>(params);
    else if (type.compare("ParallelRC") == 0)
        return std::make_shared<cap::ParallelRC>(params);
    else
        throw std::runtime_error("invalid energy storage type ``"+type+"''\n");
}

} // end namespace cap

