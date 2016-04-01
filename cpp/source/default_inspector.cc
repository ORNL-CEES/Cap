#include <cap/default_inspector.h>

namespace cap {

void DefaultInspector::inspect(EnergyStorageDevice *device)
{ }

std::map<std::string, double> DefaultInspector::get_data()
{
    return _data;
}

} // end namespace cap

