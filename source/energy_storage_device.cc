#include <cap/resistor_capacitor.h>
#ifdef WITH_DEAL_II
#include <cap/supercapacitor.h>
#endif

namespace cap
{

EnergyStorageDeviceInspector::~EnergyStorageDeviceInspector() = default;

EnergyStorageDeviceBuilder::~EnergyStorageDeviceBuilder() = default;

void EnergyStorageDeviceBuilder::register_energy_storage_device(
    std::string const &type, EnergyStorageDeviceBuilder *builder)
{
  EnergyStorageDevice::_builders[type] = builder;
}

std::map<std::string, EnergyStorageDeviceBuilder *>
    EnergyStorageDevice::_builders =
        std::map<std::string, EnergyStorageDeviceBuilder *>();

std::unique_ptr<EnergyStorageDevice>
EnergyStorageDevice::build(boost::mpi::communicator const &comm,
                           boost::property_tree::ptree const &ptree)
{
  auto const type = ptree.get<std::string>("type");
  auto const it = _builders.find(type);
  if (it != _builders.end())
    return (it->second)->build(comm, ptree);
  else
    throw std::runtime_error("invalid EnergyStorageDevice type `" + type +
                             "`\n");
}

EnergyStorageDevice::EnergyStorageDevice(
    boost::mpi::communicator const &communicator)
    : _communicator(communicator)
{
}

std::shared_ptr<EnergyStorageDevice>
buildEnergyStorageDevice(boost::mpi::communicator const &comm,
                         boost::property_tree::ptree const &ptree)
{
  std::string const type = ptree.get<std::string>("type", "unknown_type");
  if (type.compare("SeriesRC") == 0)
    return std::make_shared<cap::SeriesRC>(comm, ptree);
  else if (type.compare("ParallelRC") == 0)
    return std::make_shared<cap::ParallelRC>(comm, ptree);
  else if (type.compare("SuperCapacitor") == 0)
  {
#ifdef WITH_DEAL_II
    int const dim = ptree.get<int>("dim");
    if (dim == 2)
      return std::make_shared<cap::SuperCapacitor<2>>(comm, ptree);
    else if (dim == 3)
      return std::make_shared<cap::SuperCapacitor<3>>(comm, ptree);
    else
      throw std::runtime_error("dim=" + std::to_string(dim) +
                               " must be 2 or 3");
#else
    throw std::runtime_error("reconfigure with deal.II");
#endif
  }
  else
    throw std::runtime_error("invalid energy storage type ``" + type + "''\n");
}

EnergyStorageDevice::~EnergyStorageDevice() = default;

} // end namespace cap
