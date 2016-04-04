/* Copyright (c) 2016, the Cap authors.
 *
 * This file is subject to the Modified BSD License and may not be distributed
 * without copyright and license information. Please refer to the file LICENSE
 * for the text and further information on this license.
 */

#include <cap/resistor_capacitor.h>
#ifdef WITH_DEAL_II
#include <cap/new_supercapacitor.h>
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

EnergyStorageDevice::~EnergyStorageDevice() = default;

} // end namespace cap
