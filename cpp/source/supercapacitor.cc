/* Copyright (c) 2016, the Cap authors.
 *
 * This file is subject to the Modified BSD License and may not be distributed
 * without copyright and license information. Please refer to the file LICENSE
 * for the text and further information on this license. 
 */

#include <cap/supercapacitor.templates.h>

namespace cap
{

class SuperCapacitorBuilder : public EnergyStorageDeviceBuilder
{
public:
  SuperCapacitorBuilder()
  {
    register_energy_storage_device("SuperCapacitor", this);
  }
  std::unique_ptr<EnergyStorageDevice>
  build(boost::mpi::communicator const &comm,
        boost::property_tree::ptree const &ptree)
  {
    int const dim = ptree.get<int>("dim");
    if (dim == 2)
      return std::unique_ptr<SuperCapacitor<2>>(
          new SuperCapacitor<2>(comm, ptree));
    else if (dim == 3)
      return std::unique_ptr<SuperCapacitor<3>>(
          new SuperCapacitor<3>(comm, ptree));
    else
      throw std::runtime_error("dim=" + std::to_string(dim) +
                               " must be 2 or 3");
  }
} gloabal_SuperCapacitorBuilder;

template class SuperCapacitor<2>;
template class SuperCapacitor<3>;

} // end namespace cap
