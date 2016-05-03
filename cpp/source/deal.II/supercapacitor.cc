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
  build(boost::property_tree::ptree const &ptree,
        boost::mpi::communicator const &comm) override
  {
    int const dim = ptree.get<int>("dim");
    if (dim == 2)
      return std::make_unique<SuperCapacitor<2>>(
          SuperCapacitor<2>(ptree, comm));
    else if (dim == 3)
      return std::make_unique<SuperCapacitor<3>>(
          SuperCapacitor<3>(ptree, comm));
    else
      throw std::runtime_error("dim=" + std::to_string(dim) +
                               " must be 2 or 3");
  }
} global_SuperCapacitorBuilder;

template class SuperCapacitorInspector<2>;
template class SuperCapacitorInspector<3>;
template class SuperCapacitor<2>;
template class SuperCapacitor<3>;

} // end namespace cap
