/* Copyright (c) 2016, the Cap authors.
 *
 * This file is subject to the Modified BSD License and may not be distributed
 * without copyright and license information. Please refer to the file LICENSE
 * for the text and further information on this license. 
 */

#include <cap/deal.II/new_supercapacitor.templates.h>

namespace cap
{

class New_SuperCapacitorBuilder : public EnergyStorageDeviceBuilder
{
public:
  New_SuperCapacitorBuilder()
  {
    register_energy_storage_device("New_SuperCapacitor", this);
  }
  std::unique_ptr<EnergyStorageDevice>
  build(boost::mpi::communicator const &comm,
        boost::property_tree::ptree const &ptree)
  {
    int const dim = ptree.get<int>("dim");
    if (dim == 2)
      return std::unique_ptr<New_SuperCapacitor<2>>(
          new New_SuperCapacitor<2>(comm, ptree));
    else if (dim == 3)
      return std::unique_ptr<New_SuperCapacitor<3>>(
          new New_SuperCapacitor<3>(comm, ptree));
    else
      throw std::runtime_error("dim=" + std::to_string(dim) +
                               " must be 2 or 3");
  }
} global_New_SuperCapacitorBuilder;

template class New_SuperCapacitor<2>;
template class New_SuperCapacitor<3>;

} // end namespace cap
