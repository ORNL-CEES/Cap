/* Copyright (c) 2016, the Cap authors.
 *
 * This file is subject to the Modified BSD License and may not be distributed
 * without copyright and license information. Please refer to the file LICENSE
 * for the text and further information on this license.
 */

#ifndef CAP_DEFAULT_INSPECTOR_H
#define CAP_DEFAULT_INSPECTOR_H

#include <cap/energy_storage_device.h>

namespace cap
{

class DefaultInspector : public EnergyStorageDeviceInspector
{
public:
  void inspect(EnergyStorageDevice *device) override;
  std::map<std::string, double> get_data();

protected:
  std::map<std::string, double> _data;
};

} // end namespace cap

#endif // CAP_DEFAULT_INSPECTOR_H
