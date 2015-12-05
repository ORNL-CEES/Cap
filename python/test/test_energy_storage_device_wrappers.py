from pycap import PropertyTree,EnergyStorageDevice
from numpy import array,pi,linalg,inf
import unittest

valid_device_input=[
    "series_rc.info",
    "parallel_rc.info",
#    "super_capacitor.info",
    ]

class capEnergyStorageDeviceWrappersTestCase(unittest.TestCase):
    def test_energy_storage_device_factory(self):
        # valid input to buid an energy storage device
        for filename in valid_device_input:
            ptree=PropertyTree()
            ptree.parse_info(filename)
            device=EnergyStorageDevice(ptree)
        # invalid device will throw an exception
        ptree=PropertyTree()
        ptree.put_string('type', 'InvalidDevice')
        self.assertRaises(RuntimeError,EnergyStorageDevice,ptree)
    def test_energy_storage_device_not_copyable(self):
        device_database=PropertyTree()
        device_database.parse_xml('device.xml')
        device=EnergyStorageDevice(device_database.get_child('device'))
        from copy import copy,deepcopy
        self.assertRaises(RuntimeError,copy,device)
        self.assertRaises(RuntimeError,deepcopy,device)

if __name__ == '__main__':
    unittest.main()
