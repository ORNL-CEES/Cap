from pycap import PropertyTree,EnergyStorageDevice
from numpy import array,pi,linalg,inf
import unittest

class capEnergyStorageDeviceWrappersTestCase(unittest.TestCase):
    def test_energy_storage_device_build_database(self):
        device_database=PropertyTree()
        # build a series rc equivalent circuit database
        device_database.put_string('type','SeriesRC')
        device_database.put_double('capacitance',3.0)
        device_database.put_double('series_resistance',100.0)
        device=EnergyStorageDevice(device_database)
        # change it to a parallel rc database
        device_database.put_string('type','ParallelRC')
        # database is incomplete at this time so building a device will throw
        self.assertRaises(RuntimeError,EnergyStorageDevice,device_database)
        # add the missing information to the database
        device_database.put_double('parallel_resistance',4000.0)
        # build the new energy device
        device=EnergyStorageDevice(device_database)
    def test_energy_storage_device_not_copyable(self):
        device_database=PropertyTree()
        device_database.parse_xml('device.xml')
        device=EnergyStorageDevice(device_database.get_child('device'))
        from copy import copy,deepcopy
        self.assertRaises(RuntimeError,copy,device)
        self.assertRaises(RuntimeError,deepcopy,device)

if __name__ == '__main__':
    unittest.main()
