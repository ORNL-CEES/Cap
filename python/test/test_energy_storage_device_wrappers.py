# Copyright (c) 2016, the Cap authors.
#
# This file is subject to the Modified BSD License and may not be distributed
# without copyright and license information. Please refer to the file LICENSE
# for the text and further information on this license. 

from pycap import PropertyTree, EnergyStorageDevice
from mpi4py import MPI
import unittest

comm = MPI.COMM_WORLD

valid_device_input = [
    "series_rc.info",
    "parallel_rc.info",
#    "super_capacitor.info",
    ]


class capEnergyStorageDeviceWrappersTestCase(unittest.TestCase):
    def test_energy_storage_device_factory(self):
        # valid input to buid an energy storage device
        for filename in valid_device_input:
            ptree = PropertyTree()
            ptree.parse_info(filename)
            EnergyStorageDevice(ptree, comm)
        # invalid device will throw an exception
        ptree = PropertyTree()
        ptree.put_string('type', 'InvalidDevice')
        self.assertRaises(RuntimeError, EnergyStorageDevice, ptree, comm)

    def test_energy_storage_device_not_copyable(self):
        device_database = PropertyTree()
        device_database.parse_xml('device.xml')
        device = EnergyStorageDevice(device_database.get_child('device'), comm)
        from copy import copy, deepcopy
        self.assertRaises(RuntimeError, copy, device)
        self.assertRaises(RuntimeError, deepcopy, device)

    def test_inspect_device(self):
        ptree = PropertyTree()
        ptree.parse_info('super_capacitor.info')
        device = EnergyStorageDevice(ptree, comm)
        # method inspect() takes no argument and returns a dictionary
        data = device.inspect()
        self.assertTrue(isinstance(data, dict))
        for key in [
            'surface_area',
            'mass',
            'capacitance',
            'electrode_thickness',
            'geometric_area',
            ]:
            self.assertTrue(key in data)
        print(data)

if __name__ == '__main__':
    unittest.main()
