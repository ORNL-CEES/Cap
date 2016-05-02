# Copyright (c) 2016, the Cap authors.
#
# This file is subject to the Modified BSD License and may not be distributed
# without copyright and license information. Please refer to the file LICENSE
# for the text and further information on this license.

from pycap import PropertyTree, EnergyStorageDevice
from mpi4py import MPI
import unittest

valid_device_input = [
    "series_rc.info",
    "parallel_rc.info",
    "super_capacitor.info",
]


class capEnergyStorageDeviceWrappersTestCase(unittest.TestCase):

    def test_energy_storage_device_factory(self):
        # valid input to buid an energy storage device
        for filename in valid_device_input:
            ptree = PropertyTree()
            ptree.parse_info(filename)
            EnergyStorageDevice(ptree)
            EnergyStorageDevice(ptree, comm=MPI.COMM_WORLD)
        # invalid device will throw an exception
        ptree = PropertyTree()
        ptree.put_string('type', 'InvalidDevice')
        self.assertRaises(RuntimeError, EnergyStorageDevice, ptree)

    def test_energy_storage_device_not_copyable(self):
        ptree = PropertyTree()
        ptree.parse_info('series_rc.info')
        device = EnergyStorageDevice(ptree)
        from copy import copy, deepcopy
        self.assertRaises(RuntimeError, copy, device)
        self.assertRaises(RuntimeError, deepcopy, device)

    def test_inspect_device(self):
        ptree = PropertyTree()
        ptree.parse_info('super_capacitor.info')
        device = EnergyStorageDevice(ptree, comm=MPI.COMM_WORLD)
        # method inspect() takes no argument and returns a dictionary
        data = device.inspect()
        self.assertTrue(isinstance(data, dict))
        for key in [
            'anode_electrode_interfacial_surface_area',
            'anode_electrode_mass_of_active_material',
            'anode_electrode_double_layer_capacitance',
            'anode_electrode_thickness',
            'cathode_electrode_interfacial_surface_area',
            'cathode_electrode_mass_of_active_material',
            'cathode_electrode_double_layer_capacitance',
            'cathode_electrode_thickness',
            'geometric_area',
        ]:
            self.assertTrue(key in data)
        print(data)

    def test_sanity(self):
        # valid input to buid an energy storage device
        for filename in valid_device_input:
            ptree = PropertyTree()
            ptree.parse_info(filename)
            device = EnergyStorageDevice(ptree)
            # these are basic sanity check to ensure that the device responds
            # in an appropriate manner when operating conditions are imposed
            dt = 0.1  # time_step in seconds
            I = 5e-3  # current in amperes
            device.evolve_one_time_step_constant_current(dt, I)
            self.assertAlmostEqual(device.get_current(), I)
            U = 1.1  # voltage in volts
            device.evolve_one_time_step_constant_voltage(dt, U)
            self.assertAlmostEqual(device.get_voltage(), U)

if __name__ == '__main__':
    unittest.main()
