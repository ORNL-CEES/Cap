# Copyright (c) 2016, the Cap authors.
#
# This file is subject to the Modified BSD License and may not be distributed
# without copyright and license information. Please refer to the file LICENSE
# for the text and further information on this license.

from pycap import PropertyTree, EnergyStorageDevice
from mpi4py import MPI
import unittest
import os

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
        # check 3D supercapacitor
        ptree = PropertyTree()
        ptree.parse_info(filename)
        ptree.put_int('dim', 3)
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
        # method inspect() takes an optional argument and returns a dictionary
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

    def test_postprocessor_inspect(self):
        ptree = PropertyTree()
        ptree.parse_info('series_rc.info')
        device = EnergyStorageDevice(ptree)
        self.assertRaises(RuntimeError, device.inspect, 'postprocessor')

        ptree = PropertyTree()
        ptree.parse_info('super_capacitor.info')
       # ptree.put_int('dim', 3)
        device = EnergyStorageDevice(ptree)
        data = device.inspect('postprocessor')
        self.assertTrue(isinstance(data, dict))
        self.assertEqual(len(data), 0)

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

    def test_checkpoint_restart(self):
        # check rc devices
        for filename in valid_device_input[0:2]:
            ptree = PropertyTree()
            ptree.parse_info(filename)
            device = EnergyStorageDevice(ptree)
            dt = 0.1  # time_step in seconds
            I = 5e-3  # current in amperes
            device.evolve_one_time_step_constant_current(dt, I)
            device.save('rc_device.txt')
            new_device = EnergyStorageDevice(ptree)
            new_device.load('rc_device.txt')
            self.assertEqual(device.get_voltage(), new_device.get_voltage())
            self.assertEqual(device.get_current(), new_device.get_current())
        # check supercapacitor
        coarse_mesh = 'coarse_mesh.z'
        device_state = 'device.z'
        ptree = PropertyTree()
        ptree.parse_info('super_capacitor.info')
        # check that forgetting to set checkpoint raises an error during the
        # restart
        ptree.put_string('geometry.coarse_mesh_filename', coarse_mesh)
        device = EnergyStorageDevice(ptree)
        device.save(device_state)
        new_ptree = PropertyTree()
        new_ptree.parse_info('super_capacitor.info')
        new_ptree.put_string('geometry.type', 'restart')
        new_ptree.put_string('geometry.coarse_mesh_filename', coarse_mesh)
        new_device = EnergyStorageDevice(new_ptree)
        self.assertRaises(RuntimeError, new_device.load, device_state)
        # Turn on checkpoint
        ptree.put_bool('geometry.checkpoint', True)
        device = EnergyStorageDevice(ptree)
        dt = 0.1  # time_step in seconds
        I = 5e-3  # current in amperes
        device.evolve_one_time_step_constant_current(dt, I)
        device.save(device_state)
        new_device.load(device_state)
        self.assertEqual(device.get_voltage(), new_device.get_voltage())
        self.assertEqual(device.get_current(), new_device.get_current())
        # clean files
        os.remove(coarse_mesh)
        os.remove(device_state)


if __name__ == '__main__':
    unittest.main()
