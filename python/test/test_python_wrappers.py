import pycap
import numpy
import unittest

class capEnergyStorageDeviceWrappersTestCase(unittest.TestCase):
    def test_energy_storage_device_build_database(self):
        device_database=pycap.PropertyTree()
        # build a series rc equivalent circuit database
        device_database.put_string('type','SeriesRC')
        device_database.put_double('capacitance',3.0)
        device_database.put_double('series_resistance',100.0)
        # change it to a parallel rc database
        device_database.put_string('type','ParallelRC')
        # database is incomplete at this time so building a device will throw
        self.assertRaises(RuntimeError,pycap.EnergyStorageDevice,device_database)
        # add the missing information to the database
        device_database.put_double('parallel_resistance',4000.0)
        # build the new energy device
        device=pycap.EnergyStorageDevice(device_database)
        # build esd from get child
        device_database.put_string('child.type','SeriesRC')
        device_database.put_double('child.capacitance',3.0)
        device_database.put_double('child.series_resistance',100.0)
        child_database=device_database.get_child('child')
        device=pycap.EnergyStorageDevice(child_database)
    def test_measure_impedance(self):
        # build an energy storage device
        device_database=pycap.PropertyTree()
        device_database.parse_xml('device.xml')
        device=pycap.EnergyStorageDevice(device_database.get_child('device'))
        # measure its impedance
        eis_database=pycap.PropertyTree()
        eis_database.parse_xml('eis.xml')
        data=pycap.ElectrochemicalImpedanceSpectroscopyData()
        data.impedance_spectroscopy(device,eis_database)
        frequency=numpy.array(data.get_frequency())
        impedance=numpy.array(data.get_complex_impedance())
        w=2*numpy.pi*frequency
        # the following assume the device is a parallel rc circuit
        self.assertEqual(device_database.get_string('device.type'),'ParallelRC')
        C=device_database.get_double('device.capacitance')
        Rs=device_database.get_double('device.series_resistance')
        Rl=device_database.get_double('device.parallel_resistance')
        # compute theoretical impedance
        Z=Rs+Rl/(1+1j*Rl*C*w)
        # check that the relative error is less than .01%
        error_norm=numpy.linalg.norm((impedance-Z)/impedance,ord=numpy.inf)
        percent_tolerance=1.0e-2
        self.assertLessEqual(error_norm,percent_tolerance)
    def test_energy_storage_device_not_copyable(self):
        device_database=pycap.PropertyTree()
        device_database.parse_xml('device.xml')
        device=pycap.EnergyStorageDevice(device_database.get_child('device'))
        from copy import copy,deepcopy
        self.assertRaises(RuntimeError,copy,device)
        self.assertRaises(RuntimeError,deepcopy,device)

if __name__ == '__main__':
    unittest.main()
