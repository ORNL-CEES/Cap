from pycap import PropertyTree,EnergyStorageDevice
from pycap import ElectrochemicalImpedanceSpectroscopyData
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
# DEPRECATED
    def test_measure_impedance(self):
        # build an energy storage device
        device_database=PropertyTree()
        device_database.parse_xml('device.xml')
        device=EnergyStorageDevice(device_database.get_child('device'))
        # measure its impedance
        eis_database=PropertyTree()
        eis_database.parse_xml('eis.xml')
        data=ElectrochemicalImpedanceSpectroscopyData()
        data.impedance_spectroscopy(device,eis_database)
        f=array(data.get_frequency())
        Z_computed=array(data.get_complex_impedance())
        if not f:
            return
        # the following assume the device is a parallel rc circuit
        self.assertEqual(device_database.get_string('device.type'),'ParallelRC')
        C=device_database.get_double('device.capacitance')
        Rs=device_database.get_double('device.series_resistance')
        Rl=device_database.get_double('device.parallel_resistance')
        # compute theoretical impedance
        Z_exact=Rs+Rl/(1+1j*Rl*C*2*pi*f)
        # check that the relative error is less than .01%
        error_norm=linalg.norm((Z_conputed-Z_exact)/Z_computed,ord=inf)
        percent_tolerance=1.0e-2
        self.assertLessEqual(error_norm,percent_tolerance)

if __name__ == '__main__':
    unittest.main()
