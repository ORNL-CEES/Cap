from pycap import PropertyTree, EnergyStorageDevice
from pycap import Charge
from pycap import initialize_data
from mpi4py import MPI
import unittest

comm = MPI.COMM_WORLD
filename = 'series_rc.info'
ptree = PropertyTree()
ptree.parse_info(filename)
device = EnergyStorageDevice(comm, ptree)


class capChargeTestCase(unittest.TestCase):
    def test_charge_constant_current(self):
        ptree = PropertyTree()
        ptree.put_string('charge_mode', 'constant_current')
        ptree.put_double('charge_current', 10e-3)
        ptree.put_string('charge_stop_at_1', 'voltage_greater_than')
        ptree.put_double('charge_voltage_limit', 1.4)
        ptree.put_double('time_step', 0.2)
        charge = Charge(ptree)
        data = initialize_data()
        charge.run(device, data)
        self.assertAlmostEqual(data['current'][0], 10e-3)
        self.assertAlmostEqual(data['current'][-1], 10e-3)
        self.assertGreaterEqual(data['voltage'][-1], 1.4)
        self.assertAlmostEqual(data['time'][1]-data['time'][0], 0.2)

    def test_charge_constant_voltage(self):
        ptree = PropertyTree()
        ptree.put_string('charge_mode', 'constant_voltage')
        ptree.put_double('charge_voltage', 1.4)
        ptree.put_string('charge_stop_at_1', 'current_less_than')
        ptree.put_double('charge_current_limit', 1e-6)
        ptree.put_string('charge_stop_at_2', 'time')
        ptree.put_double('charge_max_duration', 60)
        ptree.put_double('time_step', 0.2)
        charge = Charge(ptree)
        data = initialize_data()
        charge.run(device, data)
        self.assertTrue(data['time'][-1] >= 60 or
                        abs(data['current'][-1]) <= 1e-6)
        self.assertAlmostEqual(data['voltage'][-1], 1.4)


if __name__ == '__main__':
    unittest.main()
