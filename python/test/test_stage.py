# Copyright (c) 2016, the Cap authors.
#
# This file is subject to the Modified BSD License and may not be distributed
# without copyright and license information. Please refer to the file LICENSE
# for the text and further information on this license. 

from pycap import PropertyTree, EnergyStorageDevice
from pycap import Stage, MultiStage
from pycap import initialize_data
from mpi4py import MPI
import unittest

comm = MPI.COMM_WORLD
filename = 'series_rc.info'
ptree = PropertyTree()
ptree.parse_info(filename)
device = EnergyStorageDevice(ptree, comm)


class capStageTestCase(unittest.TestCase):
    def test_constant_current_charge_for_given_time(self):
        ptree = PropertyTree()
        ptree.put_string('mode', 'constant_current')
        ptree.put_double('current', 5e-3)
        ptree.put_string('end_criterion', 'time')
        ptree.put_double('duration', 15.0)
        ptree.put_double('time_step', 0.1)
        stage = Stage(ptree)
        data = initialize_data()
        steps = stage.run(device, data)
        self.assertEqual(steps, 150)
        self.assertEqual(steps, len(data['time']))
        self.assertAlmostEqual(data['time'][-1], 15.0)
        self.assertAlmostEqual(data['current'][-1], 5e-3)

    def test_force_discharge(self):
        ptree = PropertyTree()
        ptree.put_string('mode', 'constant_voltage')
        ptree.put_double('voltage', 0.0)
        ptree.put_string('end_criterion', 'current_less_than')
        ptree.put_double('current_limit', 1e-5)
        ptree.put_double('time_step', 1.0)
        stage = Stage(ptree)
        data = initialize_data()
        steps = stage.run(device, data)
        self.assertGreaterEqual(steps, 1)
        self.assertEqual(steps, len(data['time']))
        self.assertAlmostEqual(data['voltage'][-1], 0.0)
        self.assertLessEqual(data['current'][-1], 1e-5)

    def test_time_steps(self):
        ptree = PropertyTree()
        ptree.put_int('stages', 2)
        ptree.put_int('cycles', 1)
        ptree.put_double('time_step', 1.0)
        ptree.put_string('stage_0.mode', 'hold')
        ptree.put_string('stage_0.end_criterion', 'time')
        ptree.put_double('stage_0.duration', 2.0)
        ptree.put_string('stage_1.mode', 'rest')
        ptree.put_string('stage_1.end_criterion', 'time')
        ptree.put_double('stage_1.duration', 1.0)
        ptree.put_double('stage_1.time_step', 0.1)
        multi = MultiStage(ptree)
        data = initialize_data()
        steps = multi.run(device, data)
        self.assertEqual(steps, 12)
        self.assertEqual(steps, len(data['time']))
        self.assertAlmostEqual(data['time'][5], 2.4)
        self.assertAlmostEqual(data['voltage'][0], data['voltage'][1])
        self.assertAlmostEqual(data['current'][3], 0.0)

if __name__ == '__main__':
    unittest.main()
