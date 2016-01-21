from pycap import PropertyTree, EnergyStorageDevice, EndCriterion
from numpy import nan as NaN
from mpi4py import MPI
import unittest

comm = MPI.COMM_WORLD
filename = 'series_rc.info'
ptree = PropertyTree()
ptree.parse_info(filename)
device = EnergyStorageDevice(comm,ptree)


class capEndCriterionTestCase(unittest.TestCase):
    def test_time_limit(self):
        ptree=PropertyTree()
        ptree.put_string('end_criterion', 'time')
        ptree.put_double('duration', 15)
        time_limit = EndCriterion.factory(ptree)
        time_limit.reset(0.0, device)
        self.assertFalse(time_limit.check(2.0, device))
        self.assertTrue(time_limit.check(15.0, device))
        self.assertTrue(time_limit.check(60.0, device))

    def test_voltage_limit(self):
        ptree = PropertyTree()
        ptree.put_double('voltage_limit', 1.7)
        # upper limit
        ptree.put_string('end_criterion', 'voltage_greater_than')
        voltage_limit = EndCriterion.factory(ptree)
        voltage_limit.reset(5.0, device)
        device.evolve_one_time_step_constant_voltage(0.2, 1.3)
        self.assertFalse(voltage_limit.check(0.0, device))
        self.assertFalse(voltage_limit.check(60.0, device))
        device.evolve_one_time_step_constant_voltage(0.2, 1.7)
        self.assertTrue(voltage_limit.check(45.0, device))
        device.evolve_one_time_step_constant_voltage(0.2, 2.1)
        self.assertTrue(voltage_limit.check(45.0, device))
        # lower limit
        ptree.put_string('end_criterion', 'voltage_less_than')
        voltage_limit = EndCriterion.factory(ptree)
        voltage_limit.reset(0.0, device)
        device.evolve_one_time_step_constant_voltage(0.2, 1.3)
        self.assertTrue(voltage_limit.check(0.0, device))
        device.evolve_one_time_step_constant_voltage(0.2, 1.7)
        self.assertTrue(voltage_limit.check(45.0, device))
        device.evolve_one_time_step_constant_voltage(0.2, 2.1)
        self.assertFalse(voltage_limit.check(45.0, device))

    def test_current_limit(self):
        # lower
        ptree = PropertyTree()
        ptree.put_string('end_criterion', 'current_less_than')
        ptree.put_double('current_limit', -5e-3)
        self.assertRaises(RuntimeError, EndCriterion.factory, ptree)
        ptree.put_double('current_limit', 0.0)
        self.assertRaises(RuntimeError, EndCriterion.factory, ptree)
        ptree.put_double('current_limit', 5e-3)
        current_limit = EndCriterion.factory(ptree)
        device.evolve_one_time_step_constant_current(5.0, 0.0)
        self.assertTrue(current_limit.check(NaN, device))
        device.evolve_one_time_step_constant_current(5.0, 0.002)
        self.assertTrue(current_limit.check(NaN, device))
        device.evolve_one_time_step_constant_current(5.0, -0.001)
        self.assertTrue(current_limit.check(180.0, device))
        device.evolve_one_time_step_constant_current(5.0, 0.005)
        self.assertTrue(current_limit.check(180.0, device))
        device.evolve_one_time_step_constant_current(5.0, 0.007)
        self.assertFalse(current_limit.check(180.0, device))
        device.evolve_one_time_step_constant_current(5.0, -15e3)
        self.assertFalse(current_limit.check(180.0, device))
        # upper
        ptree.put_string('end_criterion', 'current_greater_than')
        ptree.put_double('current_limit', -5e-3)
        self.assertRaises(RuntimeError, EndCriterion.factory, ptree)
        ptree.put_double('current_limit', 0.0)
        self.assertRaises(RuntimeError, EndCriterion.factory, ptree)
        ptree.put_double('current_limit', 5e-3)
        current_limit = EndCriterion.factory(ptree)
        device.evolve_one_time_step_constant_current(5.0, -1e-3)
        self.assertFalse(current_limit.check(NaN, device))
        device.evolve_one_time_step_constant_current(5.0, 0.002)
        self.assertFalse(current_limit.check(NaN, device))
        device.evolve_one_time_step_constant_current(5.0, 0.005)
        self.assertTrue(current_limit.check(NaN, device))
        device.evolve_one_time_step_constant_current(5.0, -0.2)
        self.assertTrue(current_limit.check(NaN, device))
        device.evolve_one_time_step_constant_current(5.0, 3.0)
        self.assertTrue(current_limit.check(NaN, device))

    def test_compound_criterion(self):
        ptree = PropertyTree()
        ptree.put_string('end_criterion', 'compound')
        ptree.put_string('criterion_0.end_criterion', 'time')
        ptree.put_double('criterion_0.duration', 5.0)
        ptree.put_string('criterion_1.end_criterion', 'voltage_greater_than')
        ptree.put_double('criterion_1.voltage_limit', 2.0)
        # no default value for now
        self.assertRaises(RuntimeError, EndCriterion.factory, ptree)
        ptree.put_string('logical_operator', 'bad_operator')
        self.assertRaises(RuntimeError, EndCriterion.factory, ptree)
        ptree.put_string('logical_operator', 'or')
        compound_criterion = EndCriterion.factory(ptree)
        compound_criterion.reset(0.0, device)
        device.evolve_one_time_step_constant_voltage(0.1, 1.0)
        self.assertFalse(compound_criterion.check(3.0, device))
        self.assertTrue(compound_criterion.check(5.0, device))
        device.evolve_one_time_step_constant_voltage(0.1, 2.0)
        self.assertTrue(compound_criterion.check(3.0, device))
        self.assertTrue(compound_criterion.check(5.0, device))
        ptree.put_string('logical_operator', 'and')
        compound_criterion = EndCriterion.factory(ptree)
        compound_criterion.reset(0.0, device)
        device.evolve_one_time_step_constant_voltage(0.1, 1.0)
        self.assertFalse(compound_criterion.check(3.0, device))
        self.assertFalse(compound_criterion.check(5.0, device))
        device.evolve_one_time_step_constant_voltage(0.1, 2.0)
        self.assertFalse(compound_criterion.check(3.0, device))
        self.assertTrue(compound_criterion.check(5.0, device))
        ptree.put_string('logical_operator', 'xor')
        compound_criterion = EndCriterion.factory(ptree)
        compound_criterion.reset(0.0, device)
        device.evolve_one_time_step_constant_voltage(0.1, 1.0)
        self.assertFalse(compound_criterion.check(3.0, device))
        self.assertTrue(compound_criterion.check(5.0, device))
        device.evolve_one_time_step_constant_voltage(0.1, 2.0)
        self.assertTrue(compound_criterion.check(3.0, device))
        self.assertFalse(compound_criterion.check(5.0, device))

    def test_never_statisfied(self):
        ptree = PropertyTree()
        ptree.put_string('end_criterion', 'none')
        never_statisfied = EndCriterion.factory(ptree)
        never_statisfied.reset(0.0, device)
        self.assertFalse(never_statisfied.check(NaN, device))

    def test_always_statisfied(self):
        ptree = PropertyTree()
        ptree.put_string('end_criterion', 'skip')
        always_statisfied = EndCriterion.factory(ptree)
        always_statisfied.reset(0.0, device)
        self.assertTrue(always_statisfied.check(NaN, device))

    def test_invalid_end_criterion(self):
        ptree = PropertyTree()
        ptree.put_string('end_criterion', 'bad_name')
        self.assertRaises(RuntimeError, EndCriterion.factory, ptree)

    def test_constructor(self):
        self.assertRaises(TypeError, EndCriterion)
        self.assertRaises(RuntimeError, EndCriterion, PropertyTree())

    def test_overload(self):
        class BadCriterion(EndCriterion):
            def __init__(self):
                pass
        bad_criterion = BadCriterion()
        self.assertRaises(NotImplementedError, bad_criterion.reset, 0.0,
                          device)
        self.assertRaises(NotImplementedError, bad_criterion.check, NaN,
                          device)

if __name__ == '__main__':
    unittest.main()
