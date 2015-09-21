from pycap import PropertyTree,EnergyStorageDevice,EndCriterion
from numpy import nan as NaN
import unittest

device_database=PropertyTree()
device_database.parse_xml('device.xml')
device=EnergyStorageDevice(device_database)

class capEndCriterionTestCase(unittest.TestCase):
    def test_time_limit(self):
        ptree=PropertyTree()
        ptree.put_string('end_criterion','time')
        ptree.put_double('duration',15)
        time_limit=EndCriterion.factory(ptree)
        time_limit.reset(0.0,device)
        self.assertFalse(time_limit.check(2.0,device))
        self.assertTrue(time_limit.check(15.0,device))
        self.assertTrue(time_limit.check(60.0,device))
    def test_voltage_limit(self):
        ptree=PropertyTree()
        ptree.put_double('voltage_limit',1.7)
        # upper limit
        ptree.put_string('end_criterion','voltage_greater_than')
        voltage_limit=EndCriterion.factory(ptree)
        voltage_limit.reset(5.0,device)
        device.evolve_one_time_step_constant_voltage(0.2,1.3)
        self.assertFalse(voltage_limit.check(0.0,device))
        self.assertFalse(voltage_limit.check(60.0,device))
        device.evolve_one_time_step_constant_voltage(0.2,1.7)
        self.assertTrue(voltage_limit.check(45.0,device))
        device.evolve_one_time_step_constant_voltage(0.2,2.1)
        self.assertTrue(voltage_limit.check(45.0,device))
        # lower limit
        ptree.put_string('end_criterion','voltage_less_than')
        voltage_limit=EndCriterion.factory(ptree)
        voltage_limit.reset(0.0,device)
        device.evolve_one_time_step_constant_voltage(0.2,1.3)
        self.assertTrue(voltage_limit.check(0.0,device))
        device.evolve_one_time_step_constant_voltage(0.2,1.7)
        self.assertTrue(voltage_limit.check(45.0,device))
        device.evolve_one_time_step_constant_voltage(0.2,2.1)
        self.assertFalse(voltage_limit.check(45.0,device))
    def test_current_limit(self):
        # lower
        ptree=PropertyTree()
        ptree.put_string('end_criterion','current_less_than')
        ptree.put_double('current_limit',-5e-3)
        self.assertRaises(RuntimeError,EndCriterion.factory,ptree)
        ptree.put_double('current_limit',0.0)
        self.assertRaises(RuntimeError,EndCriterion.factory,ptree)
        ptree.put_double('current_limit',5e-3)
        current_limit=EndCriterion.factory(ptree)
        device.evolve_one_time_step_constant_current(5.0,0.0)
        self.assertTrue(current_limit.check(NaN,device))
        device.evolve_one_time_step_constant_current(5.0,0.002)
        self.assertTrue(current_limit.check(NaN,device))
        device.evolve_one_time_step_constant_current(5.0,-0.001)
        self.assertTrue(current_limit.check(180.0,device))
        device.evolve_one_time_step_constant_current(5.0,0.005)
        self.assertTrue(current_limit.check(180.0,device))
        device.evolve_one_time_step_constant_current(5.0,0.007)
        self.assertFalse(current_limit.check(180.0,device))
        device.evolve_one_time_step_constant_current(5.0,-15e3)
        self.assertFalse(current_limit.check(180.0,device))
        # upper
        ptree.put_string('end_criterion','current_greater_than')
        ptree.put_double('current_limit',-5e-3)
        self.assertRaises(RuntimeError,EndCriterion.factory,ptree)
        ptree.put_double('current_limit',0.0)
        self.assertRaises(RuntimeError,EndCriterion.factory,ptree)
        ptree.put_double('current_limit',5e-3)
        current_limit=EndCriterion.factory(ptree)
        device.evolve_one_time_step_constant_current(5.0,-1e-3)
        self.assertFalse(current_limit.check(NaN,device))
        device.evolve_one_time_step_constant_current(5.0,0.002)
        self.assertFalse(current_limit.check(NaN,device))
        device.evolve_one_time_step_constant_current(5.0,0.005)
        self.assertTrue(current_limit.check(NaN,device))
        device.evolve_one_time_step_constant_current(5.0,-0.2)
        self.assertTrue(current_limit.check(NaN,device))
        device.evolve_one_time_step_constant_current(5.0,3.0)
        self.assertTrue(current_limit.check(NaN,device))
    def test_invalid_end_criterion(self):
        ptree=PropertyTree()
        ptree.put_string('end_criterion','bad_name')
        self.assertRaises(RuntimeError,EndCriterion.factory,ptree)
    def test_constructor(self):
        self.assertRaises(RuntimeError,EndCriterion)
        self.assertRaises(TypeError,EndCriterion,PropertyTree())
    def test_overload(self):
        class BadCriterion(EndCriterion):
            def __init__(self):
                pass
        bad_criterion=BadCriterion()
        self.assertRaises(NotImplementedError,bad_criterion.reset,0.0,device)
        self.assertRaises(NotImplementedError,bad_criterion.check,NaN,device)

if __name__ == '__main__':
    unittest.main()
