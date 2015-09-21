from pycap import PropertyTree,EnergyStorageDevice,Stage,initialize_data
import unittest

device_database=PropertyTree()
device_database.parse_xml('device.xml')
device=EnergyStorageDevice(device_database)

class capStageTestCase(unittest.TestCase):
    def test_nothing(self):
        ptree=PropertyTree()
        ptree.put_string('type','constant_current')
        ptree.put_double('current',5e-3)
        ptree.put_string('end_criterion','time')
        ptree.put_double('duration',15.0)
        ptree.put_double('time_step',0.1)
        print 'this should throw'
        stage=Stage(ptree)
        data=initialize_data()
        steps=stage.run(device,data)
        self.assertEqual(data['time'][-1],15.0)
        self.assertEqual(steps,150)
        self.assertEqual(device.get_current(),5e-3)

if __name__ == '__main__':
    unittest.main()
