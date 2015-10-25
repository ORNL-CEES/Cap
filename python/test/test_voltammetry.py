from pycap import PropertyTree,EnergyStorageDevice,CyclicVoltammetry
from pycap import initialize_data
from numpy import array,testing,linspace
import unittest

ptree=PropertyTree()
ptree.put_string('type','SeriesRC')
ptree.put_double('series_resistance',40e-3)
ptree.put_double('capacitance',3)
device=EnergyStorageDevice(ptree)

class capRampTestCase(unittest.TestCase):
    def test_no_name(self):
        ptree=PropertyTree()   
        ptree.put_double('initial_voltage',0)
        ptree.put_double('final_voltage',0)
        ptree.put_double('scan_limit_1',1)
        ptree.put_double('scan_limit_2',0)
        ptree.put_double('step_size',0.1)
        ptree.put_double('scan_rate',1)
        ptree.put_int   ('cycles',1)
        cv=CyclicVoltammetry(ptree)
        try:
            cv.run(device)
        except:
            self.fail('calling run without data should not raise')
        data=initialize_data()
        cv.run(device,data)
        voltage=array([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0],dtype=float)
        time=linspace(0,2,21)
        try:
            testing.assert_array_almost_equal(data['voltage'],voltage)
            testing.assert_array_almost_equal(data['time'   ],time   )
        except AssertionError as e:
            print e
            self.fail()
        


if __name__ == '__main__':
    unittest.main()
