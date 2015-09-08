import pycap
import numpy
import unittest


def evolve_constant_voltage(device,steps,time_step,constant_voltage):
    current=[]
    voltage=[]
    time=[]
    for i in range(steps):
        device.evolve_one_time_step_constant_voltage(time_step, constant_voltage)
        time.append(i*time_step)
        current.append(pycap.get_current(device))
        voltage.append(pycap.get_voltage(device))
    data=numpy.ndarray(shape=(steps,3),dtype=float)
    data[:,0]=time
    data[:,1]=current
    data[:,2]=voltage
    return data

class constantChargeTestCase(unittest.TestCase):
    def test_is_five_prime(self):
        device_database=pycap.PropertyTree()
        device_database.parse_xml('device.xml')
        device=pycap.EnergyStorageDevice(device_database)
        steps=1000
        time_step=0.01
        constant_voltage=2.0
        data=evolve_constant_voltage(device,steps,time_step,constant_voltage)
        self.assertEqual(data.shape,(steps,3))
        self.assertAlmostEqual(not pycap.get_current(device),0.0)
        self.assertAlmostEqual(pycap.get_voltage(device),2.0)

class boostPropertyTreePythonWrappersTestCase(unittest.TestCase):
    def test_get_array(self):
        ptree=pycap.PropertyTree()
        ptree.put_string('array_double','3.14,1.41')
        array_double=ptree.get_array_double('array_double')
        self.assertEqual(array_double,[3.14,1.41])
        ptree.put_string('array_int','1,2,3')
        array_int=ptree.get_array_int('array_int')
        self.assertEqual(array_int,[1,2,3])
        
    def test_property_tree(self):
        # ptree as container to store int, double, string, and bool
        ptree=pycap.PropertyTree()
        ptree.put_int('dim',3)
        self.assertEqual(ptree.get_int('dim'),3)
        ptree.put_double('path.to.pi',3.14)
        self.assertEqual(ptree.get_double('path.to.pi'),3.14)
        ptree.put_string('good.news','it works')
        self.assertEqual(ptree.get_string('good.news'),'it works')
        ptree.put_bool('is.that.a.good.idea',False)                   
        self.assertEqual(ptree.get_bool('is.that.a.good.idea'),False)
        # child nodes on nested trees
        ptree.put_string('child.name','clement')
        ptree.put_int('child.age',-2)
        child=ptree.get_child('child')
        self.assertEqual(child.get_string('name'),'clement')
        self.assertEqual(child.get_int('age'),-2)
        # property tree will throw if the specified path does not exist
        def throw_exception_bad_path(ptree):
            ptree.get_int('path.does.not.exist')
        self.assertRaises(RuntimeError, throw_exception_bad_path,ptree)
        # or if the translation fails
        def throw_exception_bad_data(ptree):
            ptree.put_string('some.path.to.a.string','not a double')
            ptree.get_double('some.path.to.a.string')
        self.assertRaises(RuntimeError, throw_exception_bad_data,ptree)

class capEnergyStorageDeviceWrappersTestCase(unittest.TestCase):
    def test_energy_storage_device_parse_xml_input(self):
        device_database=pycap.PropertyTree()
        device_database.parse_xml('device.xml')
        device=pycap.EnergyStorageDevice(device_database)
    def test_energy_storage_device_build_database(self):
        def build_energy_storage_device(device_database):
            device=pycap.EnergyStorageDevice(device_database)
            return device
        device_database=pycap.PropertyTree()
        # build a series rc equivalent circuit database
        device_database.put_string('device.type','SeriesRC')
        device_database.put_double('device.capacitance',3.0)
        device_database.put_double('device.series_resistance',100.0)
        # change it to a parallel rc database
        device_database.put_string('device.type','ParallelRC')
        # database is incomplete at this time so building a device will throw
        self.assertRaises(RuntimeError,build_energy_storage_device,device_database)
        def throw_exception_incomplete_database(device_database):
            device=pycap.EnergyStorageDevice(device_database)
        # add the missing information to the database
        device_database.put_double('device.parallel_resistance',4000.0)
        # build the new energy device
        device=pycap.EnergyStorageDevice(device_database)
        # build esd from get child
        device_database.put_string('child.device.type','SeriesRC')
        device_database.put_double('child.device.capacitance',3.0)
        device_database.put_double('child.device.series_resistance',100.0)
        child_database=device_database.get_child('child')
        device=pycap.EnergyStorageDevice(child_database)
    def test_measure_impedance(self):
        # build an energy storage device
        device_database=pycap.PropertyTree()
        device_database.parse_xml('device.xml')
        device=pycap.EnergyStorageDevice(device_database)
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
        device=pycap.EnergyStorageDevice(device_database)
        from copy import copy,deepcopy
        self.assertRaises(RuntimeError,copy,device)
        self.assertRaises(RuntimeError,deepcopy,device)

if __name__ == '__main__':
    unittest.main()

#matplotlib.pyplot.plot(time,current)
#matplotlib.pyplot.show()
