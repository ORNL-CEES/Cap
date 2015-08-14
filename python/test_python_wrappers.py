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
        device=pycap.EnergyStorageDevice("device.xml")
        steps=1000
        time_step=0.01
        constant_voltage=2.0
        data=evolve_constant_voltage(device,steps,time_step,constant_voltage)
        self.assertEqual(data.shape,(steps,3))
        self.assertAlmostEqual(not pycap.get_current(device),0.0)
        self.assertAlmostEqual(pycap.get_voltage(device),2.0)

if __name__ == '__main__':
    unittest.main()

#matplotlib.pyplot.plot(time,current)
#matplotlib.pyplot.show()
