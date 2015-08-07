import pycap
import matplotlib.pyplot
device=pycap.EnergyStorageDevice("super_capacitor.xml")
current=[]
voltage=[]
time=[]
time_step=0.01
for i in range(100):
    device.evolve_one_time_step_constant_voltage(time_step, 2.0)
    time.append(i*time_step)
    current.append(pycap.get_current(device))
    voltage.append(pycap.get_voltage(device))

#matplotlib.pyplot.plot(time,current)
#matplotlib.pyplot.show()
