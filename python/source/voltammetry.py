__all__=['CyclicVoltammetry','plot_cyclic_voltammogram']

from operator import lt,gt,sub,add
from matplotlib import pyplot
from .data_helpers import report_data

def plot_cyclic_voltammogram(data,figure=None,ls='r-o'):
    current=data['current']
    voltage=data['voltage']
    plot_linewidth=3
    label_fontsize=30
    tick_fontsize=20
    if figure:
        pyplot.figure(figure.number)
    else :
        pyplot.figure(figsize=(16,12))
    pyplot.figure(figsize=(16,12))
    pyplot.plot(voltage,current,ls,lw=plot_linewidth)
    pyplot.xlabel(r'$\mathrm{Voltage\ [V]}$',fontsize=label_fontsize)
    pyplot.ylabel(r'$\mathrm{Current\ [A]}$',fontsize=label_fontsize)
    pyplot.gca().get_xaxis().set_tick_params(labelsize=tick_fontsize)
    pyplot.gca().get_yaxis().set_tick_params(labelsize=tick_fontsize)

def ramp(device,data,voltage_limit,scan_rate,step_size):
    step=0
    initial_voltage=device.get_voltage()
    final_voltage  =voltage_limit
    if initial_voltage>final_voltage:
        compare=gt
        update=sub
    else:
        compare=lt
        update=add
    time_step=step_size/scan_rate
    if data:
        time=data['time'][-1]
    else:
        time=0
    voltage=initial_voltage
    while compare(voltage,final_voltage):
        step+=1
        voltage=update(voltage,step_size)
        time+=time_step
        device.evolve_one_time_step_changing_voltage(time_step,voltage)
        report_data(data,time,device)
    return step

class CyclicVoltammetry:
    def __init__(self,ptree):
        self.cycles         =ptree.get_int   ('cycles'         )
        self.scan_limit_1   =ptree.get_double('scan_limit_1'   )
        self.scan_limit_2   =ptree.get_double('scan_limit_2'   )
        self.initial_voltage=ptree.get_double('initial_voltage')
        self.final_voltage  =ptree.get_double('final_voltage'  )
        self.scan_rate      =ptree.get_double('scan_rate'      )
        self.step_size      =ptree.get_double('step_size'      )
    def run(self,device,data=None):
        device.evolve_one_time_step_changing_voltage(self.step_size/self.scan_rate,self.initial_voltage)
        if data:
            report_data(data,0.0,device)
        steps=0
        for cycle in range(self.cycles):
            steps+=ramp(device,data,self.scan_limit_1,self.scan_rate,self.step_size)
            steps+=ramp(device,data,self.scan_limit_2,self.scan_rate,self.step_size)    
        steps+=ramp(device,data,self.final_voltage,self.scan_rate,self.step_size)
        return steps
