__all__=[
    'initialize_data',
    'report_data',
    'plot_data',
    'TimeEvolution',
    'EndCriterion',
]

from matplotlib import pyplot
from numpy import array,append

def initialize_data():
    return {
        'time'   :array([],dtype=float),
        'current':array([],dtype=float),
        'voltage':array([],dtype=float),
    }

def report_data(data,time,device):
    data['time'   ]=append(data['time'   ],time                )
    data['current']=append(data['current'],device.get_current())
    data['voltage']=append(data['voltage'],device.get_voltage())

def plot_data(data):
    time   =data['time'   ]
    current=data['current']
    voltage=data['voltage']
    plot_linewidth=3
    label_fontsize=20
    f,axarr=pyplot.subplots(2,sharex=True,figsize=(16,12))
    axarr[0].plot(time,1e3*current,lw=plot_linewidth)
    axarr[0].set_ylabel(r'$\mathrm{Current\ [mA]}$',fontsize=label_fontsize)
    axarr[1].plot(time,voltage,lw=plot_linewidth)
    axarr[1].set_ylabel(r'$\mathrm{Voltage\  [V]}$',fontsize=label_fontsize)
    axarr[1].set_xlabel(r'$\mathrm{Time\     [s]}$',fontsize=label_fontsize)
    pyplot.show()

class TimeEvolution:
    def factory(ptree):
        mode=ptree.get_string('type')
    
        if   mode=='constant_voltage' or mode=='potentiostatic':
            constant_voltage=ptree.get_double('voltage')
            def evolve_one_time_step_constant_voltage(device,time_step):
                device.evolve_one_time_step_constant_voltage(time_step,constant_voltage)
            return evolve_one_time_step_constant_voltage

        elif mode=='constant_current' or mode=='galvanostatic':
            constant_current=ptree.get_double('current')
            def evolve_one_time_step_constant_current(device,time_step):
                device.evolve_one_time_step_constant_current(time_step,constant_current)
            return evolve_one_time_step_constant_current

        elif mode=='constant_power'  :
            constant_power  =ptree.get_double('power'  )
            def evolve_one_time_step_constant_power  (device,time_step):
                device.evolve_one_time_step_constant_power  (time_step,constant_power  )
            return evolve_one_time_step_constant_power

        elif mode=='constant_load'   :
            constant_load   =ptree.get_double('load'   )
            def evolve_one_time_step_constant_load   (device,time_step):
                device.evolve_one_time_step_constant_load   (time_step,constant_load   )
            return evolve_one_time_step_constant_load

        elif mode=='hold'            :
            def evolve_one_time_step_hold            (device,time_step):
                device.evolve_one_time_step_constant_voltage(time_step,device.get_voltage())
            return evolve_one_time_step_hold
        
        elif mode=='rest'            :
            def evolve_one_time_step_rest            (device,time_step):
                device.evolve_one_time_step_constant_current(time_step,0.0)
            return evolve_one_time_step_rest

        else:
            raise NameError("invalid TimeEvolution mode '"+mode+"'")
        
    factory=staticmethod(factory)



from operator import le,ge

class EndCriterion:
    def check(self,time,device):
        raise NotImplementedError
    def reset(self,time,device):
        raise NotImplementedError
    def factory(ptree):
        type=ptree.get_string('end_criterion')
        if   type=='time':
            return TimeLimit(ptree)
        elif type=='voltage_greater_than':
            return VoltageLimit(ptree,ge)
        elif type=='voltage_less_than'   :
            return VoltageLimit(ptree,le)
        elif type=='current_greater_than':
            return CurrentLimit(ptree,ge)
        elif type=='current_less_than'   :
            return CurrentLimit(ptree,le)
        else:
            raise NameError("invalid EndCriterion type '"+type+"'")

    factory=staticmethod(factory)

class TimeLimit(EndCriterion):
    def __init__(self,ptree):
        self.duration=ptree.get_double('duration')
    def check(self,time,device):
        return time-self.tick>=self.duration
    def reset(self,time,device):
        self.tick=time

class VoltageLimit(EndCriterion):
    def __init__(self,ptree,compare):
        self.voltage_limit=ptree.get_double('voltage_limit')
        self.compare=compare
    def check(self,time,device):
        return self.compare(device.get_voltage(),self.voltage_limit)
    def reset(self,time,device):
        pass

class CurrentLimit(EndCriterion):
    def __init__(self,ptree,compare):
        self.current_limit=ptree.get_double('current_limit')
        self.compare=compare
    def check(self,time,device):
        return not self.compare(device.get_current(),self.current_limit)
    def reset(self,time,device):
        pass
