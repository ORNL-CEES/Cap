__all__=['measure_ragone_chart','plot_ragone_chart']

from numpy import trapz,count_nonzero,array,append,power
from matplotlib import pyplot
from ._pycap import PropertyTree
from .charge_discharge import Charge,Discharge
from .data_helpers import initialize_data,save_data

def plot_ragone_chart(data):
    power =data['power' ]
    energy=data['energy']
    plot_linewidth=3
    label_fontsize=20
    fig=pyplot.figure(figsize=(16,12))
    pyplot.plot(power,energy,'r-o',lw=plot_linewidth)
    pyplot.xscale('log')
    pyplot.yscale('log')
    pyplot.xlabel(r'$\mathrm{Power\  [W]}$',fontsize=label_fontsize)
    pyplot.ylabel(r'$\mathrm{Energy\ [J]}$',fontsize=label_fontsize)
    pyplot.show()



def run_discharge(device,ptree):

    data=initialize_data()
    
    # (re)charge the device
    initial_voltage=ptree.get_double('initial_voltage')

    charge_database=PropertyTree()
    charge_database.put_string('charge_mode'                        ,'constant_current')
    charge_database.put_double('charge_current'                     ,0.5)
    charge_database.put_string('charge_stop_at_1'                   ,'voltage_greater_than')
    charge_database.put_double('charge_voltage_limit'               ,initial_voltage)
    charge_database.put_bool  ('charge_voltage_finish'              ,True)
    charge_database.put_double('charge_voltage_finish_current_limit',1e-6)
    charge_database.put_double('charge_voltage_finish_max_time'     ,600)
    charge_database.put_double('charge_rest_time'                   ,0)
    charge_database.put_double('time_step'                          ,0.1)
    
    charge=Charge(charge_database)
    steps=charge.run(device,data)

    data['time']-=data['time'][-1]

    # discharge at constant power
    discharge_power=ptree.get_double('discharge_power')
    final_voltage  =ptree.get_double('final_voltage'  )
    time_step      =ptree.get_double('time_step'      )

    discharge_database=PropertyTree()
    discharge_database.put_string('discharge_mode',        'constant_power')
    discharge_database.put_double('discharge_power'        ,discharge_power)
    discharge_database.put_string('discharge_stop_at_1'    ,'voltage_less_than')
    discharge_database.put_double('discharge_voltage_limit',final_voltage)
    discharge_database.put_double('discharge_rest_time'    ,10*time_step)
    discharge_database.put_double('time_step'              ,time_step)
    
    discharge=Discharge(discharge_database)
    steps=discharge.run(device,data)
    
    return data



def examine_discharge(data):
    time   =data['time'   ]
    current=data['current']
    voltage=data['voltage']
    power=current[:]*voltage[:]
    mask=time[:]<=0
    energy_in =trapz(power[mask],time[mask])
    mask=time[:]>=0
    energy_out=trapz(power[mask],time[mask])
    return [energy_in,energy_out]



def measure_ragone_chart(device,ptree,fout={}):
    discharge_power_lower_limit=ptree.get_double('discharge_power_lower_limit')
    discharge_power_upper_limit=ptree.get_double('discharge_power_upper_limit')
    steps_per_decade           =ptree.get_int   ('steps_per_decade'           )
    min_steps_per_discharge    =ptree.get_int   ('min_steps_per_discharge'    )
    max_steps_per_discharge    =ptree.get_int   ('max_steps_per_discharge'    )
    discharge_power=discharge_power_lower_limit
    ragone_chart_data={'energy':array([],dtype=float),'power':array([],dtype=float)}
    while discharge_power<=discharge_power_upper_limit:
        #print discharge_power
        ptree.put_double('discharge_power',discharge_power)
        try:
            # this loop control the number of time steps in the discharge
            for measurement in ['first','second']:
                data=run_discharge(device,ptree)
                if fout:
                    path='ragone_chart_data'
                    path+='/power='+str(discharge_power)+'W'
                    path+='/'+measurement
                    save_data(data,path,fout)
                energy_in,energy_out=examine_discharge(data)
                steps=count_nonzero(data['time']>0)
#                print '-- steps={0:d} energy in={1:.7f} out={2:.7f}'.format(steps,energy_in,energy_out)
                if steps>=min_steps_per_discharge:
                    break
                else:
                    ptree.put_double('time_step',data['time'][-1]/max_steps_per_discharge)
        except RuntimeError:
            print 'Failed to discharge at %f watt'%discharge_power
            break
        # TODO:
        ragone_chart_data['energy']=append(ragone_chart_data['energy'],-energy_out)
        ragone_chart_data['power' ]=append(ragone_chart_data['power' ],discharge_power)
        discharge_power*=power(10.0,1.0/steps_per_decade)
    return ragone_chart_data
