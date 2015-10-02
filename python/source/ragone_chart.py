__all__=['measure_performance','plot_ragone','retrieve_performance_data']

from numpy import trapz,count_nonzero,array,append,power,argsort
from matplotlib import pyplot
from ._pycap import PropertyTree
from .charge_discharge import Charge,Discharge
from .data_helpers import initialize_data,save_data

def plot_ragone(data,figure=None,ls='r-o'):
    power =data['power' ]
    energy=data['energy']
    plot_linewidth=3
    label_fontsize=30
    tick_fontsize=20
    if figure:
        pyplot.figure(figure.number)
    else :
        pyplot.figure(figsize=(16,12))
    pyplot.plot(power,energy,ls,lw=plot_linewidth)
    pyplot.xscale('log')
    pyplot.yscale('log')
    pyplot.xlabel(r'$\mathrm{Power\  [W]}$',fontsize=label_fontsize)
    pyplot.ylabel(r'$\mathrm{Energy\ [J]}$',fontsize=label_fontsize)
    pyplot.gca().get_xaxis().set_tick_params(labelsize=tick_fontsize)
    pyplot.gca().get_yaxis().set_tick_params(labelsize=tick_fontsize)

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



def measure_performance(device,ptree,fout=None,dummy=None):
    discharge_power_lower_limit=ptree.get_double('discharge_power_lower_limit')
    discharge_power_upper_limit=ptree.get_double('discharge_power_upper_limit')
    steps_per_decade           =ptree.get_int   ('steps_per_decade'           )
    min_steps_per_discharge    =ptree.get_int   ('min_steps_per_discharge'    )
    max_steps_per_discharge    =ptree.get_int   ('max_steps_per_discharge'    )
    discharge_power=discharge_power_lower_limit
    performance_data={'energy':array([],dtype=float),'power':array([],dtype=float)}
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
        performance_data['energy']=append(performance_data['energy'],-energy_out)
        performance_data['power' ]=append(performance_data['power' ],discharge_power)
        if dummy:
            dummy(performance_data)
        discharge_power*=power(10.0,1.0/steps_per_decade)
    return performance_data



def retrieve_performance_data(fin):
    path='ragone_chart_data'
    performance_data={'power':array([],dtype=float),'energy':array([],dtype=float)}
    for key in fin[path].keys():
        start=key.find('=')+1
        end=key.find('W')
        power=float(key[start:end])
        for measurement in ['second','first']:
            if measurement in fin[path][key].keys():
                discharge_data=fin[path][key][measurement]
                (energy_in,energy_out)=examine_discharge(discharge_data)
                performance_data['power' ]=append(performance_data['power' ],power)
                performance_data['energy']=append(performance_data['energy'],-energy_out)
                break
    sort=argsort(performance_data['power'])
    performance_data['power' ]=performance_data['power' ][sort]
    performance_data['energy']=performance_data['energy'][sort]
    return performance_data
