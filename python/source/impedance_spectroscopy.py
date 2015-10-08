__all__=['measure_impedance_spectrum','plot_nyquist','plot_bode',
    'fourier_analysis',
    'retrieve_impedance_spectrum']

from matplotlib import pyplot
from numpy import real,imag,log10,absolute,angle,array,append,power,sin,pi,sum,isclose,fft,mean,argsort
from warnings import warn
from .data_helpers import initialize_data,report_data,save_data

def _is_power_of_two(x):
    return (x!=0) and (x&(x-1)==0)

def plot_nyquist(data,figure=None,ls='r-s'):
    impedance=data['impedance']
    resistance=real(impedance)
    reactance =imag(impedance)
    plot_linewidth=3
    label_fontsize=30
    tick_fontsize=20
    if figure:
        pyplot.figure(figure.number)
    else :
        pyplot.figure(figsize=(14,14))
    pyplot.plot(resistance,-reactance,ls,lw=plot_linewidth)
    pyplot.axis('equal')
    pyplot.xlabel(r'$\mathrm{Resistance\ [\Omega]}$',fontsize=label_fontsize)
    pyplot.ylabel(r'$\mathrm{-Reactance\ [\Omega]}$',fontsize=label_fontsize)
    pyplot.gca().xaxis.set_tick_params(labelsize=tick_fontsize)
    pyplot.gca().yaxis.set_tick_params(labelsize=tick_fontsize)

def plot_bode(data):
    frequency=data['frequency']
    impedance=data['impedance']
    magnitude=20*log10(absolute(impedance))
    phase=angle(impedance,deg=True)
    label_fontsize=30
    tick_fontsize=20
    labelx=-0.05
    labely=0.5
    plot_linewidth=3
    f,axarr=pyplot.subplots(2,sharex=True,figsize=(16,12))
    axarr[0].plot(frequency,magnitude,'b-o',lw=plot_linewidth)
    axarr[0].set_xscale('log')
    axarr[0].set_ylabel(r'$\mathrm{Magnitude\ [dB]}$',fontsize=label_fontsize)
    axarr[0].get_yaxis().set_tick_params(labelsize=tick_fontsize)
    axarr[0].yaxis.set_label_coords(labelx,labely)
    axarr[1].plot(frequency,phase,'g-o',lw=plot_linewidth)
    axarr[1].set_xscale('log')
    axarr[1].set_ylabel(r'$\mathrm{Phase\ [Degrees]}$',fontsize=label_fontsize)
    axarr[1].set_xlabel(r'$\mathrm{Frequency\ [Hz]}$',fontsize=label_fontsize)
    axarr[1].get_yaxis().set_tick_params(labelsize=tick_fontsize)
    axarr[1].get_xaxis().set_tick_params(labelsize=tick_fontsize)
    axarr[1].yaxis.set_label_coords(labelx,labely)

def run_one_cycle(device,ptree):
    frequency=ptree.get_double('frequency')
    dc_voltage=ptree.get_double('dc_voltage')
    harmonics=array(ptree.get_array_int('harmonics'))
    ac_amplitudes=array(ptree.get_array_double('amplitudes'))
    phases=array(ptree.get_array_double('phases'))*pi/180
    steps_per_cycle=ptree.get_int('steps_per_cycle')
    cycles=ptree.get_int('cycles')
    time_step=1/(frequency*steps_per_cycle)
    time=0.0
    data=initialize_data()
    for cycle in range(cycles):
        for step in range(steps_per_cycle):
            time+=time_step
            excitation_signal=dc_voltage+sum(ac_amplitudes*sin(2*pi*harmonics*frequency*time+phases))
            device.evolve_one_time_step_changing_voltage(time_step,excitation_signal)
            report_data(data,time,device)
    return data

def measure_impedance_spectrum(device,ptree,fout=None,dummy=None):
    frequency_upper_limit=ptree.get_double('frequency_upper_limit')
    frequency_lower_limit=ptree.get_double('frequency_lower_limit')
    steps_per_decade     =ptree.get_int   ('steps_per_decade'     )
    eis_data={'frequency':array([],dtype=float),'impedance':array([],dtype=complex)}
    frequency=frequency_upper_limit
    while frequency>=frequency_lower_limit:
        #print frequency
        ptree.put_double('frequency',frequency)
        data=run_one_cycle(device,ptree)
        if fout:
            path='/eis_data/frequency='+str(frequency)+'Hz'
            save_data(data,path,fout)
        f,Z=fourier_analysis(data,ptree)
        eis_data['frequency']=append(eis_data['frequency'],f)
        eis_data['impedance']=append(eis_data['impedance'],Z)
        if dummy:
            dummy(eis_data)
        frequency/=power(10.0,1.0/steps_per_decade)
    return eis_data

def retrieve_impedance_spectrum(fin):
    path='eis_data'
    eis_data={'frequency':array([],dtype=float),'impedance':array([],dtype=complex)}
    for key in fin[path].keys():
        start=key.find('=')+1
        end=key.find('Hz')
        frequency=float(key[start:end])
        cycling_data=fin[path][key]
        f,Z=fourier_analysis(cycling_data)
        eis_data['frequency']=append(eis_data['frequency'],f)
        eis_data['impedance']=append(eis_data['impedance'],Z)
    sort=argsort(eis_data['frequency'])
    eis_data['frequency']=eis_data['frequency'][sort]
    eis_data['impedance']=eis_data['impedance'][sort]
    eis_data['frequency']=eis_data['frequency'][::-1]
    eis_data['impedance']=eis_data['impedance'][::-1]
    return eis_data

from .peak_detection import peakdet

def fourier_analysis(data,ptree=None):
    time   =data['time'   ]
    current=data['current']
    voltage=data['voltage']

    # inspect data
    n=len(time) # normalization factor for fft
    assert len(current)==n
    assert len(voltage)==n
    d=time[1]-time[0] # inverse of the sampling rate
    # check sampling spacing is the same everywhere
    for i in range(n-1):
        assert isclose(time[i+1]-time[i],d,atol=1e-10,rtol=1e-10)

    # truncate signals
    if ptree:
        steps_per_cycle=ptree.get_int('steps_per_cycle')
        cycles         =ptree.get_int('cycles'         )
        ignore_cycles  =ptree.get_int('ignore_cycles'  )
        assert cycles>ignore_cycles
        assert n==cycles*steps_per_cycle
        time   =time   [ignore_cycles*steps_per_cycle:]
        current=current[ignore_cycles*steps_per_cycle:]
        voltage=voltage[ignore_cycles*steps_per_cycle:]
    else:
        time   =time   [n/2:]
        current=current[n/2:]
        voltage=voltage[n/2:]

    n=len(time)
    assert len(current)==n
    assert len(voltage)==n

    if not _is_power_of_two(n):
        warn(
            "(cycles-ignore_cycles)*steps_per_cycles is not a "
            "power of 2 (most efficient for the fourier analysis)",
            RuntimeWarning)

    # perform the actual fourrier analaysis
    fft_current=fft.rfft(current)/n
    fft_voltage=fft.rfft(voltage)/n
    fft_frequency=fft.rfftfreq(n,d)

    # find the excited harmonics
    if ptree:
        harmonics=array(ptree.get_array_int('harmonics'))
        peak_indices=harmonics*(cycles-ignore_cycles)
    else:
        mx,mn=peakdet(absolute(fft_voltage),mean(absolute(fft_current)))
        peak_indices=int(mx[:,0])
        mx,mn=peakdet(absolute(fft_voltage),mean(absolute(fft_current)))
        assert peak_indices==mx[:,0]

    frequency=fft_frequency[peak_indices]
    impedance=fft_voltage[peak_indices]/fft_current[peak_indices]

    return [frequency,impedance]
