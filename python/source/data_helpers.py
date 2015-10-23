__all__=[
    'initialize_data',
    'report_data',
    'save_data',
    'plot_data',
    'open_file_in_write_mode',
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

def save_data(data,path,fout):
    fout[path+'/time'   ]=data['time'   ]
    fout[path+'/current']=data['current']
    fout[path+'/voltage']=data['voltage']

def plot_data(data):
    time   =data['time'   ]
    current=data['current']
    voltage=data['voltage']
    plot_linewidth=3
    label_fontsize=30
    tick_fontsize=20
    f,axarr=pyplot.subplots(2,sharex=True,figsize=(16,12))
    axarr[0].plot(time,current,lw=plot_linewidth)
    axarr[0].set_ylabel(r'$\mathrm{Current\  [A]}$',fontsize=label_fontsize)
#    axarr[0].plot(time,1e3*current,lw=plot_linewidth)
#    axarr[0].set_ylabel(r'$\mathrm{Current\ [mA]}$',fontsize=label_fontsize)
    axarr[0].get_yaxis().set_tick_params(labelsize=tick_fontsize)
    axarr[1].plot(time,voltage,lw=plot_linewidth)
    axarr[1].set_ylabel(r'$\mathrm{Voltage\  [V]}$',fontsize=label_fontsize)
    axarr[1].set_xlabel(r'$\mathrm{Time\     [s]}$',fontsize=label_fontsize)
    axarr[1].get_xaxis().set_tick_params(labelsize=tick_fontsize)
    axarr[1].get_yaxis().set_tick_params(labelsize=tick_fontsize)

from h5py import File
from sys import stdout,exit

def open_file_in_write_mode(filename):
    try:
        fout=File(filename,'w-')
    except IOError:
        print "file '%s' already exists..."%filename
        stdout.write('overwrite it? [Y/n] ')
        yes=set(['yes','y',''])
        no =set(['no' ,'n'   ])
        while True:
            answer=raw_input().lower()
            if answer in yes:
                fout=File(filename,'w')
                break
            elif answer in no:
                exit(0)
            else:
                stdout.write("Please respond with 'yes' or 'no'")
    return fout
