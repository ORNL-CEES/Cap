# Copyright (c) 2016, the Cap authors.
#
# This file is subject to the Modified BSD License and may not be distributed
# without copyright and license information. Please refer to the file LICENSE
# for the text and further information on this license.

from numpy import trapz, count_nonzero, array, append, power, argsort
from matplotlib import pyplot
from copy import copy
from .PyCap import PropertyTree
from .charge_discharge import Charge, Discharge
from .data_helpers import initialize_data, save_data
from .observer_pattern import Experiment, Observer

__all__ = ['plot_ragone', 'retrieve_performance_data',
           'RagonePlot', 'RagoneAnalysis']


def plot_ragone(data, figure=None, ls='r-o'):
    power = data['power']
    energy = data['energy']
    plot_linewidth = 3
    label_fontsize = 30
    tick_fontsize = 20
    if figure:
        pyplot.figure(figure.number)
    else:
        pyplot.figure(figsize=(16, 12))
    pyplot.plot(power, energy, ls, lw=plot_linewidth)
    pyplot.xscale('log')
    pyplot.yscale('log')
    pyplot.xlabel(r'$\mathrm{Power\  [W]}$', fontsize=label_fontsize)
    pyplot.ylabel(r'$\mathrm{Energy\ [J]}$', fontsize=label_fontsize)
    pyplot.gca().get_xaxis().set_tick_params(labelsize=tick_fontsize)
    pyplot.gca().get_yaxis().set_tick_params(labelsize=tick_fontsize)


def run_discharge(device, ptree):
    data = initialize_data()

    # (re)charge the device
    initial_voltage = ptree.get_double('initial_voltage')

    charge_database = PropertyTree()
    charge_database.put_string('charge_mode', 'constant_current')
    charge_database.put_double('charge_current', 0.5)
    charge_database.put_string('charge_stop_at_1', 'voltage_greater_than')
    charge_database.put_double('charge_voltage_limit', initial_voltage)
    charge_database.put_bool('charge_voltage_finish', True)
    charge_database.put_double('charge_voltage_finish_current_limit', 1e-6)
    charge_database.put_double('charge_voltage_finish_max_time', 600)
    charge_database.put_double('charge_rest_time', 0)
    charge_database.put_double('time_step', 0.1)

    charge = Charge(charge_database)
    charge.run(device, data)

    data['time'] -= data['time'][-1]

    # discharge at constant power
    discharge_power = ptree.get_double('discharge_power')
    final_voltage = ptree.get_double('final_voltage')
    time_step = ptree.get_double('time_step')

    discharge_database = PropertyTree()
    discharge_database.put_string('discharge_mode', 'constant_power')
    discharge_database.put_double('discharge_power', discharge_power)
    discharge_database.put_string('discharge_stop_at_1', 'voltage_less_than')
    discharge_database.put_double('discharge_voltage_limit', final_voltage)
    discharge_database.put_double('discharge_rest_time', 10 * time_step)
    discharge_database.put_double('time_step', time_step)

    discharge = Discharge(discharge_database)
    discharge.run(device, data)

    return data


def examine_discharge(data):
    time = data['time']
    current = data['current']
    voltage = data['voltage']
    power = current[:] * voltage[:]
    mask = time[:] <= 0
    energy_in = trapz(power[mask], time[mask])
    mask = time[:] >= 0
    energy_out = trapz(power[mask], time[mask])

    return [energy_in, energy_out]


def retrieve_performance_data(fin):
    path = 'ragone_chart_data'
    performance_data = {'power': array([], dtype=float),
                        'energy': array([], dtype=float)}
    for key in fin[path].keys():
        start = key.find('=') + 1
        end = key.find('W')
        power = float(key[start:end])
        for measurement in ['second', 'first']:
            if measurement in fin[path][key].keys():
                discharge_data = fin[path][key][measurement]
                (energy_in, energy_out) = examine_discharge(discharge_data)
                performance_data['power'] = append(performance_data['power'],
                                                   power)
                performance_data['energy'] = append(performance_data['energy'],
                                                    -energy_out)
                break
    sort = argsort(performance_data['power'])
    performance_data['power'] = performance_data['power'][sort]
    performance_data['energy'] = performance_data['energy'][sort]

    return performance_data


class RagonePlot(Observer):
    '''Ragone plot.

    Plots the values of specific energy versus specific power. Both axes are
    logarithmic.
    '''
    def __new__(cls, *args, **kwargs):
        return object.__new__(RagonePlot)

    def __init__(self, filename=None):
        self._figure = pyplot.figure(figsize=(14, 14))
        self._filename = filename

    def update(self, subject, *args, **kwargs):
        plot_ragone(subject._data, figure=self._figure)
        if self._filename is not None:
            pyplot.savefig(self._filename, bbox_inches='tight')


Observer._builders['RagonePlot'] = RagonePlot

# TODO: Probably want to make that class more general to do discharges at
# constant current as well.
# For batteries plot voltage vs capacity


class RagoneAnalysis(Experiment):
    '''Record a Ragone plot

    Performs a series of discharge at various rates and neasures the energy.

    Attributes
    ----------
    _discharge_power_lower_limit : float
    _discharge_upper_lower_limit : float
    _steps_per_decade : int
    _ptree : PropertyTree
    _data : dict
        Stores power and energy as numpy.array(s) of floating point numbers.

    See Also
    --------
    RagonePlot
    '''
    def __new__(cls, *args, **kwargs):
        return object.__new__(RagoneAnalysis)

    def __init__(self, ptree):
        Experiment.__init__(self)
        self._discharge_power_lower_limit = ptree.get_double(
            'discharge_power_lower_limit')
        self._discharge_power_upper_limit = ptree.get_double(
            'discharge_power_upper_limit')
        self._steps_per_decade = ptree.get_int('steps_per_decade')
        self._min_steps_per_discharge = ptree.get_int(
            'min_steps_per_discharge')
        self._max_steps_per_discharge = ptree.get_int(
            'max_steps_per_discharge')
        self._time_step_initial_guess = ptree.get_double('time_step')
        self._ptree = copy(ptree)
        self.reset()

    def reset(self):
        self._ptree.put_double('time_step', self._time_step_initial_guess)
        self._data = {
            'energy': array([], dtype=float),
            'power': array([], dtype=float)
        }

    def run(self, device, fout=None):
        discharge_power = self._discharge_power_lower_limit
        while discharge_power <= self._discharge_power_upper_limit:
            # print discharge_power
            self._ptree.put_double('discharge_power', discharge_power)
            try:
                # this loop control the number of time steps in the discharge
                for measurement in ['first', 'second']:
                    data = run_discharge(device, self._ptree)
                    if fout:
                        path = 'ragone_chart_data'
                        path += '/power=' + str(discharge_power) + 'W'
                        path += '/' + measurement
                        save_data(data, path, fout)
                    energy_in, energy_out = examine_discharge(data)
                    steps = count_nonzero(data['time'] > 0)
                    if steps >= self._min_steps_per_discharge:
                        break
                    else:
                        time_step = data['time'][-1]
                        time_step /= self._max_steps_per_discharge
                        self._ptree.put_double('time_step', time_step)

            except RuntimeError:
                print('Failed to discharge at {0} watt'.format(
                    discharge_power))
                break
            self._data['energy'] = append(self._data['energy'],
                                          -energy_out)
            self._data['power'] = append(self._data['power'],
                                         discharge_power)
            discharge_power *= power(10.0, 1.0 / self._steps_per_decade)
            self.notify()


Experiment._builders['RagoneAnalysis'] = RagoneAnalysis
