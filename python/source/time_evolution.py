# Copyright (c) 2016, the Cap authors.
#
# This file is subject to the Modified BSD License and may not be distributed
# without copyright and license information. Please refer to the file LICENSE
# for the text and further information on this license.

__all__ = ['TimeEvolution']


class TimeEvolution:

    def __init__(self, ptree):
        raise RuntimeError('Use TimeEvolution.factory to construct')

    def factory(ptree):
        mode = ptree.get_string('mode')

        if mode == 'constant_voltage' or mode == 'potentiostatic':
            constant_voltage = ptree.get_double('voltage')

            def evolve_one_time_step_constant_voltage(device, time_step):
                device.evolve_one_time_step_constant_voltage(time_step,
                                                             constant_voltage)
            return evolve_one_time_step_constant_voltage

        elif mode == 'constant_current' or mode == 'galvanostatic':
            constant_current = ptree.get_double('current')

            def evolve_one_time_step_constant_current(device, time_step):
                device.evolve_one_time_step_constant_current(time_step,
                                                             constant_current)
            return evolve_one_time_step_constant_current

        elif mode == 'constant_power':
            constant_power = ptree.get_double('power')

            def evolve_one_time_step_constant_power(device, time_step):
                device.evolve_one_time_step_constant_power(time_step,
                                                           constant_power)
            return evolve_one_time_step_constant_power

        elif mode == 'constant_load':
            constant_load = ptree.get_double('load')

            def evolve_one_time_step_constant_load(device, time_step):
                device.evolve_one_time_step_constant_load(time_step,
                                                          constant_load)
            return evolve_one_time_step_constant_load

        elif mode == 'hold':
            def evolve_one_time_step_hold(device, time_step):
                device.evolve_one_time_step_constant_voltage(time_step,
                                                             device.get_voltage())
            return evolve_one_time_step_hold

        elif mode == 'rest':
            def evolve_one_time_step_rest(device, time_step):
                device.evolve_one_time_step_constant_current(time_step, 0.0)
            return evolve_one_time_step_rest

        else:
            raise RuntimeError("invalid TimeEvolution mode '" + mode + "'")

    factory = staticmethod(factory)
