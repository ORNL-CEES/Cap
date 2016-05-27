# Copyright (c) 2016, the Cap authors.
#
# This file is subject to the Modified BSD License and may not be distributed
# without copyright and license information. Please refer to the file LICENSE
# for the text and further information on this license.

from .PyCap import PropertyTree
from .stage import MultiStage

__all__ = ['Charge', 'Discharge', 'CyclicChargeDischarge']


class Charge(MultiStage):

    def __init__(self, ptree):
        other = PropertyTree()
        other.put_int('cycles', 1)
        other.put_int('stages', 3)
        time_step = ptree.get_double('time_step')
        assert time_step > 0.0
        other.put_double('time_step', time_step)
        # charge
        # time evolution
        charge_mode = ptree.get_string('charge_mode')
        voltage_finish = ptree.get_bool_with_default_value('charge_voltage_finish',
                                                           False)
        other.put_string('stage_0.mode', charge_mode)
        if charge_mode in ['constant_current', 'galvanostatic']:
            charge_current = ptree.get_double('charge_current')
            other.put_double('stage_0.current', charge_current)
            assert charge_current > 0.0
        elif charge_mode in ['constant_voltage', 'potentiostatic']:
            charge_voltage = ptree.get_double('charge_voltage')
            other.put_double('stage_0.voltage', charge_voltage)
            if voltage_finish:
                raise RuntimeError(
                    'voltage finish does not make sense under potentiostatic conditions')
        elif charge_mode == 'constant_power':
            charge_power = ptree.get_double('charge_power')
            other.put_double('stage_0.power', charge_power)
            assert charge_power > 0.0
        else:
            raise RuntimeError("invalid charge mode '" + charge_mode + "'")
        # end criterion
        other.put_string('stage_0.end_criterion', 'compound')
        other.put_string('stage_0.logical_operator', 'or')
        charge_stop_at_1 = ptree.get_string('charge_stop_at_1')
        if charge_stop_at_1 == 'voltage_greater_than':
            other.put_string('stage_0.criterion_0.end_criterion',
                             charge_stop_at_1)
            charge_voltage_limit = ptree.get_double('charge_voltage_limit')
            other.put_double('stage_0.criterion_0.voltage_limit',
                             charge_voltage_limit)
        elif charge_stop_at_1 == 'current_less_than':
            other.put_string('stage_0.criterion_0.end_criterion',
                             charge_stop_at_1)
            charge_current_limit = ptree.get_double('charge_current_limit')
            other.put_double('stage_0.criterion_0.current_limit',
                             charge_current_limit)
        else:
            raise RuntimeError("invalid charge_stop_at_1 criterion '" +
                               charge_stop_at_1 + "'")
        charge_stop_at_2 = ptree.get_string_with_default_value('charge_stop_at_2',
                                                               'none')
        other.put_string('stage_0.criterion_1.end_criterion', charge_stop_at_2)
        if charge_stop_at_2 != 'none':
            charge_max_duration = ptree.get_double('charge_max_duration')
            other.put_double('stage_0.criterion_1.duration',
                             charge_max_duration)
        # voltage finish
        other.put_string('stage_1.mode', 'constant_voltage')
        other.put_double('stage_1.voltage', 0.0)
        if voltage_finish:
            other.put_double('stage_1.voltage', charge_voltage_limit)
            charge_voltage_finish_max_time = ptree.get_double(
                'charge_voltage_finish_max_time')
            charge_voltage_finish_current_limit = ptree.get_double(
                'charge_voltage_finish_current_limit')
            other.put_string('stage_1.end_criterion', 'compound')
            other.put_string('stage_1.logical_operator', 'or')
            other.put_string('stage_1.criterion_0.end_criterion',
                             'current_less_than')
            other.put_double('stage_1.criterion_0.current_limit',
                             charge_voltage_finish_current_limit)
            other.put_string('stage_1.criterion_1.end_criterion',
                             'time')
            other.put_double('stage_1.criterion_1.duration',
                             charge_voltage_finish_max_time)
        else:
            other.put_string('stage_1.end_criterion', 'skip')
        # rest at open circuit
        other.put_string('stage_2.mode', 'rest')
        try:
            other.put_string('stage_2.end_criterion', 'time')
            charge_rest_time = ptree.get_double('charge_rest_time')
            assert charge_rest_time >= 0.0
            other.put_double('stage_2.duration', charge_rest_time)
        except RuntimeError:
            other.put_string('stage_2.end_criterion', 'skip')
        MultiStage.__init__(self, other)


class Discharge(MultiStage):

    def __init__(self, ptree):
        other = PropertyTree()
        other.put_int('cycles', 1)
        other.put_int('stages', 2)
        time_step = ptree.get_double('time_step')
        other.put_double('time_step', time_step)
        assert time_step > 0.0
        # discharge
        discharge_mode = ptree.get_string('discharge_mode')
        other.put_string('stage_0.mode', discharge_mode)
        if discharge_mode in ['constant_current', 'galvanostatic']:
            discharge_current = ptree.get_double('discharge_current')
            other.put_double('stage_0.current', -discharge_current)
            assert discharge_current > 0.0
        elif discharge_mode in ['constant_voltage', 'potentiostatic']:
            discharge_voltage = ptree.get_double('discharge_voltage')
            other.put_double('stage_0.voltage', discharge_voltage)
        elif discharge_mode == 'constant_power':
            discharge_power = ptree.get_double('discharge_power')
            other.put_double('stage_0.power', -discharge_power)
            assert discharge_power > 0.0
        elif discharge_mode == 'constant_load':
            discharge_load = ptree.get_double('discharge_load')
            other.put_double('stage_0.load', discharge_load)
            assert discharge_load > 0.0
        else:
            raise RuntimeError("invalid discharge mode '" +
                               discharge_mode + "'")
        discharge_stop_at_1 = ptree.get_string('discharge_stop_at_1')
        if discharge_stop_at_1 == 'voltage_less_than':
            other.put_string('stage_0.end_criterion', 'compound')
            other.put_string('stage_0.logical_operator', 'or')
            other.put_string('stage_0.criterion_0.end_criterion',
                             discharge_stop_at_1)
            discharge_voltage_limit = ptree.get_double(
                'discharge_voltage_limit')
            other.put_double('stage_0.criterion_0.voltage_limit',
                             discharge_voltage_limit)
        else:
            raise RuntimeError("invalid discharge_stop_at_1 criterion '" +
                               discharge_stop_at_1 + "'")
        try:
            discharge_stop_at_2 = ptree.get_string('discharge_stop_at_2')
            other.put_string('stage_0.criterion_1.end_criterion',
                             discharge_stop_at_2)
        except RuntimeError:
            other.put_string('stage_0.criterion_1.end_criterion', 'none')
        # rest at open circuit
        other.put_string('stage_1.mode', 'rest')
        try:
            other.put_string('stage_1.end_criterion', 'time')
            discharge_rest_time = ptree.get_double('discharge_rest_time')
            assert discharge_rest_time >= 0.0
            other.put_double('stage_1.duration', discharge_rest_time)
        except RuntimeError:
            other.put_string('stage_1.end_criterion', 'skip')
        MultiStage.__init__(self, other)


class CyclicChargeDischarge(MultiStage):

    def __init__(self, ptree):
        start_with = ptree.get_string('start_with')
        # TODO
        ptree.put_int("stages", 0)
        MultiStage.__init__(self, ptree)
        if start_with == 'charge':
            self.stages = [Charge(ptree), Discharge(ptree)]
        elif start_with == 'discharge':
            self.stages = [Discharge(ptree), Charge(ptree)]
        else:
            raise RuntimeError("Invalid first step '" + start_with +
                               "' in CyclicChargeDischarge")
