# Copyright (c) 2016, the Cap authors.
#
# This file is subject to the Modified BSD License and may not be distributed
# without copyright and license information. Please refer to the file LICENSE
# for the text and further information on this license.

from .time_evolution import TimeEvolution
from .end_criterion import EndCriterion
from .data_helpers import report_data

__all__ = ['Stage', 'MultiStage']


class Stage:

    def __init__(self, ptree):
        self.evolve_one_time_step = TimeEvolution.factory(ptree)
        self.end_criterion = EndCriterion.factory(ptree)
        self.time_step = ptree.get_double('time_step')

    def run(self, device, data=None):
        if data is None:
            data = {}
        if not data:
            time = 0.0
        elif len(data['time']) > 0:
            time = data['time'][-1]
        else:
            time = 0.0
        steps = 0
        self.end_criterion.reset(time, device)
        while not self.end_criterion.check(time + 0.01 * self.time_step, device):
            steps += 1
            time += self.time_step
            self.evolve_one_time_step(device, self.time_step)
            if data:
                report_data(data, time, device)

        return steps


class MultiStage(Stage):

    def __init__(self, ptree):
        self.stages = []
        for stage in range(ptree.get_int('stages')):
            child = ptree.get_child('stage_' + str(stage))
            try:
                child.get_double('time_step')
            except:
                time_step = ptree.get_double('time_step')
                child.put_double('time_step', time_step)
            self.stages.append(Stage(child))
        self.cycles = ptree.get_int('cycles')

    def run(self, device, data=None):
        if data is None:
            data = {}
        steps = 0
        for cycle in range(self.cycles):
            for stage in self.stages:
                steps += stage.run(device, data)

        return steps
