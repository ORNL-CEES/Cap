# Copyright (c) 2016, the Cap authors.
#
# This file is subject to the Modified BSD License and may not be distributed
# without copyright and license information. Please refer to the file LICENSE
# for the text and further information on this license. 

from operator import le, ge, or_, and_, xor

__all__ = ['EndCriterion']


class EndCriterion:
    def __init__(self, ptree):
        raise RuntimeError('Use EndCriterion.factory to construct')

    def check(self, time, device):
        raise NotImplementedError

    def reset(self, time, device):
        raise NotImplementedError

    def factory(ptree):
        type = ptree.get_string('end_criterion')
        if type == 'time':
            return TimeLimit(ptree)
        elif type == 'voltage_greater_than':
            return VoltageLimit(ptree, ge)
        elif type == 'voltage_less_than':
            return VoltageLimit(ptree, le)
        elif type == 'current_greater_than':
            return CurrentLimit(ptree, ge)
        elif type == 'current_less_than':
            return CurrentLimit(ptree, le)
        elif type == 'compound':
            op = ptree.get_string('logical_operator')
            if op == 'or':
                logical_operator = or_
            elif op == 'and':
                logical_operator = and_
            elif op == 'xor':
                logical_operator = xor
            else:
                raise RuntimeError("Invalid logical operator '" + op +
                                   "' in CompoundCriterion")
            return CompoundCriterion(
                EndCriterion.factory(ptree.get_child('criterion_0')),
                EndCriterion.factory(ptree.get_child('criterion_1')),
                logical_operator)
        elif type == 'none':
            return NeverSatisfied()
        elif type == 'skip':
            return AlwaysSatisfied()
        else:
            raise RuntimeError("invalid EndCriterion type '"+type+"'")

    factory = staticmethod(factory)


class TimeLimit(EndCriterion):
    def __init__(self, ptree):
        self.duration = ptree.get_double('duration')

    def check(self, time, device):
        return time-self.tick >= self.duration

    def reset(self, time, device):
        self.tick = time


class VoltageLimit(EndCriterion):
    def __init__(self, ptree, compare):
        self.voltage_limit = ptree.get_double('voltage_limit')
        self.compare = compare

    def check(self, time, device):
        return self.compare(device.get_voltage(), self.voltage_limit)

    def reset(self, time, device):
        pass


class CurrentLimit(EndCriterion):
    def __init__(self, ptree, compare):
        self.current_limit = ptree.get_double('current_limit')
        if self.current_limit <= 0.0:
            raise RuntimeError(
                "CurrentLimit end criterion check for absolute value of the "
                "current. 'current_limit' (="+str(self.current_limit)+") "
                "must be greater than zero."
            )
        self.compare = compare

    def check(self, time, device):
        return self.compare(abs(device.get_current()), self.current_limit)

    def reset(self, time, device):
        pass


class CompoundCriterion(EndCriterion):
    def __init__(self, a, b, op):
        self.criterion_0 = a
        self.criterion_1 = b
        self.logical_operator = op

    def check(self, time, device):
        return self.logical_operator(
            self.criterion_0.check(time, device),
            self.criterion_1.check(time, device))

    def reset(self, time, device):
        for end_criterion in [self.criterion_0, self.criterion_1]:
            end_criterion.reset(time, device)


class NeverSatisfied(EndCriterion):
    def __init__(self):
        pass

    def check(self, time, device):
        return False

    def reset(self, time, device):
        pass


class AlwaysSatisfied(EndCriterion):
    def __init__(self):
        pass

    def check(self, time, device):
        return True

    def reset(self, time, device):
        pass
