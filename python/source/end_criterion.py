__all__=['EndCriterion']

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
            raise RuntimeError("invalid EndCriterion type '"+type+"'")

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
