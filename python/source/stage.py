__all__=['stage','multi_stage']

from .time_evolution import TimeEvolution
from .end_criterion import EndCriterion
from .data_helpers import report_data

class Stage:
    def __init__(self,ptree):
        self.evolve_one_time_step=TimeEvolution.factory(ptree)
        self.end_criterion=EndCriterion.factory(ptree)
        self.time_step=ptree.get_double('time_step')
    def run(self,device,data={}):
        if not data:
            time=0.0
        elif len(data['time'])>0:
            time=data['time'][-1]
        else:
            time=0.0
        steps=0
        self.end_criterion.reset(time,device)
        while not self.end_criterion.check(time,device):
            steps+=1
            time+=self.time_step
            self.evolve_one_time_step(device,self.time_step)
            if data:
                report_data(data,time,device)
        return steps

class MultiStage:
    def __init__(self,ptree):
        self.stages=[]
        for stage in range(ptree.get_int('stages')):
            child=ptree.get_child('stage_'+str(stage))
            try:
                child.get_double('time_step')
            except:
                time_step=ptree.get_double('time_step')
                child.put_double('time_step',time_step)
            self.stages.append(Stage(child))
        self.cycles=ptree.get_int('cycles')
    def run(self,device,data={}):
        steps=0
        for cycle in range(self.cycles):
            for stage in self.stages:
                steps+=stage.run(device,data)
                print steps
