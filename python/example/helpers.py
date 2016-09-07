from pycap import Observer, ECLabAsciiFile
from IPython import display
from numpy import real, imag, absolute, angle
from matplotlib import pyplot
from sys import stdout, exit
from os import remove

class PrintColumns(Observer):
    def __new__(cls, *args, **kwargs):
        return object.__new__(PrintColumns)
    def __init__(self):
        self._template = u''
        for i in range(3):
            self._template += '{left}{0}:{format_spec}{right}{separator}'\
                .format(i, format_spec='{format_spec}',
                        left='{', right='}', separator='\t')
    def update(self, subject, *args, **kwargs):
        extra = '>20'
        print(self._template.format('freq/Hz',
                                    'Re(Z)/ohm',
                                    '-Im(Z)/ohm',
                                    format_spec=extra+"s"),
              file=stdout)
        n = subject._data['frequency'].size
        for i in range(n):
            f = subject._data['frequency'][i]
            Z = subject._data['impedance'][i]
            Y = 1.0 / Z
            place_holder = 255
            line = self._template.format(float(f),
                                         float(real(Z)),
                                         -float(imag(Z)),                                   
                                         format_spec=extra+'.7e')
            print(line, file=stdout)

class RefreshDisplay(Observer):
    def __new__(cls, *args, **kwargs):
        return object.__new__(RefreshDisplay)
    def update(self, subject, *args, **kwargs):
        display.clear_output(wait=True)
        display.display(pyplot.gcf())

def check_input(device, experiment):
    experiment._extra_data = device.inspect()
    dummy = ECLabAsciiFile('dummy')
    dummy.update(experiment)
    with open('dummy', 'r', encoding='latin-1') as fin:
        lines = fin.readlines()
        for line in lines[7:-1]:
            print(line.rstrip('\n'))
    remove('dummy')

    print('continue? [Y/n]')
    yes = set(['yes', 'y', ''])
    no = set(['no', 'n'])
    while True:
        answer = input().lower()
        if answer in yes:
            break
        elif answer in no:
            exit(0)
        else:
            print("Please respond with 'yes' or 'no'")