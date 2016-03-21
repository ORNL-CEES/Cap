# Copyright (c) 2016, the Cap authors.
#
# This file is subject to the Modified BSD License and may not be distributed
# without copyright and license information. Please refer to the file LICENSE
# for the text and further information on this license. 

from pycap import PropertyTree, EnergyStorageDevice, Experiment,\
                  RagoneAnalysis, RagonePlot
from pycap import measure_performance, retrieve_performance_data
from numpy import sqrt, log, inf, linalg
from mpi4py import MPI
from h5py import File
import unittest

comm = MPI.COMM_WORLD

R = 50.0e-3  # ohm
R_L = 500.0  # ohm
C = 3.0      # farad
U_i = 2.7    # volt
U_f = 1.2    # volt


def setup_expertiment():
    ragone_database = PropertyTree()
    ragone_database.put_double('discharge_power_lower_limit', 1e-2)
    ragone_database.put_double('discharge_power_upper_limit', 1e+2)
    ragone_database.put_int('steps_per_decade', 5)
    ragone_database.put_double('initial_voltage', U_i)
    ragone_database.put_double('final_voltage', U_f)
    ragone_database.put_double('time_step', 15)
    ragone_database.put_int('min_steps_per_discharge', 2000)
    ragone_database.put_int('max_steps_per_discharge', 3000)

    return ragone_database


class RagoneAnalysisTestCase(unittest.TestCase):
    def test_retrieve_data(self):
        ptree = PropertyTree()
        ptree.put_string('type', 'SeriesRC')
        ptree.put_double('series_resistance', 50e-3)
        ptree.put_double('capacitance', 3)
        device = EnergyStorageDevice(ptree)

        ptree = PropertyTree()
        ptree.put_string('type', 'RagoneAnalysis')
        ptree.put_double('discharge_power_lower_limit', 1e-1)
        ptree.put_double('discharge_power_upper_limit', 1e+1)
        ptree.put_int('steps_per_decade', 1)
        ptree.put_double('initial_voltage', 2.1)
        ptree.put_double('final_voltage', 0.7)
        ptree.put_double('time_step', 1.5)
        ptree.put_int('min_steps_per_discharge', 20)
        ptree.put_int('max_steps_per_discharge', 30)
        ragone = Experiment(ptree)
       
        with File('trash.hdf5', 'w') as fout:
            ragone.run(device, fout)
        performance_data = ragone._data

        fin = File('trash.hdf5', 'r')
        retrieved_data = retrieve_performance_data(fin)
        fin.close()
        # a few digits are lost when power is converted to string
        self.assertLess(linalg.norm(performance_data['power'] -
                                    retrieved_data['power'], inf), 1e-12)
        self.assertEqual(linalg.norm(performance_data['energy'] -
                                     retrieved_data['energy'], inf), 0.0)

        # TODO: probably want to move this into its own test
        ragoneplot = RagonePlot("ragone.png")
        ragoneplot.update(ragone)

    def testSeriesRC(self):
        # setup experiment
        ptree = PropertyTree()
        ptree.put_double('discharge_power_lower_limit', 1e-2)
        ptree.put_double('discharge_power_upper_limit', 1e+2)
        ptree.put_int('steps_per_decade', 3)
        ptree.put_double('initial_voltage', U_i)
        ptree.put_double('final_voltage', U_f)
        ptree.put_double('time_step', 15)
        ptree.put_int('min_steps_per_discharge', 2000)
        ptree.put_int('max_steps_per_discharge', 3000)
        ragone = RagoneAnalysis(ptree)
        # make series RC equivalent circuit
        device_database = PropertyTree()
        device_database.put_string('type', 'SeriesRC')
        device_database.put_double('series_resistance', R)
        device_database.put_double('parallel_resistance', R_L)
        device_database.put_double('capacitance', C)
        device = EnergyStorageDevice(device_database)
        # setup experiment and measure
        ragone.run(device)
        performance_data = ragone._data
        # extract data
        E_computed = performance_data['energy']
        P = performance_data['power']
        # compute the exact solution
        U_0 = U_i/2+sqrt(U_i**2/4-R*P)
        E_exact = C/2*(-R*P*log(U_0**2/U_f**2)+U_0**2-U_f**2)
        # ensure the error is small
        max_percent_error = 100*linalg.norm((E_computed-E_exact)/E_computed,
                                            inf)
        self.assertLess(max_percent_error, 0.1)

    def testParallelRC(self):
        # make parallel RC equivalent circuit
        device_database = PropertyTree()
        device_database.put_string('type', 'ParallelRC')
        device_database.put_double('series_resistance', R)
        device_database.put_double('parallel_resistance', R_L)
        device_database.put_double('capacitance', C)
        device = EnergyStorageDevice(device_database, comm)
        # setup experiment and measure
        ragone_database = setup_expertiment()
        performance_data = measure_performance(device, ragone_database)
        # extract data
        E_computed = performance_data['energy']
        P = performance_data['power']
        # compute the exact solution
        U_0 = U_i/2+sqrt(U_i**2/4-R*P)
        tmp = (U_f**2/R_L+P*(1+R/R_L))/(U_0**2/R_L+P*(1+R/R_L))
        E_exact = C/2*(-R_L*P*log(tmp)-R*R_L/(R+R_L)*P*log(tmp*U_0**2/U_f**2))
        # ensure the error is small
        max_percent_error = 100*linalg.norm((E_computed-E_exact)/E_computed,
                                            inf)
        self.assertLess(max_percent_error, 0.1)

if __name__ == '__main__':
    unittest.main()
