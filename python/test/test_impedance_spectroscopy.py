# Copyright (c) 2016, the Cap authors.
#
# This file is subject to the Modified BSD License and may not be distributed
# without copyright and license information. Please refer to the file LICENSE
# for the text and further information on this license. 

from pycap import PropertyTree, EnergyStorageDevice, Experiment,\
                  ECLabAsciiFile
from pycap import measure_impedance_spectrum, retrieve_impedance_spectrum,\
                  fourier_analysis, initialize_data
from numpy import inf, linalg, real, imag, pi, log10, absolute, angle
from numpy import array, testing
from warnings import catch_warnings, simplefilter
from h5py import File
import unittest

R = 50.0e-3   # ohm
R_L = 500.0   # ohm
C = 3.0       # farad

ptree = PropertyTree()
ptree.put_string('type', 'ElectrochemicalImpedanceSpectroscopy')
ptree.put_double('frequency_upper_limit', 1e+4)
ptree.put_double('frequency_lower_limit', 1e-6)
ptree.put_int('steps_per_decade', 3)
ptree.put_int('steps_per_cycle', 1024)
ptree.put_int('cycles', 2)
ptree.put_int('ignore_cycles', 1)
ptree.put_double('dc_voltage', 0)
ptree.put_string('harmonics', '3')
ptree.put_string('amplitudes', '5e-3')
ptree.put_string('phases', '0')
eis = Experiment(ptree)

def setup_expertiment():
    eis_database = PropertyTree()
    eis_database.put_double('frequency_upper_limit', 1e+4)
    eis_database.put_double('frequency_lower_limit', 1e-6)
    eis_database.put_int('steps_per_decade', 3)
    eis_database.put_int('steps_per_cycle', 1024)
    eis_database.put_int('cycles', 2)
    eis_database.put_int('ignore_cycles', 1)
    eis_database.put_double('dc_voltage', 0)
    eis_database.put_string('harmonics', ' 3')
    eis_database.put_string('amplitudes', ' 5e-3')
    eis_database.put_string('phases', '0')

    return eis_database


class capImpedanceSpectroscopyTestCase(unittest.TestCase):
    def test_fourier_analysis(self):
        ptree = PropertyTree()
        ptree.put_int('steps_per_cycle', 3)
        ptree.put_int('cycles', 1)
        ptree.put_int('ignore_cycles', 0)
        ptree.put_string('harmonics', '1')
        # uninitialized data
        data = {}
        self.assertRaises(KeyError, fourier_analysis, data, ptree)
        # empty data
        data = initialize_data()
        self.assertRaises(IndexError, fourier_analysis, data, ptree)
        # bad data
        data['time'] = array([1, 2, 3], dtype=float)
        data['current'] = array([4, 5, 6], dtype=float)
        data['voltage'] = array([7, 8], dtype=float)
        self.assertRaises(AssertionError, fourier_analysis, data, ptree)
        # poor data (size not a power of 2)
        data['voltage'] = array([7, 8, 9], dtype=float)
        with catch_warnings():
            simplefilter("error")
            self.assertRaises(RuntimeWarning, fourier_analysis, data, ptree)
        # data unchanged after analyze
        dummy = array([1, 2, 3, 4, 5, 6, 7, 8], dtype=float)
        data['time'] = dummy
        data['current'] = dummy
        data['voltage'] = dummy
        # ptree needs to be updated
        self.assertRaises(AssertionError, fourier_analysis, data, ptree)
        ptree.put_int('steps_per_cycle', 4)
        ptree.put_int('cycles', 2)
        ptree.put_int('ignore_cycles', 0)
        fourier_analysis(data, ptree)
        try:
            testing.assert_array_equal(data['time'], dummy)
            testing.assert_array_equal(data['current'], dummy)
            testing.assert_array_equal(data['voltage'], dummy)
        except AssertionError:
            self.fail('data should not be changed by the fourier analyzis')

    def test_retrieve_data(self):
        # TODO: no need to use a device here
        device_database = PropertyTree()
        device_database.put_string('type', 'SeriesRC')
        device_database.put_double('series_resistance', R)
        device_database.put_double('capacitance', C)
        device = EnergyStorageDevice(device_database)
        eis_database  = setup_expertiment()
        eis_database.put_int('steps_per_decade', 1)
        eis_database.put_int('steps_per_cycle', 64)
        eis_database.put_int('cycles', 2)
        eis_database.put_int('ignore_cycles', 1)
        fout = File('trash.hdf5', 'w')
        spectrum_data = measure_impedance_spectrum(device, eis_database, fout)
        fout.close()
        fin = File('trash.hdf5', 'r')
        retrieved_data = retrieve_impedance_spectrum(fin)
        fin.close()
        print(spectrum_data['impedance']-retrieved_data['impedance'])
        print(retrieved_data)
        self.assertEqual(linalg.norm(spectrum_data['frequency'] -
                                     retrieved_data['frequency'], inf), 0.0)
        # not sure why we don't get equality for the impedance
        self.assertLess(linalg.norm(spectrum_data['impedance'] -
                                    retrieved_data['impedance'], inf), 1e-10)

    def test_verification_with_equivalent_circuit(self):
        # setup equivalent circuit database
        device_database = PropertyTree()
        device_database.put_double('series_resistance', R)
        device_database.put_double('parallel_resistance', R_L)
        device_database.put_double('capacitance', C)
        # analytical solutions
        Z = {}
        Z['SeriesRC'] = lambda f : R+1/(1j*C*2*pi*f)
        Z['ParallelRC'] = lambda f : R+R_L/(1+1j*R_L*C*2*pi*f)
        for device_type in ['SeriesRC', 'ParallelRC']:
            # create a device
            device_database.put_string('type', device_type)
            device = EnergyStorageDevice(device_database)
            # setup experiment and measure
            eis.reset()
            eis.run(device)
            f = eis._data['frequency']
            Z_computed = eis._data['impedance']
            # compute the exact solution
            Z_exact = Z[device_type](f)
            # ensure the error is small
            max_phase_error_in_degree = linalg.norm(
                angle(Z_computed)*180/pi - angle(Z_exact)*180/pi,
                inf)
            max_magniture_error_in_decibel = linalg.norm(
                20*log10(absolute(Z_exact)) - 20*log10(absolute(Z_computed)),
                inf)
            print(device_type)
            print('-- max_phase_error_in_degree = {0}'.format(max_phase_error_in_degree))
            print('-- max_magniture_error_in_decibel = {0}'.format(max_magniture_error_in_decibel))
            self.assertLessEqual(max_phase_error_in_degree, 1)
            self.assertLessEqual(max_magniture_error_in_decibel, 0.2)

    def test_export_eclab_ascii_format(self):
            pass

if __name__ == '__main__':
    unittest.main()
