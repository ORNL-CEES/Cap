from pycap import PropertyTree,EnergyStorageDevice
from pycap import measure_impedance_spectrum,retrieve_impedance_spectrum,fourier_analysis,initialize_data
from numpy import sqrt,log,inf,linalg,real,imag,pi,log10,absolute,angle
from numpy import array,testing
from warnings import catch_warnings,simplefilter
import unittest

R  =50.0e-3 # ohm
R_L=500.0   # ohm
C  =3.0     # farad

def setup_expertiment():
    eis_database=PropertyTree()
    eis_database.put_double('frequency_upper_limit',1e+4)
    eis_database.put_double('frequency_lower_limit',1e-6)
    eis_database.put_int   ('steps_per_decade',3)
    eis_database.put_int   ('steps_per_cycle',1024)
    eis_database.put_int   ('cycles',2)
    eis_database.put_int   ('ignore_cycles',1)
    eis_database.put_double('dc_voltage',0)
    eis_database.put_string('harmonics','3')
    eis_database.put_string('amplitudes','5e-3')
    eis_database.put_string('phases','0')
    return eis_database

class capImpedanceSpectroscopyTestCase(unittest.TestCase):
    def testFourierAnalysis(self):
        ptree=PropertyTree()
        ptree.put_int('steps_per_cycle',3)
        ptree.put_int('cycles',1)
        ptree.put_int('ignore_cycles',0)
        ptree.put_string('harmonics','1')
        # uninitialized data
        data={}
        self.assertRaises(KeyError,fourier_analysis,data,ptree)
        # empty data
        data=initialize_data()
        self.assertRaises(IndexError,fourier_analysis,data,ptree)
        # bad data
        data['time'   ]=array([1,2,3],dtype=float)
        data['current']=array([4,5,6],dtype=float)
        data['voltage']=array([7,8],dtype=float)
        self.assertRaises(AssertionError,fourier_analysis,data,ptree)
        # poor data (size not a power of 2)
        data['voltage']=array([7,8,9],dtype=float)
        with catch_warnings():
            simplefilter("error")
            self.assertRaises(RuntimeWarning,fourier_analysis,data,ptree)
        # data unchanged after analyze
        dummy=array([1,2,3,4,5,6,7,8],dtype=float)
        data['time'   ]=dummy
        data['current']=dummy
        data['voltage']=dummy
        # ptree needs to be updated
        self.assertRaises(AssertionError,fourier_analysis,data,ptree)
        ptree.put_int('steps_per_cycle',4)
        ptree.put_int('cycles',2)
        ptree.put_int('ignore_cycles',0)
        fourier_analysis(data,ptree)
        try:
            testing.assert_array_equal(data['time'   ],dummy)
            testing.assert_array_equal(data['current'],dummy)
            testing.assert_array_equal(data['voltage'],dummy)
        except AssertionError:
             self.fail('data should not be changed by the fourier analyzis')
    def testRetrieveData(self):
        try:
            from h5py import File
        except ImportError:
            print 'module h5py not found'
            return
        device_database=PropertyTree()
        device_database.put_string('type','SeriesRC')
        device_database.put_double('series_resistance',R)
        device_database.put_double('capacitance'      ,C)
        device=EnergyStorageDevice(device_database)
        eis_database=setup_expertiment()
        eis_database.put_int   ('steps_per_decade',1)
        eis_database.put_int   ('steps_per_cycle',64)
        eis_database.put_int   ('cycles',2)
        eis_database.put_int   ('ignore_cycles',1)
        fout=File('trash.hdf5','w')
        spectrum_data=measure_impedance_spectrum(device,eis_database,fout)
        fout.close()
        fin=File('trash.hdf5','r')
        retrieved_data=retrieve_impedance_spectrum(fin)
        fin.close()
        print spectrum_data['impedance']-retrieved_data['impedance']
        print retrieved_data
        self.assertEqual(linalg.norm(spectrum_data['frequency']-retrieved_data['frequency'],inf),0.0)
        # not sure why we don't get equality for the impedance
        self.assertLess(linalg.norm(spectrum_data['impedance']-retrieved_data['impedance'],inf),1e-12)
    def testSeriesRC(self):
        # make series RC equivalent circuit
        device_database=PropertyTree()
        device_database.put_string('type','SeriesRC')
        device_database.put_double('series_resistance',R)
        device_database.put_double('capacitance'      ,C)      
        device=EnergyStorageDevice(device_database)
        # setup experiment and measure 
        eis_database=setup_expertiment()
        spectrum_data=measure_impedance_spectrum(device,eis_database)
        # extract data
        f         =spectrum_data['frequency']
        Z_computed=spectrum_data['impedance']
        R_computed=real(Z_computed)
        X_computed=imag(Z_computed)
        M_computed=20*log10(absolute(Z_computed))
        P_computed=angle(Z_computed)*180/pi
        # compute the exact solution
        Z_exact=R+1/(1j*C*2*pi*f)
        R_exact=real(Z_exact)
        X_exact=imag(Z_exact)
        M_exact=20*log10(absolute(Z_exact))
        P_exact=angle(Z_exact)*180/pi
        # ensure the error is small
        max_phase_error_in_degree=linalg.norm(P_computed-P_exact,inf)
        max_magniture_error_in_decibel=linalg.norm(M_computed-M_exact,inf)
        print 'max_phase_error_in_degree =',max_phase_error_in_degree
        print 'max_magniture_error_in_decibel =',max_magniture_error_in_decibel
        self.assertLessEqual(max_phase_error_in_degree,1)
        self.assertLessEqual(max_magniture_error_in_decibel,0.2)
    def testParallelRC(self):
        # make parallel RC equivalent circuit
        device_database=PropertyTree()
        device_database.put_string('type','ParallelRC')
        device_database.put_double('series_resistance'  ,R  )
        device_database.put_double('parallel_resistance',R_L)
        device_database.put_double('capacitance'        ,C  )      
        device=EnergyStorageDevice(device_database)
        # setup experiment and measure 
        eis_database=setup_expertiment()
        spectrum_data=measure_impedance_spectrum(device,eis_database)
        # extract data
        f         =spectrum_data['frequency']
        Z_computed=spectrum_data['impedance']
        M_computed=20*log10(absolute(Z_computed))
        P_computed=angle(Z_computed)*180/pi
        # compute the exact solution
        Z_exact=R+R_L/(1+1j*R_L*C*2*pi*f)
        M_exact=20*log10(absolute(Z_exact))
        P_exact=angle(Z_exact)*180/pi
        # ensure the error is small
        max_phase_error_in_degree=linalg.norm(P_computed-P_exact,inf)
        max_magniture_error_in_decibel=linalg.norm(M_computed-M_exact,inf)
        print 'max_phase_error_in_degree =',max_phase_error_in_degree
        print 'max_magniture_error_in_decibel =',max_magniture_error_in_decibel
        self.assertLessEqual(max_phase_error_in_degree,1)
        self.assertLessEqual(max_magniture_error_in_decibel,0.2)

if __name__ == '__main__':
    unittest.main()
