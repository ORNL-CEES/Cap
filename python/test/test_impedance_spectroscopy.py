from pycap import PropertyTree,EnergyStorageDevice
from pycap import measure_impedance_spectrum
from numpy import sqrt,log,inf,linalg,real,imag,pi,log10,absolute,angle
import unittest

R  =50.0e-3 # ohm
R_L=500.0   # ohm
C  =3.0     # farad

def setup_expertiment():
    eis_database=PropertyTree()
    eis_database.put_double('frequency_upper_limit',1e+4)
    eis_database.put_double('frequency_lower_limit',1e-6)
    eis_database.put_int   ('steps_per_decade',10)
    eis_database.put_int   ('steps_per_cycle',1024)
    eis_database.put_double('dc_voltage',0)
    eis_database.put_string('harmonics','3')
    eis_database.put_string('amplitudes','5e-3')
    eis_database.put_string('phases','0')
    return eis_database

class capImpedanceSpectroscopyTestCase(unittest.TestCase):
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
