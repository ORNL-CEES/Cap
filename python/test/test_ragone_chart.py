from pycap import PropertyTree,EnergyStorageDevice
from pycap import measure_ragone_chart
from numpy import sqrt,log,inf,linalg
import unittest

R  =50.0e-3 # ohm
R_L=500.0   # ohm
C  =3.0     # farad
U_i=2.7     # volt
U_f=1.2     # volt

def setup_expertiment():
    ragone_database=PropertyTree()
    ragone_database.put_double('discharge_power_lower_limit',1e-2)
    ragone_database.put_double('discharge_power_upper_limit',1e+2)
    ragone_database.put_int   ('steps_per_decade'           ,5)
    ragone_database.put_double('initial_voltage'            ,U_i)
    ragone_database.put_double('final_voltage'              ,U_f)
    ragone_database.put_double('time_step'                  ,15)
    ragone_database.put_int   ('min_steps_per_discharge'    ,2000)
    ragone_database.put_int   ('max_steps_per_discharge'    ,3000)
    return ragone_database

class capRagoneChartTestCase(unittest.TestCase):
    def testSeriesRC(self):
        # make series RC equivalent circuit
        device_database=PropertyTree()
        device_database.put_string('type','SeriesRC')
        device_database.put_double('series_resistance',R)
        device_database.put_double('capacitance'      ,C)      
        device=EnergyStorageDevice(device_database)
        # setup experiment and measure 
        ragone_database=setup_expertiment()
        ragone_chart_data=measure_ragone_chart(device,ragone_database)
        # extract data
        E_computed=ragone_chart_data['energy']
        P         =ragone_chart_data['power' ]
        # compute the exact solution
        U_0=U_i/2+sqrt(U_i**2/4-R*P)
        E_exact=C/2*(-R*P*log(U_0**2/U_f**2)+U_0**2-U_f**2)
        # ensure the error is small
        max_percent_error=100*linalg.norm((E_computed-E_exact)/E_computed,inf)
        self.assertLess(max_percent_error,0.1)
    def testParallelRC(self):
        # make parallel RC equivalent circuit
        device_database=PropertyTree()
        device_database.put_string('type','ParallelRC')
        device_database.put_double('series_resistance'  ,R  )
        device_database.put_double('parallel_resistance',R_L)
        device_database.put_double('capacitance'        ,C  )      
        device=EnergyStorageDevice(device_database)
        # setup experiment and measure 
        ragone_database=setup_expertiment()
        ragone_chart_data=measure_ragone_chart(device,ragone_database)
        # extract data
        E_computed=ragone_chart_data['energy']
        P         =ragone_chart_data['power' ]
        # compute the exact solution
        U_0=U_i/2+sqrt(U_i**2/4-R*P)
        tmp=(U_f**2/R_L+P*(1+R/R_L))/(U_0**2/R_L+P*(1+R/R_L))
        E_exact=C/2*(-R_L*P*log(tmp)-R*R_L/(R+R_L)*P*log(tmp*U_0**2/U_f**2))
        # ensure the error is small
        max_percent_error=100*linalg.norm((E_computed-E_exact)/E_computed,inf)
        self.assertLess(max_percent_error,0.1)

if __name__ == '__main__':
    unittest.main()

