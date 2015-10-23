
.. code:: python

    from pycap import PropertyTree,EnergyStorageDevice,CyclicVoltammetry
    
    # setup the experiment
    ptree=PropertyTree()
    ptree.put_int   ('cycles',2)
    ptree.put_double('initial_voltage',0.0)
    ptree.put_double('final_voltage',0.0)
    ptree.put_double('scan_limit_1',2.4)
    ptree.put_double('scan_limit_2',-0.5)
    ptree.put_double('scan_rate',100e-3)
    ptree.put_double('step_size',5e-3)
    
    cv=CyclicVoltammetry(ptree)
    
    # build an energy storage device
    ptree=PropertyTree()
    ptree.put_string('type','SeriesRC')
    ptree.put_double('capacitance',3)
    ptree.put_double('series_resistance',50e-3)
    device=EnergyStorageDevice(ptree)
    
    from pycap import initialize_data,report_data,plot_data
    from pycap import plot_cyclic_voltammogram
    
    # run the experiment and visualize the measured data
    data=initialize_data()
    steps=cv.run(device,data)
    
    print "%d steps"%steps
    
    %matplotlib inline
    plot_data(data)
    plot_cyclic_voltammogram(data)



.. parsed-literal::

    2324 steps



.. image:: cv_files/cv_2_1.png


.. image:: cv_files/cv_2_3.png

