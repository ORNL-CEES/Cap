Equivalent circuits
===================

Series RC
---------

.. figure:: series_rc.png

A resistor and a capacitor are connected in series (denoted :math:`\mathrm{ESR}` 
and :math:`\mathrm{C}` in the figure above).
::
   <device>
       <type>              SeriesRC </type>
       <series_resistance>  50.0e-3 </series_resistance> <!-- ohm -->
       <capacitance>         3.0    </capacitance>       <!-- farad -->
   </device>
   
Above is the database to build a 3 F capacitor in series with a :math:`50 \mili\ohm` resistance.

Parallel RC
-----------

.. figure:: parallel_rc.png
   
An extra resistance is placed in parallel of the capacitor. It can be instantiated by the following database.
::
   <device>
       <type>                ParallelRC </type>
       <parallel_resistance>     2.5e+6 </parallel_resistance> <!-- ohm -->
       <series_resistance>      50.0e-3 </series_resistance>   <!-- ohm -->
       <capacitance>             3.0    </capacitance>         <!-- farad -->
   </device>
type has been changed from SeriesRC to ParallelRC. A 2.5 MΩ leakage resistance is specified.

Using pycap
-----------

.. testcode::

    >>> import pycap
    >>> device_database=pycap.PropertyTree()
    >>> device_database.parse_xml('device.xml')
    >>> device=pycap.EnergyStorageDevice(device_database)
