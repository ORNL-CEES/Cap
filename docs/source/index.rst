.. cap documentation master file, created by
   sphinx-quickstart on Fri Sep 11 15:54:53 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to cap's documentation!
===============================

Contents:

.. toctree::
   :maxdepth: 2

   faq

.. math::
   \sum_k \frac{2\pi}{k}

Series RC
---------

.. figure:: series_rc.png

A resistor and a capacitor are connected in series (denoted $ESR$ and $C$ in the figure above).
::

   <device>
       <type>              SeriesRC </type>
       <series_resistance>  50.0e-3 </series_resistance> <!-- ohm -->
       <capacitance>         3.0    </capacitance>       <!-- farad -->
   </device>
   
Above is the database to build a 3 F capacitor in series with a 50 mΩ resistance.

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

.. testcode::

    >>> import pycap
    >>> device_database=pycap.PropertyTree()
    >>> device_database.parse_xml('device.xml')
    >>> device=pycap.EnergyStorageDevice(device_database)

Grid table:

+------------+------------+-----------+ 
| Header 1   | Header 2   | Header 3  | 
+============+============+===========+ 
| body row 1 | column 2   | column 3  | 
+------------+------------+-----------+ 
| body row 2 | Cells may span columns.| 
+------------+------------+-----------+ 
| body row 3 | Cells may  | - Cells   | 
+------------+ span rows. | - contain | 
| body row 4 |            | - blocks. | 
+------------+------------+-----------+

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

