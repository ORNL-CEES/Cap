Equivalent circuits
===================

Series RC
---------

.. figure:: series_rc.png

A resistor and a capacitor are connected in series (denoted :math:`\mathrm{ESR}` 
and :math:`\mathrm{C}` in the figure above).

.. literalinclude:: series_rc.xml
   :language: xml
   
Above is the database to build a :math:`\mathrm{3\ F}` capacitor in series with a 
:math:`50\ \mathrm{m\Omega}` resistance.

Parallel RC
-----------

.. figure:: parallel_rc.png
   
An extra resistance is placed in parallel of the capacitor. It can be instantiated 
by the following database.

.. literalinclude:: parallel_rc.xml
   :language: xml

``type`` has been changed from ``SeriesRC`` to ``ParallelRC``. A :math:`2.5\ \mathrm{M\Omega}` 
leakage resistance is specified.

