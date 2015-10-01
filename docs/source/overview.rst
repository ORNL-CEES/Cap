Getting started
===============

Overview
--------

`cap <https://github.com/dalg24/cap>`_ is a library for modelling energy
storage devices.
Its core is implemented in C++ but Python wrappers are also available.

.. testcode::

    >>> from pycap import PropertyTree,EnergyStorageDevice
    >>> input_database=PropertyTree()
    >>> input_database.parse_xml('super_capacitor.xml')
    >>> device=EnergyStorageDevice(input_database.get_child('device'))
