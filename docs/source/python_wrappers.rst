Python wrappers
===============

.. testcode::

    >>> import pycap
    >>> device_database=pycap.PropertyTree()
    >>> device_database.parse_xml('device.xml')
    >>> device=pycap.EnergyStorageDevice(device_database)

