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


Alternative to the full install procedure
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
All you need is a working installation of Docker.

The following command starts a Docker container with a Jupyter Notebook server listening for HTTP connections on port 8888.
It has pycap already installed on it and comes with a few notebooks as example.
::

    docker run -d -p 8888:8888 dalg24/cap

Open your web browser and follow http://192.168.99.100:8888/tree/cap-notebooks
