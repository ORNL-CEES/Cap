Energy storage devices
======================

.. figure:: energy_storage_device.png

Only the electrical current :math:`I` and voltage :math:`U` of the device are
measurable. Several operating conditions are possibles. One may want to
impose:

    - The voltage :math:`U` across the device.
    - The electrical current :math:`I` that flows through it.
    - The load :math:`R=U/I` the device is subject to.
    - The power :math:`P=UI`.

:class:`EnergyStorageDevice`
----------------------------

The class :py:class:`pycap.EnegyStorageDevice` is an abstract representation
for an energy storage device. It can evolve in time at various operating
conditions and return the voltage drop across itself and the electrical
current that flows through it.

.. py:class:: pycap.EnergyStorageDevice(ptree)

    Wrappers for the abstract representation of an energy storage device.

    :param pycap.PropertyTree ptree: Tree to construct the appropriate energy
        storage device.

    .. py:method:: evolve_one_time_step_constant_voltage(time_step,voltage)

        Impose the voltage and evolve one time step.

        :param float time_step: The time step in units of seconds.
        :param float voltage: The voltage in units of volts.

    .. py:method:: evolve_one_time_step_constant_current(time_step,current)

        Impose the current and evolve one time step.

        :param float time_step: The time step in units of seconds.
        :param float current: The electric current in units of amperes.

    .. py:method:: evolve_one_time_step_constant_load(time_step,load)

        Impose the load and evolve one time step.

        :param float time_step: The time step in units of seconds.
        :param float load: The load in units of ohms.

    .. py:method:: evolve_one_time_step_constant_power(time_step,power)

        Impose the power and evolve one time step.

        :param float time_step: The time step in units of seconds.
        :param float power: The power in units of watts.

    .. py:method:: get_voltage()

        Measure the voltage drop across the device. 

        :return: The voltage across the device in units of volts.
        :rtype: float

    .. py:method:: get_current()

        Measure the electrical current that flows through the device.

        :return: The electrical current in units of amperes.
        :rtype: float

    .. py:staticmethod:: compute_equivalent_circuit(ptree)

        Compute the equivalent circuit to a supercapacitor.

        :param pycap.PropertyTree ptree: The tree to build a supercapacitor.
        :return: The tree to build the equivalent circuit.
        :rtype: pycap.PropertyTree

.. note::

    The class methods :py:meth:`evolve_one_time_step_constant_xxx` have their
    conterpart :py:meth:`evolve_one_time_step_changing_xxx` where the
    operating condition is imposed on the device as a ramp function rather
    than a step function.

The rest of this section describes the energy storage devices that are
available in `cap`, namely:

    - Equivalent circuits

    - Supercapacitors

.. include:: equivalent_circuits.rst

.. include:: super_capacitor.rst

.. include:: battery.rst
