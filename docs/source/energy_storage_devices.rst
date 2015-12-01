Energy storage devices
======================

The class :py:class:`pycap.EnegyStorageDevice` is an abstract representation
for an energy storage device. It can evolve in time at various operating
conditions and return the voltage drop across itself and the electrical
current that flows through it.

.. py:class:: pycap.EnergyStorageDevice(ptree)

    :param pycap.PropertyTree ptree: Tree to construct the appropriate energy
        storage device

    Wrappers for the abstract representation of an energy storage device.

    .. py:method:: evolve_one_time_step_constant_voltage(...)

        evolve_one_time_step_constant_voltage( (float)time_step, (float)voltage) -> None :
            Impose the voltage and evolve one time step.

    .. py:method:: evolve_one_time_step_constant_current(...)

        evolve_one_time_step_constant_current( (float)time_step, (float)current) -> None :
            Impose the current and evolve one time step.

    .. py:method:: evolve_one_time_step_constant_load(...)

        evolve_one_time_step_constant_load( (float)time_step, (float)load) -> None :
            Impose the load and evolve one time step.

    .. py:method:: evolve_one_time_step_constant_power(...)

        evolve_one_time_step_constant_power( (float)time_step, (float)power) -> None :
            Impose the power and evolve one time step.

    .. py:method:: get_voltage(...)

        get_voltage( ) -> float :
             Return the voltage drop across the device in units of volts. 

    .. py:method:: get_current(...)

        get_current( ) -> float :
             Return the electrical current that flows through the device in
             units of amperes.

    .. py:staticmethod:: compute_equivalent_circuit(...)

        compute_equivalent_circuit( (pycap.PropertyTree)ptree) -> pycap.PropertyTree :
             Return the tree to build an equivalent circuit.

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
