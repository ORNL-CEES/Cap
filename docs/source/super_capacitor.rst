Supercapacitors
---------------

``type`` is set to ``SuperCapacitor``.
``dim`` is used to select two- or three-dimensional simulations.

.. code::

    device {
        type SuperCapacitor
        dim 2
        geometry {
            [...]
        }
        material_properties {
            [...]
        }
    }


Geometry
^^^^^^^^

.. code::

    geometry {
        type supercapacitor

        anode_collector_thickness    5.0e-4 ; [centimeter]
        anode_electrode_thickness   50.0e-4 ; [centimeter]
        separator_thickness         25.0e-4 ; [centimeter]
        cathode_electrode_thickness 50.0e-4 ; [centimeter]
        cathode_collector_thickness  5.0e-4 ; [centimeter]
        geometric_area              25.0e-2 ; [square centimeter]
    }

The thickness of each layer in the sandwich (anode collector, anode electrode,
separator, cathode electrode, cathode current collector) can be adjusted
independently from one another. The specified cross-sectional area applies to
the whole stack.

.. figure:: sandwich.png
    :figwidth: 400px
    :align: center

    Schematic representation of the supercapacitor conventional sandwich-like
    configuration.
    1: anode electrode,
    2: separator,
    3: cathode electrode,
    4: anode collector,
    5: cathode collector.


Governing equations
^^^^^^^^^^^^^^^^^^^

.. |solid_current| replace::
    :math:`i_1 = -\sigma \nabla \Phi_1`

.. |liquid_current| replace::
    :math:`i_2 = -\kappa \nabla \Phi_2`

.. |interfacial_current| replace::
    :math:`-\nabla \cdot i_1 = \nabla \cdot i_2 = a i_n`

+------------------------------+-----------------------+------------------------------+
| collector                    | electrode             | separator                    |
+==============================+=======================+==============================+
|                              |                       |                              |
| |solid_current|              | |solid_current|       | |liquid_current|             |
|                              |                       |                              |
| :math:`\nabla \cdot i_1 = 0` | |liquid_current|      | :math:`\nabla \cdot i_2 = 0` |
|                              |                       |                              |
|                              | |interfacial_current| |                              |
|                              |                       |                              |
+------------------------------+-----------------------+------------------------------+


.. |alias1| replace::
    :math:`-\sigma \left.\frac{\partial \Phi_1}{\partial n}\right|_c
    = -\sigma \left.\frac{\partial \Phi_1}{\partial n}\right|_e`

.. |alias2| replace::
    :math:`-\kappa \left.\frac{\partial \Phi_2}{\partial n}\right|_e
    = -\kappa \left.\frac{\partial \Phi_2}{\partial n}\right|_s`

.. |alias3| replace::
    :math:`0 = -\kappa \left. \frac{\partial \Phi_2}{\partial n} \right|_e`

.. |alias4| replace::
    :math:`-\sigma \left. \frac{\partial \Phi_1}{\partial n} \right|_e = 0`

+-------------------------------+--------------------------------+
| collector-electrode interface |  electrode-separator interface |
+===============================+================================+
|                               |                                |
| |alias3|                      | |alias2|                       |
|                               |                                |
| |alias1|                      | |alias4|                       |
|                               |                                |
+-------------------------------+--------------------------------+

+----------------------------------------------------------------+
| boundary collector tab                                         |
+================================================================+
|:math:`\Phi_1 = U`                                              |
|                                                                |
|or                                                              |
|                                                                |
|:math:`-\sigma \frac{\partial \Phi_1}{\partial n} = I/S`        |
|                                                                |
|or                                                              |
|                                                                |
|:math:`-\sigma \frac{\partial \Phi_1}{\partial n} \Phi_1 = P/S` |
|                                                                |
|or                                                              |
|                                                                |
|:math:`-\sigma \frac{\partial \Phi_1}{\partial n} R S = \Phi_1` |
|                                                                |
+----------------------------------------------------------------+

Ignoring the influence of the electrolyte concentration, the current density
in the matrix and solution phases can be expressed by Ohm’s law as

.. math::

    i_1 = -\sigma \nabla \Phi_1

    i_2 = -\kappa \nabla \Phi_2

:math:`i` and :math:`\Phi` represent current density and potential; subscript
indices :math:`1` and :math:`2` denote respectively the solid and the liquid
phases. :math:`\sigma` and :math:`\kappa` are the matrix and solution phase
conductivities.

The total current density is given by :math:`i = i_1 + i_2`. Conservation of
charge dictates that

.. math::

    -\nabla \cdot i_1 = \nabla \cdot i_2 = a i_n

where :math:`a` is the interfacial area per unit volume and the current
transferred from the matrix phase to the electrolyte :math:`i_n` is the sum of
the double-layer the faradaic currents

.. math::

    i_n = C \frac{\partial}{\partial t} \left(\Phi_1 - \Phi_2\right)
        + i_0 \left( e^{\frac{\alpha_a F}{RT}\eta}
            - e^{-\frac{\alpha_c F}{RT}\eta} \right)

:math:`C` is the double-layer capacitance. :math:`i_0` is the exchange current
density, :math:`\alpha_a` and :math:`\alpha_c` the anodic and cathodic charge
transfer coefficients, respectively. :math:`F`, :math:`R`, and :math:`T` stand
for Faraday’s constant, the universal gas constant and temperature.
:math:`\eta` is the overpotential relative to the equilibrium potential
:math:`U_{eq}`

.. math::

    \eta = \Phi_1 - \Phi_2 - U_{eq}


Material properties
^^^^^^^^^^^^^^^^^^^

.. code::

    material_properties {
        anode {
            type           porous_electrode
            matrix_phase   electrode_material
            solution_phase electrolyte
        }
        cathode {
            type           porous_electrode
            matrix_phase   electrode_material
            solution_phase electrolyte
        }
        separator {
            type           permeable_membrane
            matrix_phase   separator_material
            solution_phase electrolyte
        }
        collector {
            type           current_collector
            metal_foil     collector_material
        }

        separator_material {
            void_volume_fraction             0.6       ;
            tortuosity_factor                1.29      ;
            pores_characteristic_dimension   1.5e-7    ; [centimeter]
            pores_geometry_factor            2.0       ;
            mass_density                     3.2       ; [gram per cubic centimeter]
            heat_capacity                    1.2528e3  ; [joule per kilogram kelvin]
            thermal_conductivity             0.0019e2  ; [watt per meter kelvin]
        }
        electrode_material {
            differential_capacitance         3.134     ; [microfarad per square centimeter]
            exchange_current_density         7.463e-10 ; [ampere per square centimeter]
            void_volume_fraction             0.67      ;
            tortuosity_factor                2.3       ;
            pores_characteristic_dimension   1.5e-7    ; [centimeter]
            pores_geometry_factor            2.0       ;
            mass_density                     2.3       ; [gram per cubic centimeter]
            electrical_resistivity           1.92      ; [ohm centimeter]
            heat_capacity                    0.93e3    ; [joule per kilogram kelvin]
            thermal_conductivity             0.0011e2  ; [watt per meter kelvin]
        }
        collector_material {
            mass_density                     2.7       ; [gram per cubic centimeter]
            electrical_resistivity          28.2e-7    ; [ohm centimeter]
            heat_capacity                    2.7e3     ; [joule per kilogram kelvin]
            thermal_conductivity           237.0       ; [watt per meter kelvin]
        }
        electrolyte {
            mass_density                     1.2       ; [gram per cubic centimeter]
            electrical_resistivity           1.49e3    ; [ohm centimeter]
            heat_capacity                    0.0       ; [joule per kilogram kelvin]
            thermal_conductivity             0.0       ; [watt per meter kelvin]
        }
    }

The specific surface area per unit volume :math:`a` is estimated using

.. math::

    a = \frac{(1+\zeta)\varepsilon}{r}

where :math:`\zeta` is the pore's geometry factor (:math:`\zeta=2` for
spheres, :math:`1` for cylinders, and :math:`0` for slabs) and :math:`r` is
the pore's characteristic dimension.
[M. W. Verbrugge and B. J. Koch, J. Electrochem. Soc., 150, A374 2003]

The solution electrical conductivity :math:`\kappa` incorporates the effect
of porosity and tortuosity

.. math::

    \kappa = \frac{\kappa_\infty \varepsilon}{\Gamma}

where :math:`\kappa_\infty` is the liquid phase (free solution) conductivity,
:math:`\varepsilon` is the void volume fraction, and :math:`\kappa` is the
tortuosity factor.

The solid phase conductivity is also corrected for porosity (and tortuosity???)

.. math::

    \sigma = \sigma_\infty (1-\varepsilon)
