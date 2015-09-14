Supercapacitors
===============

``type`` is set to ``SuperCapacitor``.
``dim`` is used to select two- or three-dimensional simulations.

.. code:: xml

    <device>
        <type>SuperCapacitor</type>
        <dim>2</dim>
        <geometry>
            [...]
        </geometry>
        <material_properties>
            [...]
        </material_properties>
    </device>


Geometry
--------

.. literalinclude:: super_capacitor.xml
    :language: xml
    :start-after: <!-- geometry begin -->
    :end-before: <!-- geometry end -->

``mesh_file`` give the path to the triangulation. The dimension has to match
``dim`` or an exception will be thrown.
The width of each layer in the sandwich (anode collector, anode electrode,
separator, cathode electrode, cathode current collector) can be adjusted
independently from one another. The overall sandwich height and depth (in
3-D) can be changed as well.

.. figure:: sandwich.png
    :figwidth: 400px
    :align: center

    This is the caption.

Material properties
-------------------

.. literalinclude:: super_capacitor.xml
    :language: xml
    :start-after: <!-- material_properties begin -->
    :end-before: <!-- material_properties end -->

