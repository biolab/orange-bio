PIPA
====

.. image:: ../../../widgets/icons/PIPA.png
   :alt: widget icon
   
Signals
-------

Inputs:
   - None

Outputs:
   - Selected microarrays (ExampleTable)
        Selected microarrays from PIPA. Each annotated column contains results of a single microarray experiment.

Description
-----------

.. image:: PIPA.*
   :alt: PIPA widget

Examples
--------

Any of your schemas should probably start with the File_ widget. In the schema below, 
the widget is used to read the data that is then sent to both `Data Table`_ widget and 
to widget that displays `Attribute statistics`_.

.. image:: PIPA_schema.*
   :alt: Example schema with PIPA widget
   
.. image:: PIPA_datatable.*
   :alt: Output of PIPA in a data table.
 
