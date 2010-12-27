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
        Selected experiments. Each annotated column contains results of a single experiment or, if the corresponding option is chosen, the average of multiple replicates.

Description
-----------

.. image:: PIPA.*
   :alt: PIPA widget

The PIPA widget lists accessible experiments, which can be filtered with the "Search" box at the top. The selected experiments will appear on the output when the "Commit" button is clicked. You can connect the output of the PIPA widget to any Orange widget which accepts ExampleTable as input. The widget will automatically save (cache) downloaded data and you will therefore be able to analyse them offline.

To select multiple experiments click them while holding the "Control" key. For frequent combinations of selections use the "Experiment Sets" feature: select experiments and click on the "+" button.

The logarithmic transformation is computed as the binary logarithm of the (value + 1). If username and passwords are not given, only the public experiments will be accessible.

Examples
--------

In the schema below we connected PIPA to Data Table and GO Enrichment Analysis widgets.

.. image:: PIPA_schema.*
   :alt: Example schema with PIPA widget

The Data Table widget below contains the output from the PIPA widget. Each column contains gene expressions of a single experiment. The labels are shown in the table header.
  
.. image:: PIPA_datatable.*
   :alt: Output of PIPA in a data table.
 
