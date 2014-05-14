.. _GEO Data Sets:

GEO Data Sets
=============

.. image:: ../../orangecontrib/bio/widgets/icons/GEODataSets.svg
   :alt: GEO Data Sets widget icon
   :class: widget-category-bioinformatics widget-icon

Provides access to data sets from gene expression omnibus
(`GEO DataSets <http://www.ncbi.nlm.nih.gov/geo/>`_).
   
Signals
-------

Inputs:
   - (None)

Outputs:
   - :obj:`Data`
         Attribute-valued data set created in the widget.


Description
-----------

Data sets from GEO (officially,
`GEO DataSets <http://www.ncbi.nlm.nih.gov/geo/>`_) is a data base of gene
expression curated profiles maintained by `NCBI <http://www.ncbi.nlm.nih.gov/>`_
and included in
`Gene Expression Omnibus <http://www.ncbi.nlm.nih.gov/geo/>`_). This Orange
widget provides an access to all its data sets and outputs a data table
that can be further mined in Orange. For convenience, each data set selected
for further processing is stored locally. Upon its first access the widget
loads the data from GEO, but then, in any further queries, supports the
offline access.

.. image:: images/geodatasets-stamped.png
   :alt: GEO Data Sets widget

.. rst-class:: stamp-list

   1. Information on the number of data sets in GEO data base and the number of
      data sets that are available locally and were already loaded to a local
      computer.
   #. Type the text to be matched with description of the data sets. Only the
      matching data sets will be displayed.
   #. Available data sets matching the search string in the :obj:`Filter`.
      Locally available (previously downloaded) data sets are marked with a
      bullet. Click on an ID (e.g., GDS817) to open an original GEO
      web page with a data set description, or on a PubMed ID for an original
      publication. Select a line with a data set to display additional
      information (:obj:`Description`, :obj:`Sample Annotations`).
   #. Rows in the output table can include either gene profiles (columns
      represent samples) or sample profiles (columns represent microarray
      spots or genes).
   #. Click on the row of the data set to select it and display description
      of the data set and its annotations.
   #. Press :obj:`Commit` to push the data of the selected data set
      to the output of the widget. This
      action will also download the data from the GEO's web site if this
      is not already present on the local computer.
   #. Textual description of the data set.
   #. Sample annotations. Selection of annotations is important when saving
      samples as rows, as these define the class of the data instance to be
      mined by supervised methods (e.g., classification).

Examples
--------

