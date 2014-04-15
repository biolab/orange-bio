.. _Databases:

Databases
=========

.. image:: ../../orangecontrib/bio/widgets/icons/Databases.svg
   :alt: Databases widget icon
   :class: widget-category-bioinformatics widget-icon

Updates local systems biology databases, like gene ontologies, annotations,
gene names, protein interaction networks, and similar.
   
Signals
-------

Inputs:
   - (None)

Outputs:
   - (None)


.. _my-reference-label:

Description
-----------

Many widgets in Orange bioinformatics add-on rely on
information on genes, gene sets, pathways, and alike. This information is
stored on your local computer when the widget requires them for the first
time. The corresponding data comes from different web resources, and is
either preprocessed and then stored on Orange server, or accessed directly
from a dedicated web site.

Orange does not change the data on your local
computer, and with time this becomes different to the newest version of the
online data sets. Databases widget can update the data on your local machine,
and can also be used to manage (remove or add) any locally stored
systems biology data set.

.. image:: images/databases-stamped.png
   :alt: Databases widget

.. rst-class:: stamp-list

   1. A list of available data bases matched by the search string in the
      filter.
   #. Checked if database is on your local computer. To load a new database
      that is not yet present on local computer, find the data base and
      click on the (unchecked) check box.
   #. This database needs requires an update, there is a newer version
      on the server.
   #. This database does not require un update, the version on the local
      computer is the same as the version on the server.
   #. Type any text to filter display only the matching databases. Matched
      are either the name of the data base or any of its tags.
   #. Update all the local databases.
   #. Download all the databases currently displayed in the table. If a
      search string in filter is defined, only the matching databases will
      be loaded.
   #. Info on a listed databases, reporting on local data bases and
      those on the server.
   #. Some databases may require access key.

To get a more detailed information on the particular database that requires
an update, hover on its :obj:`Update` button.

.. image:: images/databases-hover.png
   :alt: Databases widget
