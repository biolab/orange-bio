.. py:currentmodule:: Orange.bio.go
.. py:module:: Orange.bio.go

=========================
Gene Ontology (:mod:`go`)
=========================


Provides access to `Gene Ontology`_ and its gene annotations.

.. _Gene Ontology: http://geneontology.org/


.. autoclass:: Orange.bio.go.Ontology(filename=None, progress_callback=None, rev=None)
   :members:
      defined_slims_subsets, named_slims_subset, set_slims_subset, slims_for_term,
      extract_super_graph, extract_sub_graph, __getitem__, __len__,
      __iter__, __contains__

   Ontology supports a subset of the Mapping protocol:

   >>> term_ids = list(ontology)
   >>> term = ontology[term_ids[0]]


.. autoclass:: Orange.bio.go.Term

   .. attribute:: id

   	  The term id.

   .. attribute:: namespace

      The namespace of the term.

   .. attribute:: def_

      The term definition (Note the use of trailing underscore
      to avoid conflict with a python keyword).

   .. attribute:: is_a

      List of term ids this term is a subterm of (parent terms).

   .. attribute:: related

      List of (rel_type, term_id) tuples with rel_type specifying
      the relationship type with term_id.


.. autoclass:: Orange.bio.go.Annotations(filename_or_organism=None, ontology=None, genematcher=None, progress_callback=None, rev=None)
   :members:
   :member-order: bysource
   :exclude-members:
      set_ontology, load, Load, parse_file, ParseFile, RemapGenes,
      AddAnnotation, DrawEnrichmentGraph, GetAllAnnotations, GetAllGenes,
      GetAnnotatedTerms, GetEnrichedTerms, GetGeneNamesTranslator,
      GetOntology, SetOntology, aliasMapper, allAnnotations,
      geneAnnotations, geneNames, geneNamesDict, termAnnotations


.. autoclass:: Orange.bio.go.AnnotationRecord
   :members:
       from_string,
   	   DB, DB_Object_ID, DB_Object_Symbol, Qualifier, GO_ID,
       DB_Reference, Evidence_Code, With_From, Aspect, DB_Object_Name,
       DB_Object_Synonym, DB_Object_Type, Taxon, Date, Assigned_By,
       Annotation_Extension, Gene_Product_Form_ID


Example
-------

Load the ontology and print out some terms::

	from Orange.bio import go
	ontology = go.Ontology()
	term = ontology["GO:0097194"] # execution phase of apoptosis

	# print a term
	print term

	# access fields by name
	print term.id, term.name
	# note the use of underscore due to a conflict with a python def keyword
	print term.def_


Searching the annotation (part of :download:`code/go_gene_annotations.py`)

.. literalinclude:: code/go_gene_annotations.py

Term enrichment (part of :download:`code/go_enrichment.py`)

.. literalinclude:: code/go_enrichment.py
   :lines: 6-
