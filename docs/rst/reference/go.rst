.. py:currentmodule:: Orange.bio.go

=========================
Gene Ontology (:mod:`go`)
=========================

Provides access to `Gene Ontology`_ and it's gene annotations.

.. _Gene Ontology: http://geneontology.org/


.. autoclass:: Orange.bio.go.Ontology(filename=None, progress_callback=None, rev=None)
   :members:
   :member-order: bysource
   :exclude-members:
      Load, ExtractSubGraph, ExtractSuperGraph, GetDefinedSlimsSubsets,
      GetSlimTerms, GetSlimsSubset, GetTermDepth, ParseFile, SetSlimSubsets,
      aliasMapper, reverseAliasMapper, slimsSubset, parse_file, load


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


.. autoclass:: Orange.bio.go.Term
   :members:


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
