from obiBioMart import *

## Printing attribute configurations 

connection = BioMartConnection("http://www.biomart.org/biomart/martservice")
registry = connection.registry()
#for schema in registry.virtual_schemas()[:1]:
#    for database in schema.marts()[:1]:
#        for dataset in database.datasets()[:2]:
#            for attrTree in dataset.configuration().attributes():
#                if not getattr(attrTree, "hidden", "false") == "true":
#                    print dataset.name, "has attribute", getattr(attrTree, "displayName", "<unknown>")
                    
## Printing dataset attributes

database = registry["ensembl"]
dataset = database["hsapiens_gene_ensembl"]

for attr in dataset.attributes():
    print attr

for filter in dataset.filters():
    print filter
                    
query = BioMartQuery(connection, dataset="hsapiens_gene_ensembl", attributes=["ensembl_transcript_id", "chromosome_name"], 
                     filters=[("chromosome_name", ["22"])])
print query.get_count()

print query.run()

query = BioMartQuery(connection)
query.set_dataset("hsapiens_gene_ensembl")
query.add_filter("chromosome_name", ["22"])
query.add_attribute("ensembl_transcript_id")
query.add_attribute("chromosome_name")
query.add_attribute("uniprot_swissprot")
print query.get_count()
print query.run()