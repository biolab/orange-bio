##interval:7
from common import *
from Orange.bio import obiGene, obiTaxonomy
from gzip import GzipFile

tmpdir = os.path.join(environ.buffer_dir, "tmp_NCBIGene_info")
try:
    os.mkdir(tmpdir)
except Exception, ex:
    pass

gene_info_filename = os.path.join(tmpdir, "gene_info")
gene_history_filename = os.path.join(tmpdir, "gene_history")

obiGene.NCBIGeneInfo.get_geneinfo_from_ncbi(gene_info_filename)
obiGene.NCBIGeneInfo.get_gene_history_from_ncbi(gene_history_filename)

info = open(gene_info_filename, "rb")
hist = open(gene_history_filename, "rb")

taxids = obiGene.NCBIGeneInfo.common_taxids()
essential = obiGene.NCBIGeneInfo.essential_taxids()

genes = dict([(taxid, []) for taxid in taxids])
for gi in info:
    if any(gi.startswith(id + "\t") for id in taxids):
        genes[gi.split("\t", 1)[0]].append(gi.strip())

history = dict([(taxid, []) for taxid in taxids])
for hi in hist:
    if any(hi.startswith(id + "\t") for id in taxids): 
        history[hi.split("\t", 1)[0]].append(hi.strip())

for taxid, genes in genes.items():
    filename = os.path.join(tmpdir, "gene_info.%s.db" % taxid)
    f = open(filename, "wb")
    f.write("\n".join(genes))
    f.flush()
    f.close()
    print "Uploading", filename
    sf_server.upload("NCBI_geneinfo", "gene_info.%s.db" % taxid, filename,
              title = "NCBI gene info for %s" % obiTaxonomy.name(taxid),
              tags = ["NCBI", "gene info", "gene_names", obiTaxonomy.name(taxid)] + (["essential"] if taxid in essential else []))
    sf_server.unprotect("NCBI_geneinfo", "gene_info.%s.db" % taxid)
    
    filename = os.path.join(tmpdir, "gene_history.%s.db" % taxid)
    f = open(filename, "wb")
    f.write("\n".join(history.get(taxid, "")))
    f.flush()
    f.close()
    print "Uploading", filename
    sf_server.upload("NCBI_geneinfo", "gene_history.%s.db" % taxid, filename,
              title = "NCBI gene history for %s" % obiTaxonomy.name(taxid),
              tags = ["NCBI", "gene info", "history", "gene_names", obiTaxonomy.name(taxid)] + (["essential"] if taxid in essential else []))
    sf_server.unprotect("NCBI_geneinfo", "gene_history.%s.db" % taxid)
