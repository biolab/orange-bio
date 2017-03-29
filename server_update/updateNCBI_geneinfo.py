from server_update import *
from server_update.tests.test_NCBIGeneInfo import NCBIGeneInfoTest
from orangecontrib.bio import obiGene, obiTaxonomy


DOMAIN = 'NCBI_geneinfo'
domain_path = sf_local.localpath(DOMAIN)
create_folder(domain_path)

gene_info_filename = os.path.join(domain_path, "gene_info")
gene_history_filename = os.path.join(domain_path, "gene_history")


taxids = obiGene.NCBIGeneInfo.common_taxids()
essential = obiGene.NCBIGeneInfo.essential_taxids()
genes = dict([(taxid, []) for taxid in taxids])
history = dict([(taxid, []) for taxid in taxids])

print("Downloading geneinfo file...")
obiGene.NCBIGeneInfo.get_geneinfo_from_ncbi(gene_info_filename)
with open(gene_info_filename, 'r') as gene_info_temp:
    print("creating genes info dict...")
    for gi in gene_info_temp:
        if any(gi.startswith(id + "\t") for id in taxids):
            genes[gi.split("\t", 1)[0]].append(gi.strip())


print("Downloading gene history file...")
obiGene.NCBIGeneInfo.get_gene_history_from_ncbi(gene_history_filename)
with open(gene_history_filename, 'r') as gene_history_temp:
    print("creating genes history dict...")
    for hi in gene_history_temp:
        if any(hi.startswith(id + "\t") for id in taxids):
            history[hi.split("\t", 1)[0]].append(hi.strip())

print("Done!")


for taxid, genes in genes.items():
    TAGS_ORGANISM = [obiTaxonomy.name(taxid)] + obiTaxonomy.shortname(taxid) + \
                    (["essential"] if taxid in essential else [])

    TAGS_INFO = ["NCBI", "gene info", "gene_names"] + TAGS_ORGANISM
    TAGS_HIST = ["NCBI", "gene info", "history", "gene_names"] + TAGS_ORGANISM

    filename = os.path.join(domain_path, "gene_info.%s.db" % taxid)
    with open(filename, 'w') as f:
        f.write("\n".join(genes))
        f.flush()
        f.close()
        print("{} created".format(filename))
        create_info_file(filename, title="NCBI gene info for %s" % obiTaxonomy.name(taxid), tags=TAGS_INFO)
        print("{}.info file created".format(filename))

    filename = os.path.join(domain_path, "gene_history.%s.db" % taxid)
    with open(filename, 'w') as f:
        f.write("\n".join(history.get(taxid, '')))
        f.flush()
        f.close()
        print("{} created".format(filename))
        create_info_file(filename, title="NCBI gene history for %s" % obiTaxonomy.name(taxid), tags=TAGS_HIST)
        print("{}.info file created".format(filename))

helper = SyncHelper(DOMAIN, NCBIGeneInfoTest)
helper.run_tests()
helper.sync_files()

helper.remove_update_folder()
