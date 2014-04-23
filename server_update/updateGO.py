##!interval=7
##!contact=ales.erjavec@fri.uni-lj.si

import urllib2
import re
import cPickle
import tarfile

from datetime import datetime
from collections import defaultdict

from common import *

from Orange.bio import obiGO, obiTaxonomy, obiGene


tmp_path = os.path.join(environ.buffer_dir, "tmp_GO")
try:
    os.makedirs(tmp_path)
except Exception:
    pass

uncompressedSize = lambda filename: sum(info.size for info in tarfile.open(filename).getmembers())

DATE_FMT_1 = "%Y-%m-%d %H:%M:%S.%f"
DATE_FMT_2 = "%Y-%m-%d %H:%M:%S"


def info_date_time_parse(time):
    """
    Parse a "datetime" field from the sf info record into a datetime.datetime.
    """
    try:
        return datetime.strptime(time, DATE_FMT_1)
    except ValueError:
        return datetime.strptime(time, DATE_FMT_2)


def http_last_modified(url):
    """
    Retrieve a "last-modified" time for the url as a datetime.datetime object.
    """
    stream = urllib2.urlopen(url)
    return datetime.strptime(stream.headers.get("Last-Modified"),
                             "%a, %d %b %Y %H:%M:%S %Z")


def list_available_organisms():
    """
    Return a list of all available GO organism codes.
    """
    source = urllib2.urlopen("http://www.geneontology.org/gene-associations/").read()
    codes = re.findall("gene_association\.([a-zA-z0-9_]+?)\.gz", source)
    return sorted(set(codes))


def sf_org_mtime(org_code):
    info = sf_server.info("GO", "gene_association.{}.tar.gz".format(org_code))
    return info_date_time_parse(info["datetime"])


def web_org_mtime(org_code):
    return http_last_modified(
        "http://www.geneontology.org/gene-associations/gene_association.{}.gz"
        .format(org_code))


def sf_ontology_mtime():
    info = sf_server.info("GO", "gene_ontology_edit.obo.tar.gz")
    return info_date_time_parse(info["datetime"])


def web_ontology_mtime():
    return http_last_modified(
        "http://www.geneontology.org/ontology/gene_ontology.obo")


if web_ontology_mtime() > sf_ontology_mtime():
    print "donwloading ontology"
    filename = os.path.join(tmp_path, "gene_ontology_edit.obo.tar.gz")
    obiGO.Ontology.DownloadOntology(filename)

    ##load the ontology to test it
    o = obiGO.Ontology(filename)
    del o
    ##upload the ontology
    print "Uploading gene_ontology_edit.obo.tar.gz"
    sf_server.upload(
        "GO", "gene_ontology_edit.obo.tar.gz", filename,
        title="Gene Ontology (GO)",
        tags=["gene", "ontology", "GO", "essential",
              "#uncompressed:%i" % uncompressedSize(filename),
              "#version:%i" % obiGO.Ontology.version]
    )
    sf_server.unprotect("GO", "gene_ontology_edit.obo.tar.gz")


orgMap = {"352472": "44689", "562": "83333", "3055": None,
          "7955": None, "11103": None, "2104": None, "4754":
          None, "31033": None, "8355": None, "4577": None}

commonOrgs = dict([(obiGO.from_taxid(id), id)
                   for id in obiTaxonomy.common_taxids()
                   if obiGO.from_taxid(id) != None])

essentialOrgs = [obiGO.from_taxid(id)
                 for id in obiTaxonomy.essential_taxids()]

exclude = ["goa_uniprot", "goa_pdb", "GeneDB_tsetse", "reactome",
           "goa_zebrafish", "goa_rat", "goa_mouse"]

updatedTaxonomy = defaultdict(set)


for org in list_available_organisms():

    if org in exclude or org not in commonOrgs:
        continue

    if web_org_mtime(org) <= sf_org_mtime(org):
        # Skip update
        continue

    print "Updating", org

    filename = os.path.join(tmp_path, "gene_association." + org + ".tar.gz")
    obiGO.Annotations.DownloadAnnotations(org, filename)

    ## Load the annotations to test them and collect all taxon ids from them
    print filename
    a = obiGO.Annotations(filename, genematcher=obiGene.GMDirect())
    taxons = set([ann.Taxon for ann in a.annotations])
    ## exclude taxons with cardinality 2
    taxons = [tax for tax in taxons if "|" not in tax]
    for tax in taxons:
        taxid = tax.split(":", 1)[-1]
        updatedTaxonomy[taxid].add(org)
    del a

    orgName = obiTaxonomy.name(commonOrgs[org])
    taxid = obiTaxonomy.taxname_to_taxid(orgName)

    print "Uploading", "gene_association." + org + ".tar.gz"
    sf_server.upload(
        "GO", "gene_association." + org + ".tar.gz", filename,
        title="GO Annotations for " + orgName,
        tags=["gene", "annotation", "ontology", "GO", orgName,
              "#uncompressed:%i" % uncompressedSize(filename),
              "#organism:" + orgName,
              "#version:%i" % obiGO.Annotations.version] +
             (["essential"] if org in essentialOrgs else []) +
             obiTaxonomy.shortname(taxid)
    )
    sf_server.unprotect("GO", "gene_association." + org + ".tar.gz")

try:
    tax = cPickle.load(open(sf_local.localpath_download("GO", "taxonomy.pickle"), "rb"))
except Exception:
    tax = {}

## Upload taxonomy if any differences in the updated taxonomy
if any(tax.get(key, set()) != updatedTaxonomy.get(key, set())
       for key in set(updatedTaxonomy)):
    tax.update(updatedTaxonomy)
    cPickle.dump(tax, open(os.path.join(tmp_path, "taxonomy.pickle"), "wb"))
    print "Uploading", "taxonomy.pickle"
    sf_server.upload(
        "GO", "taxonomy.pickle", os.path.join(tmp_path, "taxonomy.pickle"),
        title="GO taxon IDs",
        tags=["GO", "taxon", "organism", "essential",
              "#version:%i" % obiGO.Taxonomy.version])
    sf_server.unprotect("GO", "taxonomy.pickle")
