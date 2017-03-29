""" GO update """
import gzip
import io
import pickle
import re
import tarfile


from server_update import *
from collections import defaultdict
from urllib.request import urlopen
from server_update.tests.test_GO import GOTest
from orangecontrib.bio import obiGO, obiTaxonomy, obiGene


DOMAIN = 'GO'
download_path = sf_local.localpath('downloaded_files')
domain_path = sf_local.localpath(DOMAIN)
create_folder(domain_path)
create_folder(download_path)

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
    stream = urlopen(url)
    return datetime.strptime(stream.headers.get("Last-Modified"),
                             "%a, %d %b %Y %H:%M:%S %Z")


def list_available_organisms():
    """
    Return a list of all available GO organism codes.
    """
    source = urlopen("http://www.geneontology.org/gene-associations/").read()
    codes = re.findall("gene_association\.([a-zA-z0-9_]+?)\.gz", source.decode('utf-8'))
    return sorted(set(codes))


def sf_org_mtime(org_code):
    info = sf_server.info("GO", "gene_association.{}.tar.gz".format(org_code))
    if info:
        return info_date_time_parse(info["datetime"])
    else:
        # No file on the server, new file will be created
        return info_date_time_parse("2000-01-01 00:00:0.0")


def web_org_mtime(org_code):
    return http_last_modified(
        "http://www.geneontology.org/gene-associations/gene_association.{}.gz"
        .format(org_code))


def sf_ontology_mtime():
    info = sf_server.info("GO", "gene_ontology_edit.obo.tar.gz")
    if info:
        return info_date_time_parse(info["datetime"])
    else:
        # No file on the server, new file will be created
        return info_date_time_parse("2000-01-01 00:00:0.0")


def web_ontology_mtime():
    return http_last_modified(
        "http://www.geneontology.org/ontology/gene_ontology.obo")


def download_annotations(org, fpath):
    file_onserver = 'gene_association.' + org + '.gz'
    stream = urlopen('http://www.geneontology.org/gene-associations/' + file_onserver)

    compressed_file = io.BytesIO(stream.read())
    decompressed_file = gzip.GzipFile(fileobj=compressed_file)

    with open(download_path + "/gene_association." + org, 'wb') as outfile:
        outfile.write(decompressed_file.read())

    annos = obiGO.Annotations(outfile.name, genematcher=obiGene.GMDirect())

    with open(os.path.join(download_path, "gene_names.pickle"), "wb") as genes_f:
        pickle.dump(annos.gene_names, genes_f)

    with tarfile.open(fpath, 'w:gz') as tar_file:
        tar_file.add(outfile.name, 'gene_association')
        tar_file.add(genes_f.name, 'gene_names.pickle')


def download_ontology(fpath):
    download_file = 'gene_ontology_edit.obo'
    stream = urlopen('http://www.geneontology.org/ontology/' + download_file).read()
    temp_file_path = os.path.join(download_path, download_file)

    # temp files needed
    with open(temp_file_path, 'w+') as temp_file:
        temp_file.write(stream.decode())
        #shutil.copyfileobj(stream, temp_file)

    with tarfile.open(fpath, 'w:gz') as tar_file:
        tar_file.add(temp_file_path, arcname=download_file)


if web_ontology_mtime() > sf_ontology_mtime():
    FILENAME = 'gene_ontology_edit.obo.tar.gz'
    TITLE = 'Gene Ontology (GO)'
    TAGS = ['gene', 'ontology', 'GO', 'essential']
    VERSION = obiGO.Ontology.version
    file_path = os.path.join(domain_path, FILENAME)

    print("donwloading ontology...")
    download_ontology(file_path)

    create_info_file(file_path, title=TITLE, tags=TAGS, uncompressed=uncompressedSize(file_path),
                     compression='tar.gz', version=VERSION)


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

print(list_available_organisms())
for org in list_available_organisms():
    FILENAME = 'gene_association.' + org + '.tar.gz'
    FILE_PATH = os.path.join(domain_path, FILENAME)

    if org in exclude or org not in commonOrgs:
        continue

    print("Query", org)
    if 1 and web_org_mtime(org) <= sf_org_mtime(org):
        print('Update skiped, we have the latest version from geneontology.org\n')
        # Skip update
        continue

    print("Updating", org)

    download_annotations(org, FILE_PATH)

    # Load the annotations to test them and collect all taxon ids from them
    gene_association = os.path.join(download_path,  "gene_association." + org)
    a = obiGO.Annotations(gene_association, genematcher=obiGene.GMDirect())
    taxons = set([ann.Taxon for ann in a.annotations])

    # exclude taxons with cardinality 2
    taxons = [tax for tax in taxons if "|" not in tax]
    for tax in taxons:
        taxid = tax.split(":", 1)[-1]
        updatedTaxonomy[taxid].add(org)
    del a

    orgName = obiTaxonomy.name(commonOrgs[org])
    taxid = obiTaxonomy.taxname_to_taxid(orgName)

    TITLE = "GO Annotations for " + orgName
    TAGS = ["gene", "annotation", "ontology", "GO", orgName] + \
           (["essential"] if org in essentialOrgs else []) + obiTaxonomy.shortname(taxid)

    ORGANISM = orgName
    VERSION = obiGO.Annotations.version

    create_info_file(FILE_PATH, title=TITLE, tags=TAGS, version=VERSION,
                     uncompressed=uncompressedSize(FILE_PATH), compression='tar.gz')

try:
    with open(sf_local.localpath_download("GO", "taxonomy.pickle"), "rb") as f:
        tax = pickle.load(f)
except FileNotFoundError as e:
    tax = {}

# Upload taxonomy if any differences in the updated taxonomy
tax_path = os.path.join(domain_path, 'taxonomy.pickle')
TITLE = 'GO taxon IDs'
TAGS = ["GO", "taxon", "organism", "essential"]
VERSION = obiGO.Taxonomy.version

if any(tax.get(key, set()) != updatedTaxonomy.get(key, set()) for key in set(updatedTaxonomy)):
    print("takonomy was updated!")
    tax.update(updatedTaxonomy)

    with open(tax_path, "wb") as f:
        pickle.dump(tax, f)

    create_info_file(tax_path, title=TITLE, tags=TAGS, version=VERSION)

helper = SyncHelper(DOMAIN, GOTest)
helper.run_tests()
helper.sync_files()

helper.remove_update_folder()
