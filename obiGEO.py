import orange
import re
import gzip
import os.path
import orngServerFiles
import obiData
import orngMisc
import obiTaxonomy
import cPickle

def variableMean(x):
    vs = [v for v in x if v and v<>"?"]
    if len(vs) == 0:
        return "?"
    return sum(vs)/len(vs)

p_assign = re.compile(" = (.*$)")
p_tagvalue = re.compile("![a-z]*_([a-z_]*) = (.*)$")    
tagvalue = lambda x: p_tagvalue.search(x).groups()

DOMAIN = "GEO"
GDS_INFO_FILENAME = "gds_info.pickled"
FTP_NCBI = "ftp.ncbi.nih.gov"
FTP_DIR = "pub/geo/DATA/SOFT/GDS/"

class GDSInfo:
    def __init__(self):
        path = orngServerFiles.localpath(DOMAIN, GDS_INFO_FILENAME)
        if not os.path.exists(path):
            orngServerFiles.download(DOMAIN, GDS_INFO_FILENAME)
        f = file(path, "rb")
        self.info, self.excluded = cPickle.load(f)
    def keys(self): return self.info.keys()
    def items(self): return self.info.items()
    def values(self): return self.info.values()
    def clear(self): return self.info.clear()
    def __getitem__(self, key): return self.info[key]
    def __setitem__(self, key, item): self.info[key] = item
    def __len__(self): return len(self.info)
    

class GeneData:
    """Store mapping between spot id and gene."""
    def __init__(self, spot_id, gene_name, d):
        self.spot_id = spot_id
        self.gene_name = gene_name
        self.data = d

class GDS():
    """GEO DataSet class: read GEO datasets and convert them to ExampleTable."""
    def __init__(self, gdsname, verbose=False, force_download=False, cache=True):
        self.gdsname = gdsname
        self.verbose = verbose
        self.force_download = force_download
        self.filename = orngServerFiles.localpath(DOMAIN, self.gdsname + ".soft.gz")
        self._getinfo() # to get the info
        taxid = obiTaxonomy.search(self.info["sample_organism"], exact=True)
        self.info["taxid"] = taxid[0] if len(taxid)==1 else None
        self._getspotmap() # to get gene->spot and spot->gene mapping
        self.genes = self.gene2spots.keys()
        self.info["gene_count"] = len(self.genes)
        self.gdsdata = None
        self.data = None
        
    def _download(self):
        """Download GDS data file if not in local cache or forced download requested."""
        localpath = orngServerFiles.localpath(DOMAIN)
        if self.force_download or not os.path.exists(self.filename):
            ftp = obiData.FtpDownloader(FTP_NCBI, localpath, FTP_DIR)
            ftp.retrieve(self.gdsname + ".soft.gz", progressCallback=orngMisc.ConsoleProgressBar()
                         if self.verbose else None)
            if self.verbose:
                # print "Downloading %s" % self.gdsname
                print

    def _getinfo(self):
        """Parse GDS data file and return a dictionary with info."""
        getstate = lambda x: x.split(" ")[0][1:] 
        getid = lambda x: x.rstrip().split(" ")[2]
        self._download()
        f = gzip.open(self.filename)
        state = None; previous_state = None
    
        info = {"subsets" : []}
        subset = None
    
        # GDS information part
        for line in f:
            if line[0] == "^":
                previous_state = state; state = getstate(line)
                if state == "SUBSET":
                    if subset:
                        info["subsets"] += [subset]
                    subset = {"id" : getid(line)}
                if state == "DATASET":
                    info["dataset_id"] = getid(line)
                continue
            if state == "DATASET":
                if previous_state == "DATABASE":
                    tag, value = tagvalue(line)
                    info[tag] = value
                else:
                    if subset:
                        info["subsets"] += [subset]
                    break
            if state == "SUBSET":
                tag, value = tagvalue(line)
                if tag == "description" or tag == "type":
                    subset[tag] = value
                if tag == "sample_id":
                    subset[tag] = value.split(",")
        for t,v in info.items():
            if "count" in t:
                info[t] = int(v)
                
        # sample information
        state = None
        for line in f:
            if state == "header":
                info["samples"] = line.rstrip().split("\t")[2:]
                break
            if line.startswith("!dataset_table_begin"):
                state = "header"

        self.info = info

    def _getspotmap(self, include_spots=None):
        """Return gene to spot and spot to genes mapings."""
        f = gzip.open(self.filename)
        for line in f:
            if line.startswith("!dataset_table_begin"):
                break
        f.readline() # skip header
        spot2gene = {}
        gene2spots = {}
        for line in f:
            if line.startswith("!dataset_table_end"):
                break
            spot, gene = line.rstrip().split("\t")[0:2]
            if include_spots and (spot not in include_spots):
                continue 
            spot2gene[spot] = gene
            gene2spots[gene] = gene2spots.get(gene, []) + [spot]
    
        self.spot2gene = spot2gene
        self.gene2spots = gene2spots
        
    def sample_to_class(self, classes=None, missing_class_value=None):
        """Return class values for GDS samples."""
        samples = self.info["samples"]
        subsets = self.info["subsets"]
        if classes:
            subsets = [ss for ss in subsets if ss["id"] in classes or ss["description"] in classes]
        sample2class = {}
        for sample in samples:
            classval = [ss["description"] for ss in subsets if sample in ss["sample_id"]]
            classval.sort()
            sample2class[sample] = "|".join(classval) if classval else missing_class_value
        return sample2class

    def _parse_soft(self, filter_unknown=None):
        """Parse GDS data, returns data dictionary."""
        f = gzip.open(self.filename)
        mfloat = lambda x: float(x) if x<>'null' else '?'
    
        data = {}
        # find the start of the data part
        for line in f:
            if line.startswith("!dataset_table_begin"):
                break
        f.readline()
        
        # load data
        for line in f:
            if line.startswith("!dataset_table_end"):
                break
            d = line.rstrip().split("\t")
            if filter_unknown and (float(d[2:].count('null')) / len(d[2:]) > filter_unknown):
                continue 
            data[d[0]] = GeneData(d[0], d[1], [mfloat(v) for v in d[2:]])
        
        self.gdsdata = data
    
    def _to_ExampleTable(self, report_genes=True, merge_function=variableMean,
                                classes=None, missing_class_value=None, transpose=False):
        """Convert parsed GEO format to orange, save by genes or by spots."""
    
        sample2class = self.sample_to_class(classes, missing_class_value)
        cvalues = list(set(sample2class.values()))
        if None in cvalues:
            cvalues.remove(None)
            
        classvar = orange.EnumVariable(name="outcome", values=cvalues)
        orng_data = []
        if transpose: # genes as attributes, samples as data instances
            if report_genes: # save by genes
                atts = [orange.FloatVariable(name=gene) for gene in self.genes]
                domain = orange.Domain(atts, classvar)
                for (i, sampleid) in enumerate(self.info["samples"]):
                    vals = [merge_function([self.gdsdata[spot].data[i] \
                            for spot in self.gene2spots[gene]]) for gene in self.genes]
                    orng_data.append(vals + [sample2class[sampleid]])
                
            else: # save by spots
                spots = self.spot2gene.keys()
                atts = [orange.FloatVariable(name=id) for id in spots]
                domain = orange.Domain(atts, classvar)
                for (i, sampleid) in enumerate(info["samples"]):
                    orng_data.append([data[spot].data[i] for spot in spots] + [sample2class[sampleid]])
    
            return orange.ExampleTable(domain, orng_data)
    
        else: # samples as attributes, genes as data instances
            atts = [orange.FloatVariable(name=ss) for ss in self.info["samples"]]
            domain  = orange.Domain(atts, False)
    
            if report_genes: # save by genes
                domain.addmeta(orange.newmetaid(), orange.StringVariable("gene"))
                for g in self.genes:
                    orng_data.append(map(lambda *x: merge_function(x),
                                         *[self.gdsdata[spot].data for spot in self.gene2spots[g]]))
            else: # save by spots
                domain.addmeta(orange.newmetaid(), orange.StringVariable("spot"))
                spots = self.spot2gene.keys()
                orng_data = [self.gdsdata[spot].data for spot in spots]
    
            data = orange.ExampleTable(domain, orng_data)
            if report_genes:
                for i, g in enumerate(self.genes):
                    data[i]["gene"] = g
            else:
                for i, s in enumerate(spots):
                    data[i]["spot"] = s
            return data
        
    def getdata(self, report_genes=True, merge_function=variableMean,
                 classes=None, missing_class_value=None,
                 transpose=False, filter_unknown=None):
        """Load GDS data and returns a corresponding orange data set,
        spot<->gene mappings and subset info."""
        if self.verbose: print "Reading data ..."
        if not self.gdsdata:
            self._parse_soft(filter_unknown = filter_unknown)
        if filter_unknown:
            # some spots were filtered out, need to revise spot<>gene mappings
            self._getspotmap(include_spots=set(self.gdsdata.keys()))
        if self.verbose: print "Converting to example table ..."
        self.data = self._to_ExampleTable(merge_function=merge_function,
                                                 classes=classes, transpose=transpose,
                                                 report_genes=report_genes)
        return self.data

    def __str__(self):
        return "%s (%s), samples=%s, features=%s, genes=%s, subsets=%d" % \
              (self.info["dataset_id"],
               self.info["sample_organism"],
               self.info['sample_count'],
               self.info['feature_count'],
               self.info['gene_count'],
               len(self.info['subsets']),
               )

if __name__ == "__main__":
    gds = GDS("GDS10")
    data = gds.getdata(report_genes=True, transpose=False)
#    data = GDS("GDS10", report_genes=True, transpose=False)
