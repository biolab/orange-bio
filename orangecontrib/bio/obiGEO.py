from __future__ import absolute_import

import cPickle, gzip, os.path, re

import orange
from Orange.orng import orngMisc, orngServerFiles

from . import obiData, obiTaxonomy

def spots_mean(x):
    vs = [v for v in x if v and v<>"?"]
    if len(vs) == 0: return "?"
    return sum(vs)/len(vs)
def spots_median(x):
    vs = [v for v in x if v and v<>"?"]
    if len(vs) == 0: return "?"
    if len(vs) % 2:
        return sorted(vs)/(len(vs)/2)
    else:
        z = sorted(x)
        return (z[len(vs)/2-1] + z[len(vs)/2]) / 2. 
    return sum(vs)/len(vs)
def spots_min(x):
    vs = [v for v in x if v and v<>"?"]
    if len(vs) == 0: return "?"
    return min(vs)/len(vs)
def spots_max(x):
    vs = [v for v in x if v and v<>"?"]
    if len(vs) == 0: return "?"
    return max(vs)/len(vs)

p_assign = re.compile(" = (.*$)")
p_tagvalue = re.compile("![a-z]*_([a-z_]*) = (.*)$")    
tagvalue = lambda x: p_tagvalue.search(x).groups()

DOMAIN = "GEO"
GDS_INFO_FILENAME = "gds_info.pickled"
FTP_NCBI = "ftp.ncbi.nih.gov"
FTP_DIR = "pub/geo/DATA/SOFT/GDS/"

class GDSInfo:
    def __init__(self, force_update=False):
        path = orngServerFiles.localpath(DOMAIN, GDS_INFO_FILENAME)
        if not os.path.exists(path) or force_update:
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
    def __iter__(self): return iter(self.info)
    def __contains__(self, key): return key in self.info
    

class GeneData:
    """Store mapping between spot id and gene."""
    def __init__(self, spot_id, gene_name, d):
        self.spot_id = spot_id
        self.gene_name = gene_name
        self.data = d

class GDS():
    """GEO DataSet class: read GEO datasets and convert them to ExampleTable."""
    def __init__(self, gdsname, verbose=False, force_download=False):
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
        
    def sample_annotations(self, sample_type=None):
        """Return a dictionary with sample annotation."""
        annotation = {}
        for info in self.info["subsets"]:
            if sample_type and info["type"]<>sample_type:
                continue
            for id in info["sample_id"]:
                annotation.setdefault(id, {})[info["type"]]=info["description"]
        return annotation

    def sample_to_class(self, sample_type=None, missing_class_value=None):
        """Return class values for GDS samples."""
        annotations = self.sample_annotations(sample_type)
        return dict([(sample, "|".join([a for t,a in ann.items()])) for sample, ann in annotations.items()])
        
    def sample_types(self):
        """Return a set of sample types."""
        return set([info["type"] for info in self.info["subsets"]])
    
    def _parse_soft(self, remove_unknown=None):
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
            if remove_unknown and (float(d[2:].count('null')) / len(d[2:]) > remove_unknown):
#                print "ZZZ", float(d[2:].count('null')) / len(d[2:]), d[2:]
                continue 
            data[d[0]] = GeneData(d[0], d[1], [mfloat(v) for v in d[2:]])
        self.gdsdata = data
    
    def _to_ExampleTable(self, report_genes=True, merge_function=spots_mean,
                                sample_type=None, missing_class_value=None, transpose=False):
        """Convert parsed GEO format to orange, save by genes or by spots."""
        orng_data = []
        if transpose: # samples in rows
            sample2class = self.sample_to_class(sample_type, missing_class_value)
            cvalues = list(set(sample2class.values()))
            if None in cvalues:
                cvalues.remove(None)
            classvar = orange.EnumVariable(name=sample_type or "class", values=cvalues)
            if report_genes: # save by genes
                atts = [orange.FloatVariable(name=gene) for gene in self.gene2spots.keys()]
                domain = orange.Domain(atts, classvar)
                for (i, sampleid) in enumerate(self.info["samples"]):
                    vals = [merge_function([self.gdsdata[spot].data[i] \
                            for spot in self.gene2spots[gene]]) for gene in self.gene2spots.keys()]
                    orng_data.append(vals + [sample2class[sampleid]])
                
            else: # save by spots
                spots = self.spot2gene.keys()
                atts = [orange.FloatVariable(name=id) for id in spots]
                domain = orange.Domain(atts, classvar)
                for (i, sampleid) in enumerate(self.info["samples"]):
                    orng_data.append([self.gdsdata[spot].data[i] for spot in spots] + [sample2class[sampleid]])
            if missing_class_value == None:
                orng_data = [example for example in orng_data if example[-1] != None]
    
            return orange.ExampleTable(domain, orng_data)
    
        else: # genes in rows
            annotations = self.sample_annotations(sample_type)
            atts = [orange.FloatVariable(name=ss) for ss in self.info["samples"]]
            for i, a in enumerate(atts):
                a.setattr("attributes", annotations[self.info["samples"][i]])
            domain  = orange.Domain(atts, False)
    
            if report_genes: # save by genes
                domain.addmeta(orange.newmetaid(), orange.StringVariable("gene"))
                for g in self.gene2spots.keys():
                    orng_data.append(map(lambda *x: merge_function(x),
                                         *[self.gdsdata[spot].data for spot in self.gene2spots[g]]))
            else: # save by spots
                domain.addmeta(orange.newmetaid(), orange.StringVariable("spot"))
                spots = self.spot2gene.keys()
                orng_data = [self.gdsdata[spot].data for spot in spots]
    
            data = orange.ExampleTable(domain, orng_data)

            if report_genes:
                for i, g in enumerate(self.gene2spots.keys()):
                    data[i]["gene"] = g
            else:
                for i, s in enumerate(spots):
                    data[i]["spot"] = s
            return data
        
    def getdata(self, report_genes=True, merge_function=spots_mean,
                 sample_type=None, missing_class_value=None,
                 transpose=False, remove_unknown=None):
        """Load GDS data and returns a corresponding orange data set,
        spot<->gene mappings and subset info."""
        if self.verbose: print "Reading data ..."
#        if not self.gdsdata:
        self._parse_soft(remove_unknown = remove_unknown)
#        if remove_unknown:
            # some spots were filtered out, need to revise spot<>gene mappings
        self._getspotmap(include_spots=set(self.gdsdata.keys()))
        if self.verbose: print "Converting to example table ..."
        self.data = self._to_ExampleTable(merge_function=merge_function,
                                          sample_type=sample_type, transpose=transpose,
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


def _float_or_na(x):
    if x.isSpecial():
        return "?"
    return float(x)

def transpose_class_to_labels(data, attcol="sample"):
    """Converts data with genes as attributes to data with genes in rows."""
    if attcol in [v.name for v in data.domain.getmetas().values()]:
        atts = [orange.FloatVariable(str(d[attcol])) for d in data]
    else:
        atts = [orange.FloatVariable("S%d" % i) for i in range(len(data))]
    for i, d in enumerate(data):
        atts[i].setattr("class", str(d.getclass()))
    domain = orange.Domain(atts, False)
    
    newdata = []
    for a in data.domain.attributes:
        newdata.append([_float_or_na(d[a]) for d in data])

    gene = orange.StringVariable("gene")
    id = orange.newmetaid()
    new = orange.ExampleTable(domain, newdata)
    new.domain.addmeta(id, gene)
    for i, d in enumerate(new):
        d[gene] = data.domain.attributes[i].name

    return new

def transpose_labels_to_class(data, class_label=None, gene_label="gene"):
    """Converts data with genes in rows to data with genes as attributes."""
    # if no class_label (attribute type) given, guess it from the data
    if not class_label:
        l = []
        for a in data.domain.attributes:
            l.extend(a.attributes.keys())
        l = list(set(l))
        class_label = l[0]
        if len(set(l)) > 1:
            import warnings
            warnings.warn("More than single attribute label types (%s), took %s"
                          % (", ".join(l), class_label))

    if gene_label in [v.name for v in data.domain.getmetas().values()]:
        atts = [orange.FloatVariable(str(d[gene_label])) for d in data]
    else:
        atts = [orange.FloatVariable("A%d" % i) for i in range(len(data))]
        
    classvalues = list(set([a.attributes[class_label] for a in data.domain.attributes]))
    
    if all(map(lambda x: isinstance(x, (int, long, float, complex)), classvalues)):
        classvar = orange.FloatVariable(class_label)
    else:
        classvar = orange.EnumVariable(class_label, values=classvalues)
        
    domain = orange.Domain(atts + [classvar])
    
    newdata = []
    for a in data.domain.attributes:
        newdata.append([_float_or_na(d[a]) for d in data] + [a.attributes[class_label]])

    sample = orange.StringVariable("sample")
    id = orange.newmetaid()
    new = orange.ExampleTable(domain, newdata)
    new.domain.addmeta(id, sample)
    for i, d in enumerate(new):
        d[sample] = data.domain.attributes[i].name

    return new

def transpose(data):
    """Transposes data matrix, converts class information to attribute label and back"""
    if data.domain.classVar:
        return transpose_class_to_labels(data)
    else:
        return transpose_labels_to_class(data)
