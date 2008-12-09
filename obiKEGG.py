"""obiKEGG is an interface to Kyoto Encyclopedia of Genes and Genomes (http://www.genome.jp/kegg/) that allows easy access to KEGG pathway and genes data.

"""
try:
    import Image, ImageDraw, ImageMath
except:
    pass

import cStringIO
import math
import time
import os, sys, tarfile
import re

import obiData
import obiProb

from cPickle import load, loads, dump
from collections import defaultdict

try:
    import orngServerFiles
    default_database_path = orngServerFiles.localpath("KEGG")
except:
    default_database_path = os.path.join((os.path.split(__file__)[0] or "."), "data//kegg//")
    
try:
    os.mkdir(default_database_path)
except:
    pass

base_ftp_path = "ftp://ftp.genome.jp/pub/kegg/"

forceUpdate = False
    
class KEGGInterface(object):
    def __init__(self):
        try:
            from SOAPpy import WSDL
            wsdl = 'http://soap.genome.jp/KEGG.wsdl'
            self.serv = WSDL.Proxy(wsdl)
        except:
            self.serv = None
        
    def list_organisms(self):
        return dict(self.serv.list_organisms())

    def list_pathways(self, org="map"):
        return dict(self.serv.list_pathways(org))

    def get_pathways_by_genes(self, genes_list):
        return self.serv.get_pathways_by_genes(genes_list)

    def get_pathways_by_enzymes(self, enzyme_list):
        return self.serv.get_pathways_by_enzymes(enzyme_list)

    def get_pathways_by_compounds(self, compound_list):
        return self.serv.get_pathways_by_compounds(compound_list)

    def get_linked_pathways(self, pathway_id):
        return self.serv.get_linked_pathways(pathway_id)

    def get_genes_by_pathway(self, pathway_id):
        return self.serv.get_genes_by_pathway(pathway_id)

    def get_genes_by_organism(self, org, offset=1, limit=-1):
        if limit==-1:
            limit = self.get_number_of_genes_by_organism(org)
        return self.serv.get_genes_by_organism(org, offset, limit)

    def get_number_of_genes_by_organism(self, org):
        return self.serv.get_number_of_genes_by_organism(org)

    def get_enzymes_by_pathway(self, pathway_id):
        return self.serv.get_enzymes_by_pathway(pathway_id)
    
    def get_enzymes_by_compound(self, compound_id):
        return self.serv.get_enzymes_by_compound(compound_id)

    def get_compounds_by_enzyme(self, enzyme_id):
        return self.serv.get_compounds_by_enzyme(enzyme_id)

    def get_genes_by_enzyme(self, enzyme_id, org):
        return self.serv.get_genes_by_enzyme(enzyme_id, org)

    def get_enzymes_by_gene(self, gene_id):
        return self.serv.get_enzymes_by_gene(gene_id)

    def get_colored_pathway_image(self, pathway_id, objects):
        import htmllib
        class HTMLImageCollector(htmllib.HTMLParser):
            def __init__(self):
                self.images = []
            def handle_image(self, source, *args):
                self.images.append(source)

        obj_list = [ob[0] for ob in objects]
        fg_color = ["blue"]*len(objects)
        bg_color = [ob[1] for ob in objects]
        try:
            url = self.serv.get_html_of_colored_pathway_by_objects(pathway_id, obj_list, fg_color, bg_color)
            sitestr = urllib.urlopen(url).read()
            parser = HTMLImageCollector()
            parser.feed(sitestr)
            url = parser.images[-1]
            f = urllib.urlopen(url)
            imgstr = f.read()
            image = Image.open(cStringIO.StringIO(imgstr))
        except Exception, ex:
            print ex
            raise ValueError(pathway_id)
        return image
        
    def get_pathway_image(self, pathway_id):
        return self.get_pathway_image_ex(pathway_id[5:-5], pathway_id[-5:])

    def get_pathway_image_ex(self, org, pathway_num):
        filename = org+pathway_num+".gif"
        if org=="map":
            dir = "map/"
        else:
            dir = "organisms/"+org+"/"
        try:
            url = base_ftp_path+"pathway/"+dir+filename
            f = urllib.urlopen(url)
            imgstr = f.read()
            image = Image.open(cStringIO.StringIO(imgstr))
        except Exception, ex:
            print ex
            raise ValueError(org+pathway_num)
        return image

    def get_unique_gene_ids(self, org, genes):
        return genes, [], []

def _collect(list, func=None):
    return reduce(lambda a,b: a + (func and func(b) or b), list, [])

def _rel_dir(pathway_id):
    if "map" in pathway_id:
        return "pathway/map/"
    else:
        return "pathway/organisms/"+pathway_id.split(":")[-1][:-5]+"/"

def _tabspliter(file):
    return [tuple(l.split("\t")) for t in file.readlines()]

class DBEntry(object):
    cache = []
    def __init__(self, text):
        self.text = text
        self.section = {}
        self.parse(text)
        
    def parse(self, text):
        currsection = ""
        title = ""
        for line in text.split("\n"):
            if line.startswith(" "):
                currsection = currsection + line + "\n"
            elif line.split():
                if title:
                    self.section[title] = currsection
                title = line.split()[0]
                currsection = line[len(title):] + "\n"
        self.section[title] = currsection

    def get_by_lines(self, title):        
        if title in self.section:
            return [s.strip() for s in self.section[title].split("\n")]
        else:
            return []

    def get_by_list(self, title):
        if title in self.section:
            return self.section[title].split()
        else:
            return []

    def get_subsections(self, title):
        li = self.get_by_list(title)
        d = []
        for s in li:
            if s.endswith(":"):
                d.append((s[:-1], []))
            else:
                d[-1][1].append(s)
        return d

    def get_string(self, title):
        return " ".join(self.get_by_list(title))
        
class DBEnzymeEntry(DBEntry):
    cache = ["genes", "pathways", "name"]
    def get_genes(self):
        d = dict(self.get_subsections("GENES"))
        return _collect(d.items(), lambda (org,genes):[org.lower()+":"+g.split("(")[0] for g in genes])
    def get_pathways(self):
        d = self.get_by_lines("PATHWAY")
        return ["path:"+line.split()[1] for line in d if len(line.split())>=2]
    def get_name(self):
        e = self.get_by_list("ENTRY")
        return e and e[0].lower()+":"+e[1] or "unknown"
    
class DBCompoundEntry(DBEntry):
    cache = ["pathways", "enzymes", "name"]
    def get_pathways(self):
        d = self.get_by_lines("PATHWAY")
        return ["path:"+line.split()[1] for line in d if len(line.split())>=2]
    def get_enzymes(self):
        d = self.get_by_list("ENZYME")
        return ["ec:"+s.strip() for s in d]
    def get_name(self):
        e = self.get_by_list("ENTRY")
        return e and "cpd:"+e[0] or "unknown"

class DBGeneEntry(DBEntry):
    cache = ["name", "enzymes", "alt_names", "pathways", "db_links"]
    def get_name(self):
        e = self.get_by_list("ENTRY")
        return e and e[0].strip() or "unknown"

    def get_enzymes(self):
        import re
        s = self.get_string("DEFINITION")
        s = re.findall("\[EC:([1-9]+|-)\.([1-9]+|-)\.([1-9]+|-)\.([1-9]+|-)\]", s)
        return map(lambda t:"ec:"+".".join(t), s)
        
    def get_alt_names(self):
        lines = self.get_by_lines("DBLINKS")
        return reduce(list.__add__, [line.split()[1:] for line in lines if len(line.split())>=2], []) + [n.strip(",\t \n") for n in self.get_by_list("NAME")] +[self.get_name()]

    def get_db_links(self):
        lines = self.get_by_lines("DBLINKS")
        return dict([(line.split()[0].rstrip(":"), line.split()[1:]) for line in lines if len(line.split())>=2])
    
    def get_pathways(self):
        lines = self.get_by_lines("PATHWAY")
        return ["path:"+line.split()[1] for line in lines if len(line.split())>=2]

class DBOrganismEntry(DBEntry):
    cache = ["name", "annotation_type", "taxid"]
    def get_name(self):
        e = self.get_by_list("ENTRY")
        return e and e[0].strip() or "unknown"
    
    def get_annotation_type(self):
        return self.get_string("ANNOTATION")

    def get_taxid(self):
        return self.get_by_lines("TAXONOMY")[0].split(":")[-1].strip()

class DBEntryWrapper(object):
    def __init__(self, wrapped):
        for name in wrapped.cache:
            setattr(self, name, getattr(wrapped, "get_"+name)())
    def __getattr__(self, name):
        if name.startswith("get_") and name[4:] in self.__dict__:
            return lambda :self.__dict__[name[4:]]
        else:
            raise AttributeError(name)

class GenesDatabaseProxy(defaultdict):
    def __init__(self, interface, *args, **argskw):
        defaultdict.__init__(self, lambda :None, *args, **argskw)
        self.interface = interface
    def __missing__(self, key):
        self.__setattr__(key, self.interface._load_gene_database(key))
        return self.get(key)
    
class KEGGInterfaceLocal(object):
    _instanceCache = {}
    def __init__(self, update=False, local_database_path=None, download_progress_callback=None):
        self.local_database_path = local_database_path or default_database_path
        self._build_index()
        self.update = update
        self.download_progress_callback = download_progress_callback
        self._gene_alias = {}
        self._gene_alias_conflicting = {}
        self._from_pathway_to_genes = defaultdict(set)
        self._filenames = {"_enzymes":"ligand/enzyme/_enzymes.pickle",
                           "_from_gene_to_enzymes":"ligand/enzyme/_from_gene_to_enzymes.pickle",
                           "_compounds":"ligand/compound/_compounds.pickle",
                           "_from_enzyme_to_compounds":"ligand/compound/_from_enzyme_to_compounds.pickle",
                           "_genome":"genes/_genome.pickle"}
        self.downloader = obiData.FtpDownloader("ftp.genome.jp", self.local_database_path, "/pub/kegg/", numOfThreads=10)
        self._instanceCache[local_database_path] = self

    def _build_index(self):
        tarfiles = [name for name in os.listdir(self.local_database_path) if name.endswith(".tar.gz") and os.path.isfile(os.path.join(self.local_database_path, name))]
        tardirs = [name for name in os.listdir(self.local_database_path) if name.endswith(".tar.gz") and os.path.isdir(os.path.join(self.local_database_path, name))]
        self.openTarFiles = {}
        self.cachedExtractedFiles = {}
        self.inTarfileDict = dict([(os.path.normpath(name), filename) for filename in tarfiles for name in tarfile.open(os.path.join(self.local_database_path, filename)).getnames()])
        self.inTardirDict = {}
        for dir in tardirs:
             ls = list(os.walk(os.path.join(self.local_database_path, dir)))
             for path, dirs, files in ls:
                 self.inTardirDict.update(dict([(os.path.normpath(os.path.join(path, file).replace(dir, "", 1)), os.path.normpath(os.path.join(path, file))) for file in files]))

    def _rel_org_dir(self, org):
        if org=="map":
            return "map/"
        elif self._genome[org].get_annotation_type()=="manual": ##            return "pathway/organisms/"+pathway_id.split(":")[-1][:-5]+"/"
            return "organisms/"+org+"/"
        elif len(org)==4 and org.startswith("e"):
            return "organisms_est/"+org+"/"
        else:
            return "organisms_kaas/"+org+"/"
        
    def _rel_pathway_dir(self, pathway_id):
        return "pathway/"+self._rel_org_dir(pathway_id.split(":")[-1][:-5])

    def _pathway_from(self, pathway_id):
        return "kegg_reference.tar.gz" if "map" in pathway_id else "kegg_organism_%s.tar.gz" % pathway_id.split(":")[-1][:-5]

    def download_organism_data(self, org):
        files = ["pathway/map_title.tab", "genes/taxonomy", "genes/genome"]
        self.downloader.massRetrieve(files, update=self.update)
        rel_path = "pathway/"+self._rel_org_dir(org)
##        files = [rel_path+org+"_gene_map.tab", "pathway/map_title.tab", "genes/taxonomy"]
        self.downloader.retrieve(rel_path+org+"_gene_map.tab", update=self.update)
        file = self._retrieve(rel_path+org+"_gene_map.tab")
        pathway_nums = set(reduce(lambda a,b: a + b.split()[1:], file.readlines(), []))
        descr = dict(map(lambda line:tuple(line.strip().split("\t")), self._retrieve("pathway/map_title.tab").readlines()))
##        dump(descr, open(self.local_database_path+"list_pathways_map.pickle", "w"))
        ids = [org+num for num in pathway_nums]
        try:
            organisms = load(open(self.local_database_path+"list_organisms.pickle"))
        except:
            organisms = {}
        if org not in organisms:
            organisms[org] = self._taxonomy.get(org, "  ")[1]
##            dump(organisms, open(self.local_database_path+"list_organisms.pickle", "w"))
##        dump(dict([("path:"+org+num, descr[num]) for num in pathway_nums]), open(self.local_database_path+"list_pathways_"+org+".pickle","w"))
        
        ends = [".cpd", ".gene", ".gif", ".map", "_cpd.coord", "_gene.coord", ".conf"]
        files = [rel_path+id+ext for id in ids for ext in ends]
        self.downloader.massRetrieve(files, update=self.update, blocking=False)
        self.downloader.retrieve("genes/"+self._rel_org_dir(org)+self._taxonomy[org][0]+".ent", update=self.update, progressCallback=self.download_progress_callback)
        while not self.downloader.queue.empty():
            if self.download_progress_callback:
                self.download_progress_callback(min(100.0, 100.0*(float(len(files))-self.downloader.queue.qsize())/len(files)))
            time.sleep(0.1)

    def download_reference_data(self):
        rel_path = "pathway/map/"
        descr = dict(map(lambda line:tuple(line.strip().split("\t")), self._retrieve("pathway/map_title.tab").readlines()))
        ends = [".conf", ".gif", ".map"]
        files = [rel_path+"map"+pathNum+ext for pathNum in descr.keys() for ext in ends]
        self.downloader.massRetrieve(files, update=self.update, progressCallback=self.download_progress_callback)
        
    def __getattr__(self, name):
        if name=="_enzymes" or name=="_from_gene_to_enzymes" :
            self._load_enzyme_database()
            return getattr(self, name)
        elif name=="_compounds" or name=="_from_enzyme_to_compounds":
            self._load_compound_database()
            return getattr(self, name)
        elif name=="_genes":
            self._genes = GenesDatabaseProxy(self)
            return getattr(self, name)
            #self._load_genes_database()
        elif name=="_taxonomy" or name=="_genome":
            self._load_taxonomy()
            return getattr(self, name)
        else:
            raise AttributeError(name)

    def _load_pickled(self, filename=None, name=None, from_=None):
        if not filename and name:
##            return load(open(self.local_database_path+self._filenames[name]))
            return loads(self._retrieve(self._filenames[name], from_, mode="rb").read().replace("\r\n", "\n"))
        else:
##            return load(open(self.local_database_path+filename))
            return loads(self._retrieve(filename, from_, mode="rb").read().replace("\r\n", "\n"))

    def _dump_pickled(self, object, filename=None, name=None):
        if not  filename and name:
            dump(object, open(self.local_database_path+self._filenames[name], "wb"))
        else:
            dump(object, open(self.local_database_path+filename, "wb"))
    
    def _load_enzyme_database(self, from_=None):
        try:
            self._enzymes = self._load_pickled(name="_enzymes", from_=from_ or "kegg_enzyme_and_compounds.tar.gz")
        except Exception, ex:
            print ex
            enzymes = map(DBEnzymeEntry, filter(bool, self._retrieve("ligand/enzyme/enzyme", from_ or "kegg_enzyme_and_compounds.tar.gz").read().split("///\n")))
            self._enzymes = dict([(e.get_name(), DBEntryWrapper(e)) for e in enzymes])
            self._dump_pickled(self._enzymes, name="_enzymes")
        try:
            self._from_gene_to_enzymes = self._load_pickled(name="_from_gene_to_enzymes", from_=from_ or "kegg_enzyme_and_compounds.tar.gz")
        except Exception, ex:
            self._from_gene_to_enzymes = defaultdict(list)
            for id, e in self._enzymes.items():
                for g in e.get_genes():
                    self._from_gene_to_enzymes[g].append(id)
            self._dump_pickled(self._from_gene_to_enzymes, name="_from_gene_to_enzymes")
        
    def _load_compound_database(self, from_=None):
        try:
            self._compounds = self._load_pickled(name="_compounds", from_=from_ or "kegg_enzyme_and_compounds.tar.gz")
        except:
            compounds = map(DBCompoundEntry, filter(bool, self._retrieve("ligand/compound/compound", from_ or "kegg_enzyme_and_compounds.tar.gz").read().strip().split("///\n")))
            self._compounds = dict([(c.get_name(), DBEntryWrapper(c)) for c in compounds])
            self._dump_pickled(self._compounds, name="_compounds")
        try:
            self._from_enzyme_to_compounds = self._load_pickled(name="_from_enzyme_to_compounds", from_=from_ or "kegg_enzyme_and_compounds.tar.gz")
        except:
            self._from_enzyme_to_compounds = defaultdict(list)
            for id, c in self._compounds.items():
                for e in c.get_enzymes():
                    self._from_enzyme_to_compounds[e].append(id)
            self._dump_pickled(self._from_enzyme_to_compounds, name="_from_enzyme_to_compounds")

    def _load_gene_database(self, org, from_=None):
        rel_path = "genes/" + self._rel_org_dir(org)
        freshLoad = False
        try:
            self._genes[org] = self._load_pickled(rel_path + "_genes.pickle", from_=from_ or "kegg_organism_%s.tar.gz" % org)
        except Exception, ex:
##            print >> sys.stderr, ex
            genes = map(DBGeneEntry, filter(bool ,self._retrieve(rel_path+self._taxonomy[org][0] + ".ent", from_ or "kegg_organism_%s.tar.gz" % org).read().split("///\n")))
            self._genes[org] = dict([(org + ":" + g.get_name(), DBEntryWrapper(g)) for g in genes])
            self._dump_pickled(self._genes[org], rel_path + "_genes.pickle")
            freshLoad = True
        self._gene_alias[org] = {}
        self._gene_alias_conflicting[org] = set()
        for id, gene in self._genes[org].items():
            aliases = gene.get_alt_names()
            for alias in set(aliases):
                if alias in self._gene_alias[org]:
                    self._gene_alias_conflicting[org].add(alias)
                else:
                    self._gene_alias[org][alias] = id
            for p_id in gene.get_pathways():
                self._from_pathway_to_genes[p_id].add(id)
                    
        if freshLoad:
            try:
                dump(set(self._gene_alias[org].keys() + self._genes[org].keys()), open(self.local_database_path+org+"_genenames.pickle", "wb"))
            except Exception:
                pass
        return self._genes[org]

    def _load_taxonomy(self, from_=None):
        orgs = filter(lambda line:line.strip() and not line.startswith("#"), self._retrieve("genes/taxonomy", "kegg_taxonomy.tar.gz").readlines())
        d = dict([(line.split()[1].strip(), (line.split("\t")[-2].strip(), line.split("\t")[-1].strip())) for line in orgs])
        self._taxonomy = d
        try:
            self._genome = self._load_pickled(name = "_genome", from_=from_ or "kegg_taxonomy.tar.gz")
        except Exception:
            entrys = map(DBOrganismEntry, filter(bool, self._retrieve("genes/genome", from_ or "kegg_taxonomy.tar.gz").read().split("///\n")))
            self._genome = dict([(e.get_name(), DBEntryWrapper(e)) for e in entrys])
            try:
                self._dump_pickled(self._genome, name="_genome")
            except Exception:
                pass
        
    def _retrieve(self, filename, from_=None, mode="rb"):
        if forceUpdate == True or self.update == "Force update":
            self.downloader.retrieve(filename, update=self.update, progressCallback=self.download_progress_callback)
        if from_ and not os.path.exists(os.path.join(self.local_database_path, from_)):
            import orngServerFiles
            orngServerFiles.download("KEGG", from_)
            self._build_index()
        try:
            if from_ and os.path.isdir(os.path.join(self.local_database_path, from_)):
                return open(os.path.join(self.local_database_path, from_, filename), mode)
            else:
                return open(os.path.join(self.local_database_path, filename), mode)
        except Exception:
            if os.path.normpath(os.path.join(self.local_database_path, filename)) in self.inTardirDict:
                return open(self.inTardirDict[os.path.normpath(os.path.join(self.local_database_path, filename))], mode)
            elif os.path.normpath(filename) in self.inTarfileDict:
                tarFileName = self.inTarfileDict[os.path.normpath(filename)]
                if tarFileName not in self.openTarFiles:
##                    print "opening tar file " + tarFileName
                    self.openTarFiles[tarFileName] = tarfile.open(os.path.join(self.local_database_path, tarFileName))
                if (tarFileName, os.path.normpath(filename)) not in self.cachedExtractedFiles:
##                    print "extracting: ", filename
                    data = self.openTarFiles[tarFileName].extractfile(filename).read()
                    self.cachedExtractedFiles[tarFileName, os.path.normpath(filename)] = data
                return cStringIO.StringIO(self.cachedExtractedFiles[tarFileName, os.path.normpath(filename)])
##                return tarfile.open(os.path.join(self.local_database_path, self.inTarfileDict[os.path.normpath(filename)])).extractfile(filename)
            else:
                raise
    
    def list_organisms(self):
        return dict([(key, value[1]) for key, value in self._taxonomy.items()])
        #return load(open(self.local_database_path+"list_organisms.pickle"))
    
    def list_pathways(self, org="map"):
        if org=="map":
            r = map(lambda line:tuple(line.strip().split("\t")), self._retrieve("pathway/map_title.tab", "kegg_reference.tar.gz").readlines())
            return dict([("path:map"+p, desc) for p, desc in r])
        else:
            rel_path = "pathway/"+self._rel_org_dir(org)
            ids = set(_collect(self._retrieve(rel_path+org+"_gene_map.tab", "kegg_organism_%s.tar.gz" % org).readlines(), lambda line:line.split()[1:]))
            pathways = self.list_pathways("map")
            return dict([("path:"+org+id, pathways["path:map"+id]) for id in ids])
        
        #return load(open(self.local_database_path+"list_pathways_"+org+".pickle"))

    def get_linked_pathways(self, pathway_id):
        return ["path:"+p.strip() for p in self._retrieve(self._rel_pathway_dir(pathway_id)+pathway_id.split(":")[-1]+".map", self._pathway_from(pathway_id)).readlines()]

    def get_genes_by_organism(self, org):
        return self._genes[org].keys()
        #return [org+":"+g for g in _collect(self._retrieve("pathway/organisms/"+org+"/"+org+"_gene_map.tab").readlines(), lambda s:s.split()[:1])]

    def get_genes_by_pathway(self, pathway_id):
        self._genes[pathway_id.split(":")[-1][:-5]]
        return self._from_pathway_to_genes[pathway_id]
        ##        return [pathway_id.split(":")[-1][:-5]+":"+g for g in _collect(self._retrieve(self._rel_pathway_dir(pathway_id)+pathway_id.split(":")[-1]+".gene").readlines(), lambda s:s.split()[:1])]

    def get_enzymes_by_pathway(self, pathway_id):
        if pathway_id.startswith("path:map"):
            return self._retrieve(self._rel_pathway_dir(pathway_id)+pathway_id.split(":")[-1]+".enz", "kegg_reference.tar.gz").readlines()
        else:
            genes = self.get_genes_by_pathway(pathway_id)
            return list(set(_collect(map(self.get_enzymes_by_gene, genes))))

    def get_compounds_by_pathway(self, pathway_id):
        return _collect(self._retrieve(self._rel_pathway_dir(pathway_id)+pathway_id.split(":")[-1]+".cpd", self._pathway_from(pathway_id)).readlines(), lambda s:s.split()[:1])

    def get_pathways_by_genes(self, genes_list):
        genes = set(genes_list)
        orgs = set([g.split(":")[0] for g in genes_list])
        if len(orgs)!=1:
            return []
        org = orgs.pop()
        s = set()
        for gene in genes:
            pathways = self._genes[org][gene].get_pathways()
            for path in pathways:
                if genes.issubset(self.get_genes_by_pathway(path)):
                    s.add(path)
        return s
        """d = dict(_collect(self._retrieve("pathway/organisms/"+org+"/"+org+"_gene_map.tab").readlines(), lambda line:(lambda li:(org+":"+li[0], li[1:]))(line.split())))
        s = set(_collect(genes, lambda gene:d.get(gene, [])))
        return list(s)"""
        """pathways = self.list_pathways(orgs.pop())
        return filter(lambda p:genes.issubset(self.get_genes_by_pathway(p)), pathways)"""

    def get_pathways_by_enzymes(self, enzyme_list):
        pathways = enzyme_list and set(self._enzymes.get(enzyme_list[0], DBEnzymeEntry(" ")).get_pathways()) or []
        for enzyme in enzyme_list[1:]:
            pathways&=set(self._enzymes.get(enzyme, DBEnzymeEntry(" ")).get_pathways())
        return list(pathways)

    def get_pathways_by_compounds(self, compound_list):
        pathways = compound_list and set(self._compounds.get(compound_list[0], DBCompoundEntry(" ")).get_pathways()) or []
        for compound in compound_list[1:]:
            pathways&=set(self._compounds.get(compound,DBCompoundEntry(" ")).get_pathways())
        return list(pathways)

    def get_enzymes_by_compound(self, compound_id):
        if compound_id in self._compounds:
            return self._compounds[compound_id].get_enzymes()
        else:
            return []
    
    def get_compounds_by_enzyme(self, enzyme_id):
        return self._from_enzyme_to_compounds.get(enzyme_id, [])
    
    def get_genes_by_enzyme(self, enzyme_id, org=None):
        if enzyme_id in self._enzymes:
            genes = self._enzymes[enzyme_id].get_genes()
            if org:
                return filter(lambda g:g.startswith(org), genes)
            else:
                return genes
        else:
            return []
    
    def get_enzymes_by_gene(self, gene_id):
        return self._from_gene_to_enzymes.get(gene_id, [])

    def get_pathway_image(self, pathway_id):
        f = self._retrieve(self._rel_pathway_dir(pathway_id)+pathway_id.split(":")[-1]+".gif", self._pathway_from(pathway_id), mode="rb")
##        image = Image.open(self.local_database_path+_rel_dir(pathway_id)+pathway_id.split(":")[-1]+".gif")
        image = Image.open(f)
        return image.convert("RGB")

    def get_colored_pathway_image(self, pathway_id, objects):
        color = (255, 0, 0)
        image = self.get_pathway_image(pathway_id)
        #image = image.convert("RGB")
        tmp = Image.new("RGB", image.size)
        draw = ImageDraw.Draw(tmp)
        bb = self.get_bounding_box_dict(pathway_id)
        for object_id in objects:
            t = bb.get(object_id, [])
            for x1, y1, x2, y2 in t:
                draw.rectangle([x1, y1, x2, y2], outline=color)
        del draw
        i1, i2, i3 = image.split()
        t1, t2, t3 = tmp.split()
        i1 = ImageMath.eval("a+b", a=i1, b=t1)
        i2 = ImageMath.eval("a+b", a=i2, b=t2)
        i3 = ImageMath.eval("a+b", a=i3, b=t3)
        return Image.merge("RGB", (i1.convert("L"), i2.convert("L"), i3.convert("L")))

    def get_bounding_box_dict(self, pathway_id):
        org = pathway_id.split(":")[-1][:-5]
        d=[]
        if not pathway_id.split(":")[-1].startswith("map"):
            try:
                d = map(lambda line:(org+":"+line.split()[0], tuple(line.split()[1:])), self._retrieve(self._rel_pathway_dir(pathway_id)+pathway_id.split(":")[-1]+"_gene.coord", self._pathway_from(pathway_id)).readlines())
            except:
                pass
        try:
            d.extend(map(lambda line:("cpd:"+line.split()[0], tuple(line.split()[1:])), self._retrieve(self._rel_pathway_dir(pathway_id)+pathway_id.split(":")[-1]+"_cpd.coord", self._pathway_from(pathway_id)).readlines()))
        except:
            pass
        try:
            for line in self._retrieve(self._rel_pathway_dir(pathway_id)+pathway_id.split(":")[-1]+".conf", self._pathway_from(pathway_id)).readlines():
                match = re.findall("rect \(([0-9]+),([0-9]+)\) \(([0-9]+),([0-9]+)\)	/kegg/pathway/"+org+"/([a-z0-9]+)\.html", line)
                for t in match:
                    d.append(("path:"+t[-1], t[:-1]))
        except:
            pass
        d = [(id, tuple(map(int, t))) for id, t in d]
        bbDict = defaultdict(list)
        for id, bb in d:
            bbDict[id].append(bb)
        return bbDict
        
    def get_unique_gene_ids(self, org, genes, caseSensitive=True):
        if not caseSensitive:
            return self.get_unique_gene_ids_ci(org, genes)
        
        allGenes = self._genes[org]
        unique = {} #[]
        conflicting = []
        unknown = []
        for gene in genes:
            if gene in allGenes:
                unique[gene] = gene #unique.append(gene)
            elif gene in self._gene_alias_conflicting[org]:
                conflicting.append(gene)
            elif gene in self._gene_alias[org]:
                unique[self._gene_alias[org][gene]] = gene #unique.append(self._gene_alias[org][gene])
            else:
                unknown.append(gene)
        return unique, conflicting, unknown


    def constructAliasMapper(self, org):
        try:
            len(self._cigenes)
        except:
            self._cigenes = {}


        try:
            len(self._cigenes[org])
        except:
            allGenes = dict(self._genes[org]) #just to load the database

            for id in self._genes[org].keys():
                if allGenes.get(id.upper(), id)!=id:
                    conf.add(id.upper())
                else:
                    allGenes[id.upper()] = id

            conf = set(self._gene_alias_conflicting[org])

            aliasMapper = {}

            for alias, id in self._gene_alias[org].items():
                if aliasMapper.get(alias.upper(), id)!=id:
                    conf.add(alias.upper())
                else:
                    aliasMapper[alias.upper()]=id

            self._cigenes[org] = (allGenes,aliasMapper,conf)

    def get_unique_gene_ids_ci(self, org, genes):

        self.constructAliasMapper(org)

        unique = {}
        conflicting = []
        unknown = []
    
        allGenes,aliasMapper,conf = self._cigenes[org]

        for gene in genes:
            if gene.upper() in conf:
                conflicting.append(gene)
            elif gene.upper() in allGenes:
                unique[allGenes[gene.upper()]] = gene
            elif gene.upper() in aliasMapper:
                unique[aliasMapper[gene.upper()]] = gene
            else:
                unknown.append(gene)
        return unique, conflicting, unknown

    def get_ko_orthology(self):
        r = []
        f = self._retrieve("brite/ko/ko00001.keg", "kegg_orthology.tar.gz")
        for l in f.readlines():
            if not l.strip("ABCD\n"):
                continue
            if l.startswith("A"):
                r.append(KOClass(l))
            elif l.startswith("B"):
                r[-1].children.append(KOClass(l))
            elif l.startswith("C"):
                r[-1].children[-1].children.append(KOClass(l))
        return r

    def get_gene_name(self, org, gene):
        names = self._genes[org][gene].get_alt_names()
        return names[len(self._genes[org][gene].get_db_links()):-1][0]
    
class KEGGOrganism(object):
    def __init__(self, org, update=False, local_database_path=None):
        
        self.org = org
        self.local_database_path = local_database_path or default_database_path
        if self.local_database_path in KEGGInterfaceLocal._instanceCache:
            self.api = KEGGInterfaceLocal._instanceCache[self.local_database_path]
        else:
            self.api = KEGGInterfaceLocal(update, self.local_database_path)
        if org not in self.api._taxonomy:
            import obiTaxonomy as tax
            ids = tax.to_taxid(org, mapTo=[entry.get_taxid() for entry in self.api._genome.values()])
            ids = set(ids).intersection([entry.get_taxid() for entry in self.api._genome.values()])
##            names = [key for key, name in self.api._taxonomy.items() if org.lower() in name[-1].lower()]
            if not ids:
                print >> sys.stderr, "Could not find", org, "in KEGG database\nSearching in NCBI taxonomy"
                import obiTaxonomy as tax
                ids = tax.search(org)
                ids = set(ids).intersection([entry.get_taxid() for entry in self.api._genome.values()])
            if len(ids) == 0:
                raise tax.UnknownSpeciesIdentifier, org
            elif len(ids) > 1:
                raise tax.MultipleSpeciesException, ", ".join(["%s: %s" % (from_taxid(id), tax.name(id)) for id in ids])
            print >> sys.stderr, "Found", tax.name(list(ids)[0]), "(%s)" % from_taxid(list(ids)[0]) 
            self.org = from_taxid(ids.pop())

    def list_pathways(self):
        """Return a list of all organism specific pathways."""
        return self.api.list_pathways(self.org)
    
    def get_linked_pathways(self, pathway_id):
        """Return a list of all organism specific pathways that pathway with pathway_id links to."""
        return self.api.get_linked_pathways(pathway_id)

    def get_genes_by_pathway(self, pathway_id):
        """Return a list of all organism specific genes that are on the pathway with pathway_id."""
        return self.api.get_genes_by_pathway(pathway_id)

    def get_enzymes_by_pathway(self, pathway_id):
        """Return a list of all organism specific enzymes that are on the pathway with pathway_id."""
        return self.api.get_enzymes_by_pathway(pathway_id)

    def get_compounds_by_pathway(self, pathway_id):
        """Return a list of all organism specific compounds that are on the pathway with pathway_id."""
        return self.api.get_enzymes_by_pathway(pathway_id)

    def get_genes(self):
        """Return a list of all organism genes."""
        return self.api.get_genes_by_organism(self.org)

    def get_pathways_by_genes(self, genes):
        """Return a list of all organism specific pathways that contain all the genes."""
        return self.api.get_pathways_by_genes(genes)

    def get_enriched_pathways_by_genes(self, genes, reference=None, prob=obiProb.Binomial(), callback=None):
        """Return a dictionary with enriched pathways ids as keys and (list_of_genes, p_value, num_of_reference_genes) tuples as items."""
        allPathways = defaultdict(lambda :[[], 1.0, []])
        tmp_callback = self.api.download_progress_callback
        if callback:
            self.api.download_progress_callback = None
        if not reference:
            reference = self.get_genes()
        for i, gene in enumerate(genes):
            pathways = self.get_pathways_by_genes([gene])
            for pathway in pathways:
                allPathways[pathway][0].append(gene)
            if callback:
                callback(i*100.0/len(genes))
        if callback and self.api.download_progress_callback:
            tmp_callback = self.api.download_progress_callback
            self.api.download_progress_callback = None
        self.api.download_progress_callback = tmp_callback
        reference = set(reference)
        for p_id, entry in allPathways.items():
            entry[2].extend(reference.intersection(self.api.get_genes_by_pathway(p_id)))
##            entry[1] = _p(float(len(entry[2]))/len(reference), len(entry[0]), len(genes))
            entry[1] = prob.p_value(len(entry[0]), len(reference), len(entry[2]), len(genes))
        return dict([(pid, (genes, p, len(ref))) for pid, (genes, p, ref) in allPathways.items()])

    def get_pathways_by_enzymes(self, enzymes):
        """Return a list of all organism specific pathways that contain all the enzymes."""
        return self.api.get_pathways_by_enzymes(enzymes)

    def get_pathways_by_compounds(self, compounds):
        """Return a list of all organism specific pathways that contain all the compounds."""
        return self.api.get_pathways_by_compounds(compounds)

    def get_enzymes_by_compound(self, compound_id):
        """Return a list of all organism specific enzymes that are involved in a reaction with compound."""
        return self.api.get_enzymes_by_compound(compound_id)

    def get_compounds_by_enzyme(self, enzyme_id):
        """Return a list of all compounds that are involved in a reaction with the enzyme."""
        return self.api.get_compounds_by_enzyme(enzyme_id)

    def get_genes_by_enzyme(self, enzyme_id):
        """Return a list of all genes that are involved with the production of enzyme."""
        return self.api.get_genes_by_enzyme(enzyme_id, self.org)

    def get_enzymes_by_gene(self, gene_id):
        """Return a list of all enzymes that are a product of gene."""
        return self.api.get_enzymes_by_gene(gene_id)

    def get_unique_gene_ids(self, genes, caseSensitive=True):
        """Return a tuple with three elements. The first is a dictionary mapping from unique gene
        ids to gene names in genes, the second is a list of conflicting gene names and the third is a list
        of unknown genes.
        """
        return self.api.get_unique_gene_ids(self.org, genes, caseSensitive)

    def get_gene_name(self, geneId):
        """ Return all gene names for the given gene id
        """
        return self.api.get_gene_name(self.org, geneId)

class KEGGPathway(object):
    def __init__(self, pathway_id, update=False, local_database_path=None):
        self.pathway_id = pathway_id
        self.org = pathway_id.split(":")[-1][:-5]
        self.local_database_path = local_database_path or default_database_path
        if self.local_database_path in KEGGInterfaceLocal._instanceCache:
            self.api = KEGGInterfaceLocal._instanceCache[self.local_database_path]
        else:
            self.api = KEGGInterfaceLocal(update, self.local_database_path)
        if update:
            self.api.download_pathway_data(self.org)

    def get_image(self):
        """Return an PIL image of the pathway."""
        return self.api.get_pathway_image(self.pathway_id)

    def get_colored_image(self, objects):
        """Return an PIL image of the pathway with marked objects."""
        return self.api.get_colored_pathway_image(self.pathway_id, objects)

    def get_bounding_box(self, object_id):
        """Return a bounding box of the form (x1, y1, x2, y2) of object on the pathway image."""
        return self.api.get_bounding_box(self.pathway_id, object_id)

    def get_bounding_box_dict(self):
        """Return a dictionary mapping all objects on the pathways to bounding boxes (x1, y1, x2, y2) on the pathway image."""
        return self.api.get_bounding_box_dict(self.pathway_id)

    def get_genes(self):
        """Return all genes on the pathway."""
        return self.api.get_genes_by_pathway(self.pathway_id)

    def get_enzymes(self):
        """Return all enzymes on the pathway."""
        return self.api.get_enzymes_by_pathway(self.pathway_id)

    def get_compounds(self):
        """Return all compounds on the pathway."""
        return self.api.get_compounds_by_pathway(self.pathway_id)

def from_taxid(taxid):
    api = KEGGInterfaceLocal()
    org = [key for key, entry in api._genome.items() if entry.get_taxid() == taxid]
    if not org:
        raise ValueError, taxid
    else:
        return org[0]

def to_taxid(org):
    api = KEGGInterfaceLocal()
    return api._genome[org].get_taxid()

class KOClass(object):
    def __init__(self, text=None):
        self.children = []
        self.ko_class_id = "?"
        self.class_name = "?"
        if text:
            self._parse_line(text)
            
    def _parse_line(self, text):
        if text.startswith("A"):
            self.class_name = text.strip("<>AB/ \n")
        elif text.startswith("B"):
            self.class_name = text.strip("<>B/ \n")
        elif text.startswith("C"):
            self.class_name = text.strip("C \n")
            try:
                self.class_name = self.class_name[:self.class_name.index("[")]
            except:
                pass
        self.ko_class_id = self.class_name[:5]

from obiGenomicsUpdate import Update as UpdateBase

import tarfile

class Update(UpdateBase):
    def __init__(self, local_database_path=None, progressCallback=None):
        UpdateBase.__init__(self, local_database_path if local_database_path else default_database_path, progressCallback)
        self.api = KEGGInterfaceLocal("Force update", self.local_database_path, progressCallback)

    def LastUpdate(self, func, args):
        def _LastUpdate(path):
            size, time = self.api.downloader.ftpWorker.statFtp("/pub/kegg/"+path)
            return time
        if func == Update.UpdateOrganism:
            rel_path = self.api._rel_org_dir(args[0]).rstrip("/")
            return max(_LastUpdate("pathway/"+rel_path), _LastUpdate("genes/"+rel_path))
        elif func == Update.UpdateReference:
            return max(_LastUpdate("pathway/map"), _LastUpdate("pathway/map_title.tab"))
        elif func == Update.UpdateEnzymeAndCompounds:
            return max(_LastUpdate("ligand/compound/compound"), _LastUpdate("ligand/enzyme/enzyme"))
        elif func == Update.UpdateOrthology:
            return _LastUpdate("brite/ko/ko00001.keg")
        elif func == Update.UpdateTaxonomy:
            return max(_LastUpdate("genes/taxonomy"), _LastUpdate("genes/genome"))
        
    def IsUpdatable(self, func, args):
        return self.LastUpdate(func, args) > self.GetLastUpdateTime(func, args)
    
    def GetDownloadable(self):
        ret = []
        ret.extend([(Update.UpdateTaxonomy, ())] if (Update.UpdateTaxonomy, ()) not in self.shelve else [])
        ret.extend([(Update.UpdateOrthology, ())] if (Update.UpdateOrthology, ()) not in self.shelve else [])
        ret.extend([(Update.UpdateReference, ())] if (Update.UpdateReference, ()) not in self.shelve else [])
        ret.extend([(Update.UpdateEnzymeAndCompounds, ())] if (Update.UpdateEnzymeAndCompounds, ()) not in self.shelve else [])
        orgs = [org for org in self.api.list_organisms() if (Update.UpdateOrganism, (org,)) not in self.shelve]
        ret.extend([(Update.UpdateOrganism , (org,)) for org in orgs])
        return ret

##    @synchronized(updateLock)
    def UpdateOrganism(self, org):
        self.api.download_organism_data(org)
        rel_path = self.api._rel_org_dir(org)
        try:
            os.remove(os.path.join(self.local_database_path, "genes/", rel_path, "_genes.pickle"))
        except Exception:
            pass
        self.api._load_gene_database(org, from_=".//") #to parse the .ent file and create the _genes.pickle file
        try:
            os.remove(os.path.join(self.local_database_path, "genes/", rel_path, self.api._taxonomy[org][0]+".ent"))
        except Exception:
            pass
        self._update(Update.UpdateOrganism, (org,))

##    @synchronized(updateLock)
    def UpdateReference(self):
        self.api.download_reference_data()
        self._update(Update.UpdateReference, ())

##    @synchronized(updateLock)
    def UpdateEnzymeAndCompounds(self):
        self.api.downloader.massRetrieve(["ligand//compound//compound", "ligand//enzyme//enzyme"], progressCallback=self.progressCallback)
        for file in ["ligand//compound//_compounds.pickle", "ligand//enzyme//_enzymes.pickle", "ligand//enzyme//_from_gene_to_enzymes.pickle", "ligand//compound//_from_enzyme_to_compounds.pickle"]:
            try:
                os.remove(os.path.join(self.local_database_path, file))
            except Exception:
                pass
        self.api._load_compound_database(from_=".//")
        self.api._load_enzyme_database(from_=".//")
        for file in ["ligand//compound//compound", "ligand//enzyme//enzyme"]:
            try:
                os.remove(os.path.join(self.local_database_path, file))
            except Exception:
                pass
        self._update(Update.UpdateEnzymeAndCompounds, ())

    def UpdateTaxonomy(self):
        self.api.downloader.massRetrieve(["genes//taxonomy", "genes//genome"], progressCallback=self.progressCallback)
        self._update(Update.UpdateTaxonomy, ())

    def UpdateOrthology(self):
        self.api.downloader.retrieve("brite//ko//ko00001.keg","brite//ko//ko00001.keg", progressCallback=self.progressCallback)
        self._update(Update.UpdateOrthology, ())

    def GetTarballDirs(self):
        orgs = self.api.list_organisms()
        return ["pathway//organisms//"+org for org in orgs] + ["pathway//map"]


import cStringIO

class Orthology(object):
    pass

class Brite(object):
    pass

class Organism(object):
    def __init__(self, file):
        self.file = file
        self._cachedFiles = {}
        self._cacheNames = {"genes": "genes.pickle"}

    @classmethod
    def Load(cls, org):
        return Organism(os.path.join(default_database_path, org + "organism.tar.gz"))

    def _open(self, filename):
        filename = os.path.normpath(filename)
        if filename not in self._cachedFiles:
            self._cachedFiles[filename] = self._tarfile.extractfile(filename).read().replace("\r\n", "\n")
        return cStringIO.StringIO(self._cachedFiles[filename])
    
    def __getattr__(self, name):
        if name in self._cacheNames:
            setattr(self, name, load(self._open(self._cacheNames[name])))
            return getattr(self, name)
        else:
            raise AttributeError(name)

    def GetPathways(self):
        """Return a list of all organism specific pathways."""
        return ["path:" + info.name.split(".")[0] for info in self._tarfile.getmembers() if info.name.endswith(".xml")]
    
    def GetLinkedPathways(self, pathway_id):
        """Return a list of all organism specific pathways that pathway with pathway_id links to."""
        pathway = KEGGPathwayMk2(self._open(pathway_id.split(":")[-1] + ".xml"))
        return pathway.GetLinkedPathways()

    def GetGenesByPathway(self, pathway_id):
        """Return a list of all organism specific genes that are on the pathway with pathway_id."""
        pathway = KEGGPathwayMk2(self._open(pathway_id.split(":")[-1] + ".xml"))
        return pathway.GetGenes()

    def GetEnzymesByPathway(self, pathway_id):
        """Return a list of all organism specific enzymes that are on the pathway with pathway_id."""
        pathway = KEGGPathwayMk2(self._open(pathway_id.split(":")[-1] + ".xml"))
        return pathway.GetEnzymes()

    def GetCompoundsByPathway(self, pathway_id):
        """Return a list of all organism specific compounds that are on the pathway with pathway_id."""
        pathway = KEGGPathwayMk2(self._open(pathway_id.split(":")[-1] + ".xml"))
        return pathway.GetCompounds()
    
    def GetGenes(self):
        """Return a list of all organism genes."""
        return self.genes.keys()

    def GetPathwaysByGene(self, gene):
        """Return a list of all organism specific pathways that contain the gene."""
        return self.genes[gene].get_pathways()

    def GetEnrichedPathwaysByGenes(self, genes, reference=None, prob=obiProb.Binomial(), progressCallback=None):
        """Return a dictionary with enriched pathways ids as keys and (list_of_genes, p_value, num_of_reference_genes) tuples as items."""
        allPathways = defaultdict(lambda :[[], 1.0, []])
        if not reference:
            reference = self.GetGenes()
        for i, gene in enumerate(genes):
            pathways = self.GetPathwaysByGene(gene)
            for pathway in pathways:
                allPathways[pathway][0].append(gene)
            if progressCallback:
                progressCallback(i*100.0/len(genes))
        
        reference = set(reference)
        for p_id, entry in allPathways.items():
            entry[2].extend(reference.intersection(self.GetGenesByPathway(p_id)))
            entry[1] = prob.p_value(len(entry[0]), len(reference), len(entry[2]), len(genes))
        return dict([(pid, (genes, p, len(ref))) for pid, (genes, p, ref) in allPathways.items()])

    def GetPathwaysByEnzyme(self, enzyme):
        """Return a list of all organism specific pathways that contain the enzyme."""
        return self.api.get_pathways_by_enzymes(enzymes)

    def GetPathwaysByCompounds(self, compound):
        """Return a list of all organism specific pathways that contain the compound."""
        return self.api.get_pathways_by_compounds(compounds)

    def GetEnzymesByCompound(self, compound_id):
        """Return a list of all organism specific enzymes that are involved in a reaction with compound."""
        return self.api.get_enzymes_by_compound(compound_id)

    def GetCompoundsByEnzyme(self, enzyme_id):
        """Return a list of all compounds that are involved in a reaction with the enzyme."""
        return self.api.get_compounds_by_enzyme(enzyme_id)

    def GetGenesByEnzyme(self, enzyme_id):
        """Return a list of all genes that are involved with the production of enzyme."""
        return self.api.get_genes_by_enzyme(enzyme_id, self.org)

    def GetEnzymesByGene(self, gene_id):
        """Return a list of all enzymes that are a product of gene."""
        return self.api.get_enzymes_by_gene(gene_id)

    def GetUniqueGeneIds(self, genes, caseSensitive=True):
        """Return a tuple with three elements. The first is a dictionary mapping from unique gene
        ids to gene names in genes, the second is a list of conflicting gene names and the third is a list
        of unknown genes.
        """
        return self.api.get_unique_gene_ids(self.org, genes, caseSensitive)

    @staticmethod
    def UpdateOrganism(org, local, progressCallback=None):
        buffer = os.path.join(orngEnviron.bufferDir, "kegg_tmp")
        try:
            f = tarfile.open(local)
            f.extractall(buffer)
            f.close()
        except Exception:
            pass

        downloader = obiData.FtpDownloader("ftp.genome.jp", buffer, "/pub/kegg/", numOfThreads=7)
        
        if org == "map":
            rel_path = "pathway/map/"
        else:
            for path in ["pathway/organisms/", "pathway/organisms_kaas/", "pathway/organisms_est/"]:
                downloader.ftpWorker.statFtp("pub/kegg/" + path)
                if org in downloader.ftpWorker.statCache.get("pub/kegg/" + path, {}):
                    rel_path = path + org + "/"
                    
        if org != "map":
            rel_genes_path = rel_path.replace("pathway", "genes")
            downloader.ftpWorker.statFtp("pub/kegg/" + rel_genes_path)
            files = downloader.ftpWorker.statCache["pub/kegg/" + rel_genes_path].keys()
            names = [name for name in files if name.endswith(".ent")]
            if len(names) == 1:
                downloader.massRetrieve([(rel_genes_path + names[0], "genes")], blocking=False)
                
        xml_rel_path= "xml/map/" if map == org else "xml/organisms/" + org + "/"
        downloader.ftpWorker.statFtp("pub/kegg/" + xml_rel_path)
        pathways = [name.split(".")[0] for name in downloader.ftpWorker.statCache.get("pub/kegg/" + xml_rel_path, {}).keys()]
        #print pathways
        
        gif_files = [(rel_path + pathway_id + ".gif", pathway_id + ".gif") for pathway_id in pathways]
        downloader.massRetrieve(gif_files, blocking=True)
        
        xml_files = [(xml_rel_path + pathway_id + ".xml", pathway_id + ".xml") for pathway_id in pathways]
        downloader.massRetrieve(xml_files, blocking=False)
        
        downloader.wait(progressCallback)
        try:
            genes = [DBEntryWrapper(DBGeneEntry(entry)) for entry in open(os.path.join(buffer, "genes")).read().split("///\n") if entry]
            genes = dict([(org + ":" + entry.get_name(), entry) for entry in genes])
            dump(genes, open(os.path.join(buffer, "genes.pickle"), "wb"))
        except Exception:
            pass
        os.remove(os.path.join(buffer, "genes"))
        f = tarfile.open(local, "w")
        for filename in os.listdir(buffer):
            f.add(os.path.join(buffer, filename), filename)
        f.close()
        import shutil
        shutil.rmtree(buffer)
        
from xml.dom.minidom import parse, Text, Element

class Pathway(object):
    def __init__(self, file):
        self.xml = parse(file)
        self.xml_pathway = self.xml.childNodes[-1]
        self.xml_elements = [node for node in self.xml_pathway.childNodes if node.__class__ == Element]
        self.entrys = self.xml_pathway.getElementsByTagName("entry")
        self.reactions = self.xml_pathway.getElementsByTagName("reaction")
        self.relations = self.xml_pathway.getElementsByTagName("relation")
        self.genes = [entry for entry in self.entrys if entry.attributes["type"].value == "gene"]
        self.enzymes = [entry for entry in self.entrys if entry.attributes["type"].value == "enzyme"]
        self.compounds = [entry for entry in self.entrys if entry.attributes["type"].value == "compound"]

    @classmethod
    def Load(cls, pathway_id):
        o = Organism.Load(pathway_id.split(":")[-1][:-5])
        return Pathway(o._open(pathway_id+".xml"))

    def GetImage(self):
        """Return an PIL image of the pathway."""
        return Image.open(self._tarfile.extractfile(self.pathway_id.split(":")[-1] + ".gif"))

    def GetColoredImage(self, objects):
        """Return an PIL image of the pathway with marked objects."""
        raise NotImplementedError("GetColoredImage")

    def GetBoundingBox(self, object_id):
        """Return a bounding box of the form (x1, y1, x2, y2) of object on the pathway image."""
        return self.GetBoundingBoxDict()[object_id]

    def GetBoundingBoxDict(self):
        """Return a dictionary mapping all objects on the pathways to bounding boxes (x1, y1, x2, y2) on the pathway image."""
        result = {}
        for entry in self.entrys:
            for graphics in entry.getElementsByTagName("graphics"):
                result[entry.attributes["name"].value] = (graphics.attributes["type"].value,
                                                          graphics.attributes["x"].value, graphics.attributes["y"].value,
                                                          graphics.attributes["width"].value, graphics.attributes["height"].value)
        return result

    def GetLinkedPathways(self):
        return [entry.attributes["name"].value for entry in self.entrys if entry.attributes["type"].value == "map"]

    def GetGenes(self):
        """Return all genes on the pathway."""
        return [gene.attributes["name"].value for gene in self.genes]

    def GetEnzymes(self):
        """Return all enzymes on the pathway."""
        return [enzyme.attributes["name"].value for enzyme in self.enzymes]

    def GetCompounds(self):
        """Return all compounds on the pathway."""
        return [compound.attributes["name"].value for compound in self.compounds]

class UpdateMk2(UpdateBase):
    def __init__(self, local_database_path=None, progressCallback=None):
        UpdateBase.__init__(self, local_database_path if local_database_path else default_database_path, progressCallback)
        self.downloader = obiData.FtpDownloader("ftp.genome.jp", self.local_database_path, "/pub/kegg/", numOfThreads=10)
        
    def IsUpdatable(self, *args):
        return True

    def GetDownloadable(self, *args):
        pass
        
    def UpdateOrganism(self, org):
        KEGGOrganismMk2.UpdateOrganism(org, os.join(self.local_database_path, org + "_organism.tar.gz"))
    
            

if __name__=="__main__":
    
    org1 = KEGGOrganism("ddi")
    org2 = KEGGOrganism("ddi")
    org2.api = KEGGInterface()
    tests = [("get_genes", ()),
             ("get_genes_by_enzyme", ("ec:1.1.1.1",)),
             ("get_genes_by_pathway", ("path:ddi00010",)),
             ("get_pathways_by_genes", (["ddi:DDB_0191256"],)),
             ("get_pathways_by_enzymes", (["ec:1.1.1.1"],)),
             ("get_pathways_by_compounds", (["cpd:C00001"],)),
             ("get_linked_pathways", ("path:ddi00010",)),
             ("list_pathways", ()),
             ("get_compounds_by_enzyme", ("ec:1.1.1.1",)),
             ("get_compounds_by_pathway", ("path:ddi00010",)),
             ("get_enzymes_by_compound", ("cpd:C00001",)),
             ("get_enzymes_by_pathway", ("path:ddi00010",)),
             ("get_enzymes_by_gene", ("ddi:DDB_0191256",))]
    for name, args in tests:
        s1 = set(getattr(org1, name)(*args))
        s2 = set(getattr(org2, name)(*args))
        if s1 and s2:
            print name
            print s1-s2
            print s2-s1
        else:
            print name
            print "both empty"
