"""obiKEGG is an interface to Kyoto Encyclopedia of Genes and Genomes (http://www.genome.jp/kegg/) that allows easy access to KEGG pathway and genes data.

"""
import Image, ImageDraw, ImageMath
import cStringIO
import math
import time
import os
import re

import obiData

from cPickle import load, dump
from collections import defaultdict

try:
    import orngOrangeFoldersQt4
    default_database_path = os.path.join(orngOrangeFoldersQt4.__getDirectoryNames()['bufferDir'],"kegg//")
except:
    default_database_path = os.path.join((os.path.split(__file__)[0] or "."), "data//kegg//")

base_ftp_path = "ftp://ftp.genome.jp/pub/kegg/"

import htmllib
class HTMLImageCollector(htmllib.HTMLParser):
    def __init__(self):
        self.images = []
    def handle_image(self, source, *args):
        self.images.append(source)

def _image_from_file(f):
    imgstr = f.read()
    return Image.open(cStringIO.StringIO(imgstr))
    
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
##        if org:
##            return [org.lower()+":"+g.split("(")[0] for g in d.get(org.upper(), [])]
##        else:
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
    def __init__(self, update=False, local_database_path=None, download_progress_callback=None):
        self.local_database_path = local_database_path or default_database_path
        self.update = update
        self.download_progress_callback = download_progress_callback
        self._gene_alias = {}
        self._gene_alias_conflicting = {}
        self._filenames = {"_enzymes":"ligand/enzyme/_enzymes.pickle",
                           "_from_gene_to_enzymes":"ligand/enzyme/_from_gene_to_enzymes.pickle",
                           "_compounds":"ligand/compound/_compounds.pickle",
                           "_from_enzyme_to_compounds":"ligand/compound/_from_enzyme_to_compounds.pickle"}
        self.downloader = obiData.FtpDownloader("ftp.genome.jp", self.local_database_path, "/pub/kegg/", numOfThreads=10)

    def download_organism_data(self, org):
        rel_path = "pathway/organisms/"+org+"/"
        files = [rel_path+org+"_gene_map.tab", "pathway/map_title.tab", "genes/taxonomy"]
        self.downloader.massRetrieve(files, update=self.update)
        file = self._retrieve(rel_path+org+"_gene_map.tab")
        pathway_nums = set(reduce(lambda a,b: a + b.split()[1:], file.readlines(), []))
        descr = dict(map(lambda line:tuple(line.strip().split("\t")), self._retrieve("pathway/map_title.tab").readlines()))
        dump(descr, open(self.local_database_path+"list_pathways_map.pickle", "w"))
        ids = [org+num for num in pathway_nums]
        try:
            organisms = load(open(self.local_database_path+"list_organisms.pickle"))
        except:
            organisms = {}
        if org not in organisms:
            organisms[org] = self._taxonomy.get(org, "  ")[1]
            dump(organisms, open(self.local_database_path+"list_organisms.pickle", "w"))
        dump(dict([("path:"+org+num, descr[num]) for num in pathway_nums]), open(self.local_database_path+"list_pathways_"+org+".pickle","w"))
        
        ends = [".cpd", ".gene", ".gif", ".map", "_cpd.coord", "_gene.coord", ".conf"]
        files = [rel_path+id+ext for id in ids for ext in ends]
        self.downloader.massRetrieve(files, update=self.update, blocking=False)
        self.downloader.retrieve("genes/organisms/"+org+"/"+self._taxonomy[org][0]+".ent", update=self.update, progressCallback=self.download_progress_callback)
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
            return self.__dict__[name]
        elif name=="_compounds" or name=="_from_enzyme_to_compounds":
            self._load_compound_database()
            return self.__dict__[name]
        elif name=="_genes":
            self._genes = GenesDatabaseProxy(self)
            return self.__dict__[name]
            #self._load_genes_database()
        elif name=="_taxonomy":
            self._load_taxonomy()
            return self.__dict__[name]
        else:
            raise AttributeError(name)

    def _load_pickled(self, filename=None, name=None):
        if not filename and name:
            return load(open(self.local_database_path+self._filenames[name]))
        else:
            return load(open(self.local_database_path+filename))

    def _dump_pickled(self, object, filename=None, name=None):
        if not  filename and name:
            dump(object, open(self.local_database_path+self._filenames[name], "w"))
        else:
            dump(object, open(self.local_database_path+filename, "w"))
    
    def _load_enzyme_database(self):
        try:
            self._enzymes = self._load_pickled(name="_enzymes")
        except Exception, ex:
            print ex
            enzymes = map(DBEnzymeEntry, filter(bool, self._retrieve("ligand/enzyme/enzyme").read().split("///\n")))
            self._enzymes = dict([(e.get_name(), DBEntryWrapper(e)) for e in enzymes])
            self._dump_pickled(self._enzymes, name="_enzymes")
        try:
            self._from_gene_to_enzymes = self._load_pickled(name="_from_gene_to_enzymes")
        except Exception, ex:
            self._from_gene_to_enzymes = defaultdict(list)
            for id, e in self._enzymes.items():
                for g in e.get_genes():
                    self._from_gene_to_enzymes[g].append(id)
            self._dump_pickled(self._from_gene_to_enzymes, name="_from_gene_to_enzymes")
        
    def _load_compound_database(self):
        try:
            self._compounds = self._load_pickled(name="_compounds")
        except:
            compounds = map(DBCompoundEntry, filter(bool, self._retrieve("ligand/compound/compound").read().strip().split("///\n")))
            self._compounds = dict([(c.get_name(), DBEntryWrapper(c)) for c in compounds])
            self._dump_pickled(self._compounds, name="_compounds")
        try:
            self._from_enzyme_to_compounds = self._load_pickled(name="_from_enzyme_to_compounds")
        except:
            self._from_enzyme_to_compounds = defaultdict(list)
            for id, c in self._compounds.items():
                for e in c.get_enzymes():
                    self._from_enzyme_to_compounds[e].append(id)
            self._dump_pickled(self._from_enzyme_to_compounds, name="_from_enzyme_to_compounds")

    def _load_gene_database(self, org):
        try:
            self._genes[org] = self._load_pickled("genes/organisms/"+org+"/_genes.pickle")
        except Exception, ex:
            genes = map(DBGeneEntry, filter(bool ,self._retrieve("genes/organisms/"+org+"/"+self._taxonomy[org][0]+".ent").read().split("///\n")))
            self._genes[org] = dict([(org+":"+g.get_name(), DBEntryWrapper(g)) for g in genes])
            self._dump_pickled(self._genes[org], "genes/organisms/"+org+"/_genes.pickle")
        self._gene_alias[org] = {}
        self._gene_alias_conflicting[org] = set()
        for id, gene in self._genes[org].items():
            aliases = gene.get_alt_names()
            for alias in set(aliases):
                if alias in self._gene_alias[org]:
                    self._gene_alias_conflicting[org].add(alias)
                else:
                    self._gene_alias[org][alias] = id
        dump(set(self._gene_alias[org].keys() + self._genes[org].keys()), open(self.local_database_path+org+"_genenames.pickle","w"))
        return self._genes[org]

    def _load_taxonomy(self):
        orgs = filter(lambda line:line.strip() and not line.startswith("#"), self._retrieve("genes/taxonomy").readlines())
        d = dict([(line.split()[1].strip(), (line.split("\t")[-2].strip(), line.split("\t")[-1].strip())) for line in orgs])
        self._taxonomy = d
        
    def _retrieve(self, filename):
        self.downloader.retrieve(filename, update=self.update, progressCallback=self.download_progress_callback)
        return open(self.local_database_path+filename)
    
    def list_organisms(self):
        return dict([(key, value[1]) for key, value in self._taxonomy.items()])
        #return load(open(self.local_database_path+"list_organisms.pickle"))
    
    def list_pathways(self, org="map"):
        if org=="map":
            r = map(lambda line:tuple(line.strip().split("\t")), self._retrieve("pathway/map_title.tab").readlines())
            return dict([("path:map"+p, desc) for p, desc in r])
        else:
            ids = set(_collect(self._retrieve("pathway/organisms/"+org+"/"+org+"_gene_map.tab").readlines(), lambda line:line.split()[1:]))
            pathways = self.list_pathways("map")
            return dict([("path:"+org+id, pathways["path:map"+id]) for id in ids])
        
        #return load(open(self.local_database_path+"list_pathways_"+org+".pickle"))

    def get_linked_pathways(self, pathway_id):
        return ["path:"+p.strip() for p in self._retrieve(_rel_dir(pathway_id)+pathway_id.split(":")[-1]+".map").readlines()]

    def get_genes_by_organism(self, org):
        return self._genes[org].keys()
        #return [org+":"+g for g in _collect(self._retrieve("pathway/organisms/"+org+"/"+org+"_gene_map.tab").readlines(), lambda s:s.split()[:1])]

    def get_genes_by_pathway(self, pathway_id):
        return [pathway_id.split(":")[-1][:-5]+":"+g for g in _collect(self._retrieve(_rel_dir(pathway_id)+pathway_id.split(":")[-1]+".gene").readlines(), lambda s:s.split()[:1])]

    def get_enzymes_by_pathway(self, pathway_id):
        if pathway_id.startswith("path:map"):
            return self._retrieve(_rel_dir(pathway_id)+pathway_id.split(":")[-1]+".enz").readlines()
        else:
            genes = self.get_genes_by_pathway(pathway_id)
            return list(set(_collect(map(self.get_enzymes_by_gene, genes))))

    def get_compounds_by_pathway(self, pathway_id):
        return _collect(self._retrieve(_rel_dir(pathway_id)+pathway_id.split(":")[-1]+".cpd").readlines(), lambda s:s.split()[:1])

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
        f = self._retrieve(_rel_dir(pathway_id)+pathway_id.split(":")[-1]+".gif")
        image = Image.open(self.local_database_path+_rel_dir(pathway_id)+pathway_id.split(":")[-1]+".gif")
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
                d = map(lambda line:(org+":"+line.split()[0], tuple(line.split()[1:])), self._retrieve(_rel_dir(pathway_id)+pathway_id.split(":")[-1]+"_gene.coord").readlines())
            except:
                pass
        try:
            d.extend(map(lambda line:("cpd:"+line.split()[0], tuple(line.split()[1:])), self._retrieve(_rel_dir(pathway_id)+pathway_id.split(":")[-1]+"_cpd.coord").readlines()))
        except:
            pass
        try:
            for line in self._retrieve(_rel_dir(pathway_id)+pathway_id.split(":")[-1]+".conf").readlines():
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
            allGenes = self._genes[org] #just to load the database

            for id, entry in self._genes[org].items():
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
        f = self._retrieve("brite/ko/ko00001.keg")
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
        
class p_value(object):
    def __init__(self, max=1000):
        self.max = max
        self.lookup = [0]*(max+1)
        for i in xrange(2, max+1):
            self.lookup[i] = self.lookup[i-1] + math.log(i)
            
    def logbin(self, n ,r):
        return self.lookup[n] - self.lookup[n-r] - self.lookup[r]
    
    def binomial(self, n, r, p):
        if p==0.0:
            if r==0:
                return 0.0
            else:
                return 1.0
        elif p==1.0:
            if n==r:
                return 0.0
            else:
                return 1.0
        return math.exp(self.logbin(n, r) + r*math.log(p) + (n + r)*math.log(1.0-p))
    
    def __call__(self, p, mapped, all):
        return reduce(lambda sum, i: sum+self.binomial(all, i, p), range(mapped, all+1), 0.0)
    
class KEGGOrganism(object):
    def __init__(self, org, update=False, local_database_path=None):
        
        self.org = org
        self.local_database_path = local_database_path or default_database_path
        self.api = KEGGInterfaceLocal(update, self.local_database_path)

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

    def get_enriched_pathways_by_genes(self, genes, reference=None, callback=None):
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
        _p = p_value(len(genes))
        reference = set(reference)
        for p_id, entry in allPathways.items():
            entry[2].extend(reference.intersection(self.get_genes_by_pathway(p_id)))
            entry[1] = _p(float(len(entry[2]))/len(reference), len(entry[0]), len(genes))
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

class KEGGPathway(object):
    def __init__(self, pathway_id, update=False, local_database_path=None):
        self.pathway_id = pathway_id
        self.org = pathway_id.split(":")[-1][:-5]
        self.local_database_path = local_database_path or default_database_path
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
from obiGenomicsUpdate import synchronized

from threading import Lock
updateLock = Lock()

class Update(UpdateBase):
    def __init__(self, local_database_path=None, progressCallback=None):
        UpdateBase.__init__(self, local_database_path if local_database_path else default_database_path, progressCallback)
        self.api = KEGGInterfaceLocal(False, self.local_database_path, progressCallback)

    @synchronized(updateLock)
    def GetUpdatable(self):
        orgs = [org for org in self.api.list_organisms() if str((Update.UpdateOrganism, (org,))) in self.shelve]
        ret = [(Update.UpdateOrganism, "Update organism pathways and genes" , orgs)] if orgs else []
        ret.extend([(Update.UpdateReference, "Update reference pathways", [])] if str((Update.UpdateReference, ())) in self.shelve else [])
        ret.extend([(Update.UpdateEnzymeAndCompounds, "Update enzyme and compounds", [])] if str((Update.UpdateEnzymeAndCompounds)) in self.shelve else[])
        return ret
        
    @synchronized(updateLock)
    def GetDownloadable(self):
        orgs = [org for org in self.api.list_organisms() if str((Update.UpdateOrganism, (org,))) not in self.shelve]
        ret = [(Update.UpdateOrganism, "Update organism pathways and genes" , orgs)] if orgs else []
        ret.extend([(Update.UpdateReference, "Update reference pathways", [])] if str((Update.UpdateReference, ())) not in self.shelve else [])
        ret.extend([(Update.UpdateEnzymeAndCompounds, "Update enzyme and compounds", [])] if str((Update.UpdateEnzymeAndCompounds, ())) not in self.shelve else [])
        return ret

    @synchronized(updateLock)
    def UpdateOrganism(self, org):
        self.api.download_organism_data(org)
        self._update(Update.UpdateOrganism, (org,))

    @synchronized(updateLock)
    def UpdateReference(self):
        self.api.download_reference_data()
        self._update(Update.UpdateReference, ())

    @synchronized(updateLock)
    def UpdateEnzymeAndCompounds(self):
        self.api.massDownloader.retrieve(["ligand//compound//compound", "ligand/enzyme/enzyme"], progressCallback=self.progressCallback)
        self._update(Update.UpdateEnzymeAndCompounds, ())

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
