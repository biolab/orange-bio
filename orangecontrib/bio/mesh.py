"""
Module for browsing any analyzing sets annotated with MeSH ontology.
"""

import Orange
from xml.sax import make_parser
from xml.dom import minidom
from xml.sax.handler import ContentHandler
from math import log, exp
from urllib import urlopen
from sgmllib import SGMLParser
import os.path

import orange
from Orange.orng import orngServerFiles

FUZZYMETAID = -12

import Orange.bio.utils.stats
HYPERG = Orange.bio.utils.stats.Hypergeometric()

class obiMeSH(object):

    def __init__(self):
        """
        Function will initialize obiMeSH module and (in case of no MeSH ontology data) download necessary data using orngServerFiles.
        """
        # print "SVN version"
        # we check if all files from MeSH directory are localy present. if not, download them
        self.reference = None
        self.cluster = None
        self.ratio = 1
        self.statistics = None
        self.calculated = False

        self.ref_att = "Unknown"
        self.clu_att = "Unknown"
        self.solo_att = "Unknown"

        # we calculate log(i!)
        self.dataLoaded = self.__loadOntologyFromDisk()

    def expandToFuzzyExamples(self, examples, att, a, b):
        """
        Function will return new 'fuzzy' example table. Every example from the input table will get two additional meta attributes ('fuzzy set' and 'u') \
        based on 'a' and 'b' threshold (lower and higher) and attribute 'att'. Attribute 'fuzzy set' indicates name of the fuzzy set while atribute 'u' \
        reflects example's degree of membership to particular fuzzy set. Note that input examples with values of 'att' lying on the (a,b) will be expanded \
        into two fuzzy examples.
        """
        mu = orange.FloatVariable("u")
        mv = orange.StringVariable("fuzzy set")
        examples.domain.addmeta(FUZZYMETAID, mu)
        examples.domain.addmeta(FUZZYMETAID - 1, mv)
        newexamples = []
        for j in range(0, len(examples)):
            i = examples[j]
            v = float(i[att])
            if v > a and v < b:  # we have to expand this example
                newexamples.append(i)
                i["fuzzy set"] = 'yes'
                i["u"] = (v - a) / (b - a)
                examples.append(i)
                examples[-1]["fuzzy set"] = "no"
                examples[-1]["u"] = (b - v) / (b - a)
            else:
                if v > a:  # u(yes) = 1.0
                    i["fuzzy set"] = 'yes'
                    i["u"] = 1.0
                else:  # u(no) = 1.0
                    i["fuzzy set"] = 'no'
                    i["u"] = 1.0
        return examples

    def findSubset(self, examples, meshTerms, callback=None, MeSHtype='term'):
        """
        Function findSubset will return new example table containing examples which are annotated with at least one MeSH term from list 'meshTerms'.
        """
        # clone
        newdata = orange.ExampleTable(examples.domain)
        self.solo_att = self.__findMeshAttribute(examples)
        ids = list()
        l = len(examples)
        c = 0.0
        meshTerms = list(set(meshTerms))
        # we couldn't find any mesh attribute
        if self.solo_att == "Unknown":
            return newdata
        if MeSHtype == 'term':
            for i in meshTerms:
                ids.extend(self.toID[i])
        else:
            ids = meshTerms
        for e in examples:
            try:
                if callback and c % 10 == 0:
                    callback(int(c * 100.0 / l))
                c = c + 1.0
                ends = eval(e[self.solo_att].value)
            except SyntaxError:
                # print "Error in parsing ", e[self.solo_att].value
                continue
            endids = list()
            for i in ends:
                if self.toID.has_key(i):
                    endids.extend(self.toID[i])
            #allnodes = self.__findParents(endids)
            allnodes = endids
            # calculate intersection
            isOk = False
            for i in allnodes:
                if ids.count(i) > 0:
                    isOk = True
                    break

            if isOk:	  # intersection between example mesh terms and observed term group is None
                newdata.append(e)
        return newdata

    def findTerms(self, ids, idType="cid", callback=None):
        """
        Function findTerms returns a dictionary containing annotated items from the list 'ids'.
        """
        ret = dict()
        if(not self.dataLoaded):
            print "Annotation and ontology has never been loaded!"
            return ret

        if idType == "cid":
            for i in ids:
                if self.fromCID.has_key(int(i)):  # maybe it is on our database
                    ret[int(i)] = self.fromCID[int(i)]
                else:   # no it is not, let's go on the internet - use pubChemAPI
                    r = pubChemAPI()
                    l = r.getMeSHterms(i)
                    ca = r.getPharmActionList(i)
                    l.extend(ca)  # we extend it with Pharmacological actions
                    #l = self.findTermsForCID(i)
                    self.fromCID[int(i)] = l
                    ret[int(i)] = l

                    # if len(l)>0:
                    # if we found something lets save it to a file
                    #   TODO: Never use __file__ to access data, but pkg_resources so that it works also in eggs
                    #	__dataPath = os.path.join(os.path.dirname(__file__), self.path)
                    #	fileHandle = open(os.path.join(__dataPath,'cid-annotation.dat'), 'a')
                    #	for s in l:
                    #		fileHandle.write('\n' + str(i) + ';' + s )
                    #	fileHandle.close()
            return ret
        elif idType == "pmid":  # FIXME PMID annotation
            database = self.fromPMID
        else:
            return ret

        for i in ids:
            if(database.has_key(int(i))):
                ret[i] = database[int(i)]
        return ret

    def findCIDSubset(self, examples, meshTerms):
        """
        Function findCIDSubset will return new list of examples (ids) containing examples which are annotated with at least one MeSH term from \
        list 'meshTerms'.
        """
        newdata = []
        self.solo_att = self.__findMeshAttribute(examples)
        ids = list()
        meshTerms = list(set(meshTerms))
        # we couldn't find any mesh attribute
        if self.solo_att == "Unknown":
            return newdata
        for e in examples:
            try:
                ends = eval(e[self.solo_att].value)
            except SyntaxError:
                continue
            for i in meshTerms:
                if i in ends:
                    newdata.append(e['cid'].value)
                    break
        return newdata

    def findFrequentTerms(self, data, minSizeInTerm, treeData=False, callback=None):
        """
        Function findFrequentTerms iterates thru examples in data. For every example it finds appropriate annotation (MeSH term). Result of this \
        function is a dictionary where keys corespond to MeSH terms ids and a value is a number of examples annotated with particular MeSH term.
        """
        # we build a dictionary 		meshID -> [description, noReference, [cids] ]
        self.statistics = dict()
        self.calculated = False
        self.solo_att = self.__findMeshAttribute(data)
        # post processing variables
        ret = dict()
        ids = []
        succesors = dict()		# for each term id -> list of succesors
        succesors["tops"] = []
        # if we can't identify mesh attribute we return empty data structures
        if self.solo_att == "Unknown":
            if treeData:
                return succesors, ret
            else:
                return ret
        # plain frequency
        t = 0.0
        n = len(data)
        for i in data:
            t = t + 1
            if callback and t % 10 == 0:
                callback(int(t * 100 / n))
            try:
                endNodes = list(set(eval(i[self.solo_att].value)))  # for every CID we look up for end nodes in mesh. for every end node we have to find its ID
            except SyntaxError:
                # print "Error in parsing ",i[self.solo_att].value
                continue
            # we find ID of end nodes
            endIDs = []
            for k in endNodes:
                if(self.toID.has_key(k)):  # this should be always true, but anyway ...
                    endIDs.extend(self.toID[k])
                else:
                    print "Current ontology does not contain MeSH term ", k, "."
            # we find id of all parents
            #allIDs = self.__findParents(endIDs)
            allIDs = endIDs
            for k in allIDs:  # for every meshID we update statistics dictionary
                if(not self.statistics.has_key(k)):  # first time meshID
                    self.statistics[k] = 0
                self.statistics[k] += 1  # counter increased
        # post processing
        for i in self.statistics.iterkeys():
            if(self.statistics[i] >= minSizeInTerm):
                ret[i] = self.statistics[i]
                ids.append(i)
        # we can also return tree data
        if treeData:  # we return also compulsory data for tree printing
            return self.__treeData(ids), ret
        else:
            return ret

    def __treeData(self, ids):
        succesors = dict()
        succesors["tops"] = []
        # we build a list of succesors. for each node we know which are its succesors in mesh ontology
        for i in ids:
            succesors[i] = []
            for j in ids:
                if(i != j and self.__isPrecedesor(i, j)):
                    succesors[i].append(j)
        # for each node from list above we remove its indirect succesors
        # only  i -1-> j   remain
        for i in succesors.iterkeys():
            succs = succesors[i]
            second_level_succs = []
            for k in succs:
                second_level_succs.extend(succesors[k])
            for m in second_level_succs:
                if succesors[i].count(m) > 0:
                    succesors[i].remove(m)
        # we make a list of top nodes
        tops = list(ids)
        for i in ids:
            for j in succesors[i]:
                tops.remove(j)
        # we pack tops table and succesors hash
        succesors["tops"] = tops
        return succesors

    def findEnrichedTerms(self, reference, cluster, pThreshold=0.05, treeData=False, callback=None, fuzzy=False):
        """
        Function findEnrichedTerms computes MeSH term enrichment based on a 'reference' and 'cluster' sets. It returns a dictionary where \
        keys are enriched (their p-value lower that 'pThreshold') MeSH terms. Key values are lists made of several items (MeSH term id, \
        MeSH term description, number of examples from the reference set, number of examples from the cluster set, p value, fold enrichment).
        """
        self.clu_att = self.__findMeshAttribute(cluster)
        self.ref_att = self.__findMeshAttribute(reference)
        if((not self.calculated or self.reference != reference or self.cluster != cluster) and self.ref_att != "Unknown" and \
           self.clu_att != "Unknown"):  # Do have new data? Then we have to recalculate everything.
            self.reference = reference
            self.cluster = cluster
            self.__calculateAll(callback, fuzzy)
        # declarations
        ret = dict()
        ids = []
        succesors = dict()		# for each term id -> list of succesors
        succesors["tops"] = []
        # if some attributes were unknown
        if (self.clu_att == "Unknown" or self.ref_att == "Unknown"):
            if treeData:
                return succesors, ret
            else:
                return ret
        for i in self.statistics.iterkeys():
            if(self.statistics[i][2] <= pThreshold):  # or self.statistics[i][4] <= pThreshold ): #
                ret[i] = self.statistics[i]
                ids.append(i)
        if treeData:
            return self.__treeData(ids), ret
        else:
            return ret

    def printMeSH(self, data, selection=["term", "r", "c", "p"], func=None):
        """
        Function printMeSH can be used to print result (dictionary) from the findEnrichedTerms and findFrequentTerms functions.
        """
        # first we calculate additional info for printing MeSH ontology
        info = self.__treeData(data.keys())
        for i in info["tops"]:
            self.__pp(0, i, info, data, selection, funct=func)

    def __pp(self, offset, item, relations, data, selection, funct=None):
        mapping = {"term": 0, "desc": 1, "r": 2, "c": 3, "p": 4, "fold": 5, "func": 6}
        for i in range(0, offset):
            print " ",
        if type(data[item]) == list:
            pval = "%.4g" % data[item][2]
            fold = "%.4g" % data[item][3]
            print_data = [self.toName[item], self.toDesc[self.toName[item]], str(data[item][0]), str(data[item][1]), str(pval), str(fold)]
            for i in selection:
                if i != "term":
                    print i + "=" + print_data[mapping[i]],
                else:
                    print print_data[mapping[i]],
            if funct != None:
                print " ", funct(print_data[0]),
            # print self.toName[item], " r=" + str(data[item][1])  +" c="+ str(data[item][2])  ," p=" + str(pval) + " fold=" + str(fold)
            print ""
        else:
            print self.toName[item], " freq=" + str(data[item])
        for i in relations[item]:
            self.__pp(offset + 2, i, relations, data, selection, funct=funct)

    def printHtmlMeSH(self, data, selection=["term", "r", "c", "p"], func=None):
        """
        Function printHtmlMeSH if used to print results (dictinary) from the findFrequentTerms and findFrequentTerms functins in HTML format. \
        Together with the MeSH ontology it prints data like number of examples, p-values (enrichment).
        """
        # first we calculate additional info for printing MeSH ontology
        info = self.__treeData(data.keys())
        w = {"term": "'95px'", "r": "'70px'", "c": "'70px'", "p": "'95px'"}
        print "<table>\n<tr>"
        for i in selection:
            print "<th width=" + w[i] + " align='left'>" + i + "</th>"
        if func != None:
            func("header", "")
        print "</tr>\n"
        for i in info["tops"]:
            self.__htmlpp(0, i, info, data, selection, funct=func)
        print "</table>"

    def __htmlpp(self, offset, item, relations, data, selection, funct=None):
        mapping = {"term": 0, "desc": 1, "r": 2, "c": 3, "p": 4, "fold": 5, "func": 6}
        print "<tr>"
        if type(data[item]) == list:
            pval = "%.4g" % data[item][2]
            fold = "%.4g" % data[item][3]
            print_data = [self.toName[item], self.toDesc[self.toName[item]], str(data[item][0]), str(data[item][1]), str(pval), str(fold)]
            for i in selection:
                print "<td>"
                if i == "term":
                    for l in range(0, offset):
                        print "&nbsp;",
                elif i == "p":
                    print '%(#)2.3e' % {'#': float(print_data[mapping[i]])} + "</td>",
                    continue
                print print_data[mapping[i]] + " &nbsp;</td>",

            if funct != None:
                print funct(print_data[0], item),

            # print self.toName[item], " r=" + str(data[item][1])  +" c="+ str(data[item][2])  ," p=" + str(pval) + " fold=" + str(fold)
            print "</tr>"
        else:
            print self.toName[item], " freq=" + str(data[item])

        for i in relations[item]:
            self.__htmlpp(offset + 2, i, relations, data, selection, funct=funct)

    def parsePubMed(self, filename, attributes=["pmid", "title", "abstract", "mesh"], skipExamplesWithout=["mesh"]):
        """
        Function parsePubMed can be used to parse (into Orange example table) PubMed search results (in XML).
        """
        parser = make_parser()
        handler = pubMedHandler()
        parser.setContentHandler(handler)
        parser.parse(open(filename))
        atts = []
        for i in attributes:
            atts.append(orange.StringVariable(i))
        domain = orange.Domain(atts, 0)
        data = orange.ExampleTable(domain)
        print data.domain.attributes
        mapping = {"pmid": 0, "title": 1, "abstract": 2, "mesh": 3, "affilation": 4}
        for i in handler.articles:
            r = []
            skip = False
            for f in attributes:
                if skipExamplesWithout.count(f) > 0:
                    if (f == "mesh" and len(i[mapping[f]]) == 0) or str(i[mapping[f]]) == "":
                        skip = True
                r.append(str(i[mapping[f]]))
            if not skip:
                data.append(r)
        return data

    def __findParents(self, endNodes):
        """for each end node in endNodes function finds all nodes on the way up to the root"""
        res = []
        for n in endNodes:
            tmp = n
            res.append(tmp)
            for i in range(n.count(".")):
                tmp = tmp.rstrip("1234567890").rstrip(".")
                if(tmp not in res):
                    res.append(tmp)
        """for i in res:
			nn = self.toID[self.toName[i]]
			for t in nn:
				if not t in res:
					res.append(t) """
        return res

    def __findMeshAttribute(self, data):
        """
        Function tries to find attribute which contains list of MeSH terms.
        """
        # we get a list of attributes
        dom = data.domain.attributes
        for i in dom:			  # for each attribute
            if i.varType == 6:
                for k in data:		   # for each domain
                    att = str(i.name)
                    try:										 # we try to use eval()
                        r = eval(str(k[att].value))
                        if type(r) == list:		 # attribute type should be list
                            if self.dataLoaded:		 # if ontology is loaded we perform additional test
                                for i in r:
                                    if self.toID.has_key(i): return att
                            else:				   # otherwise we return list attribute
                                return att
                    except SyntaxError:
                        continue
                    except NameError:
                        continue
        print "Program was unable to determinate MeSH attribute."
        return "Unknown"

    def __isPrecedesor(self, a, b):
        """ function returns true if in Mesh ontology exists path from term id a to term id b """
        if b.count(a) > 0:
            return True
        return False

    def __calculateAll(self, callback, fuzzy):
        """calculates all statistics"""
        # we build a dictionary 		meshID -> [description, noReference,noCluster, enrichment, deprivement, [cids] ]
        # print "ok"
        self.statistics = dict()
        if fuzzy:
            n = 0
            for i in self.reference:
                n = n + float(i['u'])
            cln = 0
            for i in self.cluster:
                cln = cln + float(i['u'])

        else:
            n = len(self.reference) 										# reference size
            cln = len(self.cluster)											# cluster size
        # frequency from reference list
        r = 0.0
        for i in self.reference:
            if callback and r % 10 == 0:
                r += 1.0
                callback(int(100 * r / (n + cln)))
            try:
                endNodes = list(set(eval(i[self.ref_att].value)))  # for every CID we look up for end nodes in mesh. for every end node we have to find its ID
            except SyntaxError:					 # where was a parse error
                print "Error in parsing ", i[self.ref_att].value
                if fuzzy:
                    n = n - float(i["u"])
                else:
                    n = n - 1
                continue
            # we find ID of end nodes
            endIDs = []
            for k in endNodes:
                if(self.toID.has_key(k)):					# this should be always true, but anyway ...
                    endIDs.extend(self.toID[k])
                else:
                    print "Current ontology does not contain MeSH term ", k, "."
            # endIDs may be empty > in this case we can skip this example
            if len(endIDs) == 0:
                if fuzzy:
                    n = n - float(i["u"])
                else:
                    n = n - 1
                continue
            # we find id of all parents
            # allIDs = self.__findParents(endIDs) PATCH
            allIDs = endIDs
            for k in list(set(allIDs)):				# for every meshID we update statistics dictionary
                if(not self.statistics.has_key(k)):			# first time meshID
                    self.statistics[k] = [0, 0, 0.0, 0.0]
                if fuzzy:
                    self.statistics[k][0] += float(i["u"])
                else:
                    self.statistics[k][0] += 1  # increased noReference
        # frequency from cluster list
        r = 0.0
        for i in self.cluster:
            try:
                if callback and r % 10 == 0:
                    r += 1.0
                    callback(int(100 * r / (n + cln)))
                endNodes = list(set(eval(i[self.clu_att].value)))  # for every CID we look up for end nodes in mesh. for every end node we have to find its ID
            except SyntaxError:
                # print "Error in parsing ",i[self.clu_att].value
                if fuzzy:
                    cln = cln - float(i["u"])
                else:
                    cln = cln - 1
                continue
            # we find ID of end nodes
            endIDs = []
            for k in endNodes:
                if(self.toID.has_key(k)):
                    endIDs.extend(self.toID[k])  # for every endNode we add all corensponding meshIDs
            # endIDs may be empty > in this case we can skip this example
            if len(endIDs) == 0:
                if fuzzy:
                    cln = cln - float(i["u"])
                else:
                    cln = cln - 1
                continue
            # we find id of all parents
            #allIDs = self.__findParents(endIDs)
            allIDs = endIDs
            for k in list(set(allIDs)):  # for every meshID we update statistics dictionary
                if self.statistics.has_key(k):
                    if fuzzy:
                        self.statistics[k][1] += float(i["u"])
                    else:
                        self.statistics[k][1] += 1  # increased noCluster
        self.ratio = float(cln) / float(n)
        # enrichment
        for i in self.statistics.iterkeys():
            self.statistics[i][2] = HYPERG.p_value(int(self.statistics[i][1]), int(n), int(cln), int(self.statistics[i][0]))
            self.statistics[i][3] = float(self.statistics[i][1]) / float(self.statistics[i][0]) / self.ratio   # fold enrichment
        self.calculated = True

    def __loadOntologyFromDisk(self):
        """
        Function loads MeSH ontology and chemical annotation into internal data structures.
        """
        self.toID = dict()  # name -> [IDs] Be careful !!! One name can match many IDs!
        self.toName = dict()  # ID -> name
        self.toDesc = dict()  # name -> description
        self.fromCID = dict()  # cid -> term id
        self.fromPMID = dict()  # pmid -> term id

        d = file(orngServerFiles.localpath_download('MeSH', 'mesh-ontology.dat'))
        f = file(orngServerFiles.localpath_download('MeSH', 'cid-annotation.dat'))

        # loading ontology graph
        t = 0
        for i in d:
            t += 1
            parts = i.split("\t")  # delimiters are tabs
            if(len(parts) != 3):
                print "error reading ontology ", parts[0]
            parts[2] = parts[2].rstrip("\n\r")
            ids = parts[1].split(";")
            self.toID[parts[0]] = ids  # append additional ID
            self.toDesc[parts[0]] = parts[2]
            for r in ids:
                self.toName[r] = parts[0]
            # loading cid -> mesh
        for i in f:
            parts = i.split(";")		# delimiters are tabs
            if(len(parts) != 2):
                print "error reading ontology ", parts[0]
            parts[1] = parts[1].rstrip("\n\r")
            cid = int(parts[0])
            if self.fromCID.has_key(cid):
                self.fromCID[cid].append(parts[1])
            else:
                self.fromCID[cid] = [parts[1]]
        # loading pmid -> mesh, TODO
        # print "Current MeSH ontology contains ", t, " mesh terms."
        return True


class pubMedHandler(ContentHandler):

    def __init__(self):
        self.state = 0  # 0 start state, 1 pmid, 2 title, 3 abstract, 4 mesh
        self.articles = []
        self.pmid = "0"
        self.title = ""
        self.mesh = list()
        self.abstract = ""
        self.affiliation = ""

    def startElement(self, name, attributes):
        # print "parsam ", name
        if name == "PubmedArticle":
            self.pmid = ""
            self.abstract = ""
            self.title = ""
            self.mesh = []
        if name == "PMID":
            self.state = 1
        if name == "ArticleTitle":
            self.state = 2
        if name == "AbstractText":
            self.state = 3
        if name == "DescriptorName":
            self.state = 4
            self.mesh.append("")
        if name == "Affiliation":
            self.state = 5

    def characters(self, data):
        if self.state == 1:
            self.pmid += data
        if self.state == 2:
            self.title += data.encode("utf-8")
        if self.state == 3:
            self.abstract += data.encode("utf-8")
        if self.state == 4:
            self.mesh[-1] += data.encode("utf-8")
        if self.state == 5:
            self.affiliation += data.encode("utf-8")

    def endElement(self, name):
        # print "   koncujem ", name
        self.state = 0
        if name == "PubmedArticle":
            self.articles.append([self.pmid, self.title, self.abstract, self.mesh, self.affiliation])


class MappedMeSHParser(SGMLParser):

    def reset(self):
        self.pieces = []
        self.terms = []
        self.foundMeSH = False
        self.nextIsTerm = False
        self.endTags = ['Display', 'Write to the Help Desk', 'Previous Indexing:', 'Entry Terms:', 'Pharmacologic Action:']
        SGMLParser.reset(self)

    def unknown_starttag(self, tag, attrs):
        strattrs = "".join([' %s="%s"' % (key, value) for key, value in attrs])
        if self.foundMeSH and tag == 'a':
            self.nextIsTerm = True

    def handle_data(self, text):
        text = text.strip()
        if text == '':
            return
        if text == 'Heading Mapped to:':
            self.foundMeSH = True
        if self.endTags.count(text) > 0:
            self.foundMeSH = False
        elif self.nextIsTerm:
            self.terms.append(text)
            self.nextIsTerm = False


class PubChemMeSHParser(SGMLParser):

    def reset(self):
        self.next = 0
        self.nextLink = ''
        self.directTerms = []
        self.indirectTerms = []
        self.foundMeSH = False
        self.foundIndirectMeSH = 0
        SGMLParser.reset(self)
        # strategy as follows
        # Beetween strings "Drug and Chemical Info" and ("Pharmalogical Action" or "PubMed via MeSH" or "PubMed MeSH Keyword Summary") find hyperlinks. Based on title attribute we can distingue direct MeSH terms beetween mapped terms.

    def unknown_starttag(self, tag, attrs):
        # if self.foundMeSH:
        # print 'tag ', tag, attrs
        if tag == 'a' and len(attrs) > 2 and attrs[0][0] == 'name' and attrs[0][1] == 'gomesh':  # print attrs
            self.nextLink = attrs[1][1]

    def handle_data(self, text):
        text = text.strip()
        if text == '':
            return
        # print 'data ', text
        if self.foundMeSH:
            # print '*'
            self.directTerms.append(text)

        if self.foundIndirectMeSH == 2:
            # print "indirect term", text, self.nextLink
            self.indirectTerms.append((text, self.nextLink))
            self.foundIndirectMeSH = 0

        if self.foundIndirectMeSH == 1:
            self.foundIndirectMeSH = 2

        if text == 'Drug and Chemical Information:':
            self.foundIndirectMeSH = 1

        if text == 'Pharmacological Action' or text == 'Medication Information':
            self.foundIndirectMeSH = 0
        if text == "Chemical Classification":
            self.foundMeSH = True
        elif (text == "Safety and Toxicology" or text == "Literature" or text == "Classification" or text == "PubMed via MeSH" or text == "PubMed MeSH Keyword Summary"):
            self.foundMeSH = False
            if len(self.directTerms) > 0:
                self.directTerms.pop()


class SmilesParser(SGMLParser):

    def reset(self):
        self.nextSmile = False
        self.smiles = ""
        SGMLParser.reset(self)

    def unknown_starttag(self, tag, attrs):
        # print tag, " ", attrs
        if tag == "item":
            for (i, j) in attrs:
                if i == "name" and j == "CanonicalSmile":
                    self.nextSmile = True

    def handle_data(self, text):
        if self.nextSmile:
            self.smiles = text
            self.nextSmile = False


class CIDsParser(SGMLParser):

    def reset(self):
        self.cid = False
        self.cids = []
        SGMLParser.reset(self)

    def unknown_starttag(self, tag, attrs):
        # print tag, " ", attrs
        if tag == "id":
            self.cid = True

    def unknown_endtag(self, tag):
        # print tag,
        if tag == "id":
            self.cid = False

    def handle_data(self, text):
        # print text
        if self.cid and not text == '\n\t\t':
            # print int(text.strip())
            self.cids.append(int(text))


class FormulaParser(SGMLParser):

    def reset(self):
        self.nextFormula = False
        self.formula = ""
        SGMLParser.reset(self)

    def unknown_starttag(self, tag, attrs):
        # print tag, " ", attrs
        if tag == "item":
            for (i, j) in attrs:
                if i == "name" and j == "MolecularFormula":
                    self.nextFormula = True

    def handle_data(self, text):
        if self.nextFormula:
            self.formula = text
            self.nextFormula = False


class CIDSParser(SGMLParser):

    def reset(self):
        self.nextCID = False
        self.cids = []
        SGMLParser.reset(self)

    def unknown_starttag(self, tag, attrs):
        # print tag, " ", attrs
        if tag == "item":
            for (i, j) in attrs:
                if i == "name" and j == "int":
                    self.nextCID = True

    def handle_data(self, text):
        if self.nextCID:
            self.cids.append(text)
            self.nextCID = False


class pubChemAPI(object):

    """
    Class pubChemAPI is used to simplift the access to the pubChem database.
    """

    def __init__(self):
        self.data = ""
        self.identifier = ""
        self.database = ""

    def getSMILE(self, id, typ):
        """
        Function getSMILE takes chemicals id and its type and returns a SMILES representation.
        """
        dbs = {"sid": "pcsubstance", "cid": "pccompound"}
        if not dbs.has_key(typ):
            return "Unknown identifier"
        # maybe we already have the data ...
        if id == self.identifier and dbs[typ] == self.database:
            raw = self.data
        else:
            url2fetch = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=" + dbs[typ] + "&id=" + str(id)
            # print url2fetch
            d = urlopen(url2fetch)
            raw = d.read()
            self.data = raw
            self.identifier = id
            self.database = dbs[typ]
        pr = SmilesParser()
        pr.feed(raw)
        data = pr.smiles
        pr.close()
        return data

    def getMeSHterms(self, cid):
        """
        Functions getMeSHterms tries to find MeSH annotation for given the CID in the pubChem database.
        """
        usock = urlopen("http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi?cid=" + str(cid) + "&viewopt=PubChem&hcount=100&tree=t#MeSH")
        parser = PubChemMeSHParser()
        text = usock.read()
        # print text
        parser.feed(text)
        usock.close()
        parser.close()
        allTerms = parser.directTerms
        # print allTerms
        for (k, i) in parser.indirectTerms:
            # print k, i
            # continue
            usock = urlopen(i)
            parser = MappedMeSHParser()
            parser.feed(usock.read())
            allTerms.extend(parser.terms)
            usock.close()
            parser.close()
        # from allTerms make a set
        allTerms = list(set(allTerms))
        return allTerms

    def getMolecularFormula(self, id, typ):
        """
        Functions getMolecularFormula tries to find molecular formula for the given the chemical id and type in the pubChem database.
        """
        dbs = {"sid": "pcsubstance", "cid": "pccompound"}
        if not dbs.has_key(typ):
            return "Unknown identifier"
        # maybe we already have the data ...
        if id == self.identifier and dbs[typ] == self.database:
            raw = self.data
        else:
            url2fetch = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=" + dbs[typ] + "&id=" + str(id)
            # print url2fetch
            d = urlopen(url2fetch)
            raw = d.read()
            self.data = raw
            self.identifier = id
            self.database = dbs[typ]
        pr = FormulaParser()
        pr.feed(raw)
        data = pr.formula
        pr.close()
        return data

    def getCIDfromSID(self, sid):
        """
        Functions getCIDfromSID converts compound id to substance id.
        """
        dbs = {"sid": "pcsubstance", "cid": "pccompound"}
        if not dbs.has_key("sid"):
            return "Unknown identifier"
        # maybe we already have the data ...
        if sid == self.identifier and dbs["sid"] == self.database:
            raw = self.data
        else:
            url2fetch = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=" + dbs["sid"] + "&id=" + str(sid)
            # print url2fetch
            d = urlopen(url2fetch)
            raw = d.read()
            self.data = raw
            self.identifier = sid
            self.database = dbs["sid"]
        pr = CIDSParser()
        pr.feed(raw)
        data = pr.cids
        pr.close()
        return data

    def getPharmActionList(self, cid):
        """
        Functions getPharmActionList finds pharmacological MeSH annotation for the given CID.
        <Item Name="PharmActionList" Type="List">
        """
        dbs = {"sid": "pcsubstance", "cid": "pccompound"}
        # maybe we already have the data ...
        if cid == self.identifier and dbs["cid"] == self.database:
            raw = self.data
        else:
            url2fetch = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=" + dbs["cid"] + "&id=" + str(cid)
            d = urlopen(url2fetch)
            raw = d.read()
        # print url2fetch
        self.data = raw
        self.identifier = cid
        self.database = dbs["cid"]
        ndata = raw.replace('\t', '').replace('\n', '')
        xmldoc = minidom.parseString(ndata)
        items = xmldoc.getElementsByTagName('Item')
        data = []
        for i in items:
            if i.attributes['Name'].value == 'PharmActionList':
                terms = i.childNodes
                for t in terms:
                    data.append(str(t.childNodes[0].data))
        return data

    def getCIDs(self, name, weight):
        """
        Functions getCIDs tries to find correct chemical id based on given chemical name and weight.
        """
        offset = 0.005
        url2fetch = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pccompound&term="
        url2fetch += name.replace(" ", "%20") + "%20" + str(weight * (1 - offset)) + ":" + str(weight * (1 + offset)) + "[mw]&retmax=10000"
        # print url2fetch
        d = urlopen(url2fetch)
        raw = d.read()
        e = CIDsParser()
        e.feed(raw)
        return e.cids

    def getCIDs(self, name):
        """
        Functions getCIDs tries to find correct chemical id based on given chemical name.
        """
        # print "k"
        url2fetch = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pccompound&term="
        url2fetch += name.replace(" ", "%20") + "&retmax=10000"
        # print url2fetch
        d = urlopen(url2fetch)
        raw = d.read()
        e = CIDsParser()
        e.feed(raw)
        return e.cids

    def __strip_ml_tags__(self, in_text):
        s_list = list(in_text)
        i, j = 0, 0
        while i < len(s_list):
            if s_list[i] == '<':
                while s_list[i] != '>':
                    s_list.pop(i)
                s_list.pop(i)
            else:
                i = i + 1
        join_char = ''
        return join_char.join(s_list)
