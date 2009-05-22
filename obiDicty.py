import sys, pprint, time, re
from itertools import *
import urllib2
import orange
import socket
import os

#utility functions - from Marko's mMisc.py

def splitN(origsize, maxchunk):
    """
    Splits an integer into chunks of given size. Each created chunk
    except possibly the last one is of maximum allowed size.
    Chunks are returned in list.
    """
    l = [maxchunk]*(origsize/maxchunk)
    a = origsize % maxchunk
    if a > 0: 
        l.append(a)
    return l

def split(l, maxchunk):
    """
    Splits list l into chunks of size maxchunk. Each created chunk
    except possibly the last one is of maximum allowed size.
    """
    sizes = splitN(len(l), maxchunk)
    out = []
    tillNow = 0
    for s in sizes:
        out.append(l[tillNow:tillNow+s])
        tillNow += s
    return out       

def lloc(l,n):
    """
    List location in list of list structure.
    Enable the use of negative locations:
    -1 is the last element, -2 second last...
    """
    if n < 0:
        return len(l[0])+n
    else:
        return n

def loc(l,n):
    """
    List location.
    Enable the use of negative locations:
    -1 is the last element, -2 second last...
    """
    if n < 0:
        return len(l)+n
    else:
        return n

def nth(l,n):
    """
    Returns only nth elemnt in a list.
    """
    n = lloc(l,n)
    return [ a[n] for a in l ]

def imnth(l, ns):
    """
    Return only columns as specified in ns. Returns an generator.
    """
    ns = [ lloc(l,n) for n in ns ]
    for a in l:
        yield [ a[n] for n in ns ]

def flatten(l,r=0):
    """
    Flatten a python structure into a list. Leave strings alone.
    """
    if type(l) == type("a"):
        return [ l ]
    try: #if enumerable then flatten it's elements
        rec = [ flatten(a,r=r+1) for a in l ]
        ret = []
        for a in rec:
            ret = ret + a
        return ret
    except:
        return [ l ]

def mxrange(lr):
    """
    Multiple xranges. Can be used to traverse matrices.
    This function is very slow due to unknown number of
    parameters.

    >>> mxrange([3,5]) 
    [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2)]
    >>> mxrange([[3,5,1],[9,0,-3]])
    [(3, 9), (3, 6), (3, 3), (4, 9), (4, 6), (4, 3)]
    """
    if len(lr) == 0:
        yield ()
    else:
        #it can work with single numbers
        index = lr[0]
        if type(1) == type(index):
            index = [ index ]
        for a in range(*index):
            for b in mxrange(lr[1:]):
                yield tuple([a] + list(b))


def issequencens(x):
    """
    Is x a sequence and not string ? We say it is if it has a __getitem__ 
    method and it is not an instance of basestring.
    """
    return hasattr(x, '__getitem__') and not isinstance(x, basestring)

#end utility functions

socket.setdefaulttimeout(30)

verbose = 0

def median(l):
    if len(l) == 0:
        return None
    l = sorted(l)
    if len(l) % 2 == 1: #odd
        return l[len(l)/2]
    else: #even
        return (l[len(l)/2-1] + l[len(l)/2])/2.0

class HttpGetException(Exception): pass

def replaceChars(address):
    return address.replace(" ", "%20")

def httpGet(address, *args, **kwargs):
    if verbose: 
        print address
    address = replaceChars(address)
    t1 = time.time()
    f = urllib2.urlopen(address, *args, **kwargs)
    """
    if verbose:
        print "Info"
        print f.info()
    """
    read = f.read()
    if verbose:
        print "bytes", len(read)
    if verbose:
        print time.time() - t1
    return read

def txt2ll(s, separ=' ', lineSepar='\n'):
    return [ a.split(separ) for a in s.split(lineSepar) ]

class DBInterface(object):
 
    def __init__(self, address, *args, **kwargs):
        self.address = address
        self.args = args
        self.kwargs = kwargs

    def raw(self, request, tryN=3):
        if verbose:
            print "tryN", tryN

        if tryN == 0:
            return None

        try:
            return httpGet(self.address + request, *self.args, **self.kwargs)
        except IOError:
            return self.raw(request, tryN=tryN-1)

    def get(self, request, tryN=3):
        rawf = self.raw(request)
        if rawf == None:
            raise Exception("Connection error when contacting " + self.address + request)
        if rawf[:1] == "<": #an error occurred - starting some html input
            #TODO are there any other kinds of errors?
            if tryN > 0:
                print "trying again"
                return self.get(request, tryN=tryN-1)
            else:
                raise Exception("Error with the database")


        a = txt2ll(rawf, separ='\t')
        
        if a[-1][0] == "": #remove empty line on end
            a = a[:-1]
        return a

def _test():
    import doctest
    doctest.testmod()

def splitTableOnColumn(ll,n):
    omap = {}
    n = lloc(ll,n)
    for l in ll:
        cell = omap.get(l[n], [])
        cell.append(l)
        omap[l[n]] = cell
    return omap

def neededColumns(legend, want):
    return [ legend.index(a) for a in want ]

def onlyColumns(ll, legend, want):
    return list(imnth(ll, neededColumns(legend, want)))

def ll2dic(ll, key=0, value=1):
    """
    Converts LL to map. Key is key position, value value position
    """
    names = nth(ll, lloc(ll, key))
    if len(names) == len(set(names)):
        return dict(list(imnth(ll, [ lloc(ll, key),  lloc(ll, value) ] )))
    else:
        raise Exception("all keys are not unique")

def dic2ll(dic):

    columns = sorted(dic.values()[0].keys())
    ll = []
    for key, d in dic.items():
        ll.append([ key ] + [ d[a] for a in columns ])
    return columns,ll

def allUnique(els):
    return len(els) == len(set(els))

def reverseDic(d):
    """
    Create a reverse dictionary if a unique reverse is possible.
    """
    if allUnique(d.values()):
        return dict([ (b,a) for a,b in d.items() ])

def chainLookup(a, dics, force=[]):
    """
    Goes through a list of dictionaries. On each step try to
    transform current key to value. The values becomes the key
    for next search If unsuccessfull, query the next
    dictionary with current key. Force last translation.
    """
    force = [ loc(dics, i) for i in force ]
    for i,dic in enumerate(dics):
        if i in force: #force specified translations
            a = dic[a]
        else:
            if a in dic:
                a = dic[a]
    return a

class DatabaseConnection(object):

    """
    Type is object id
    """
    
    aoidPairs = txt2ll("""time extractions.developmental_time_point
sample biological_samples.sample
growthCond biological_samples.growth_condition
treatment biological_samples.treatment
replicate biological_sample_replicates.replicate
techReplicate chips.replicate
platform chips.chip_platform
isTimeSeries biological_sample_replicates.is_time_series""")

    obidPairs = txt2ll("""norms normalizations
samples biological_samples
replicates biological_sample_replicates
analysis analysis
experiments experiments
extractions extractions
chips chips""")

    def __init__(self, address, buffer=None):
        self.address = address
        self.db = DBInterface(address)
        self.buffer = buffer
        self.preload()

    def preload(self):

        # aoids are mappings from annotation name to annotation id
        self.aoids = ll2dic(self.__annotationTypes(), 1, 0)
        self.saoids = ll2dic(self.aoidPairs, 0, 1)
        self.aoidsr = reverseDic(self.aoids)
        self.saoidsr = reverseDic(self.saoids)

        # obids are mappings from object id to annotation id
        self.obids = ll2dic(self.__objects(), 1, 0)
        self.sobids = ll2dic(self.obidPairs, 0, 1)
        self.obidsr = reverseDic(self.obids)
        self.sobidsr = reverseDic(self.sobids)

    def aoidt(self, s):
        return chainLookup(s, [self.saoids, self.aoids], force=[-1])

    def obidt(self, s):
        return chainLookup(s, [self.sobids, self.obids], force=[-1])
 
    def aoidtr(self, s, **kwargs):
        return chainLookup(s, [self.aoidsr, self.saoidsr], **kwargs)

    def obidtr(self, s):
        return chainLookup(s, [self.obidsr, self.sobidsr])

    def pq(self, q):
        """
        Prepare query. 
        ||| separator between conditions, 
        *** denotes equivalence
        """
        o =  "|||".join([ self.aoidt(a) + "***" + b for a,b in q.items()])
        return o

    def downloadMulti(self, command, ids, chunk=100, transformfn=None, separatefn=None):
        """
        Results in the same order as in ids.
        """

        sids = split(ids,chunk)
    
        def bufcommand():
            if transformfn:
                return "TRANS " + command
            else:
                return command

        for i,sidp in enumerate(sids):

            buffered = []
            unbuffered = []
        
            for a in sidp:
                if self.inBuffer(bufcommand().replace("$MULTI$", a)):
                    buffered.append(a)
                else:
                    unbuffered.append(a)

            res = []
            legend = []

            if len(unbuffered) > 0:
                res, legend = self.sq(command.replace("$MULTI$", ",".join(unbuffered)),\
                    buffer=False)
            else:
                # get legend from buffer also
                legend = self.fromBuffer(bufcommand().replace("$MULTI$", buffered[0]))[0]

            #split on different values of the first column - first attribute

            if not separatefn:
                antss = splitTableOnColumn(res, 0)
            else:
                legend, antss = separatefn([legend]+res)

            #if transform before saving is requested, do it
            if transformfn:
                nantss = {}
                nlegend = None
                for a,b in antss.items():
                    nb, nlegend = transformfn(b, legend)
                    nantss[a] = nb
                legend = nlegend
                antss = nantss
 
            #here save buffer
            for a,b in antss.items():
                self.toBuffer(bufcommand().replace("$MULTI$", a), [ legend ] + b)

            antssb = dict([ (b, self.fromBuffer(bufcommand().replace("$MULTI$", b))[1:]) for b in buffered ])
            antss.update(antssb)

            #put results in order
            tl = []
            for ci in sidp:
                yield antss[ci], legend

    def geneInfo(self):
        res,legend = self.sq("action=gene_info")
        return res, legend

    def annotationOptions(self, ao=None, onlyDiff=False, **kwargs):
        """
        Returns annotation options for given query. Returns all possible 
        annotations if the query is omitted.

        If ao is choosen, only result
        """
        params = ""
        if len(kwargs) > 0: params += "&query=" + self.pq(kwargs)
        if ao: params += "&annotation_object_id=" +  self.aoidt(ao)
        res,legend = self.sq("action=get_annotation_options%s" % (params), bufadd=self.address)
        res = onlyColumns(res, legend, ['annotation_object_id', 'value' ])

        #join columns with the annotation object id
        joined = {}
        for key,v in res:
            key = self.aoidtr(key)
            cur = joined.get(key, [])
            cur.append(v)
            joined[key] = cur

        if onlyDiff:
            joined = dict([ (a,b) for a,b in joined.items() if len(b)>1 ])

        return dict([ (a, sorted(b)) for a,b in joined.items() ])

    
    def inBuffer(self, addr):
        if self.buffer:
            return self.buffer.contains(self.address + addr)
        else:
            return False

    def fromBuffer(self, addr):
        return self.buffer.get(self.address + addr)

    def toBuffer(self, addr, cont):
        if self.buffer:
            return self.buffer.add(self.address + addr, cont)

    def bufferFun(self, bufkey, fn, *args, **kwargs):
        """
        If bufkey is already present in buffer, return its contents.
        If not, run function with arguments and save its result
        into the buffer.
        """
        if self.inBuffer(bufkey):
            res = self.fromBuffer(bufkey)
        else:
            res = fn(*args, **kwargs)
            self.toBuffer(bufkey, res)
        return res

    def sq(self, s1, buffer=True, bufadd=""):
        if buffer:
            res = self.bufferFun(bufadd + s1, self.db.get, s1)
        else:
            res = self.db.get(s1)
        return res[1:],res[0]

    def annotation(self, type, id):
        return list(self.annotations(type, [ id ]))[0]

    def meaningfulAnnot(self, name):
        if name in self.saoids:
            return True
        else:
            return False

    def keepOnlyMeaningful(self, annot):
        """
        Keep only meaningful annotations
        """
        if type(annot) == type({}):
            return dict( [ (a,b) for a,b in annot.items() \
                if self.meaningfulAnnot(a) ] )
        else:
            return [ [ a,b ] for a,b in annot \
                if self.meaningfulAnnot(a) ]


    def annotations(self, type, ids=None, all=False):
        """
        Returns a generator returning annotations for specified type and ids. 
        If ids are left blank, all annotations are outputed. Annotations are in the same order
        as input ids.
        If all is True, all annotations are kept, else keep only "meaningful".
        """
        
        inputids = False
        if ids != None:
            inputids = True
            antss = self.downloadMulti(
                "action=get_annotations&ids=$MULTI$&object_id=%s" 
                % (self.obidt(type)), ids)
        else:
            res,legend = self.sq(
                "action=get_annotations&object_id=%s"
                % (self.obidt(type)))
            antss = splitTableOnColumn(res, 0)
            ids = nth(antss.items(),0)
            antss = zip(nth(antss.items(),1), [ legend ]*len(antss))

        for ants in izip(antss,ids):
            (res, legend), id = ants
            res2 = onlyColumns(res, legend, ['name', 'value'])
            res2 = [ [ self.aoidtr(a),b ] for a,b in res2 ]
            if not all:
                res2 = self.keepOnlyMeaningful(res2)
            if inputids:
                yield res2
            else:
                yield (id, res2)

    def search(self, type, **kwargs):
        """
        Break search for multiple values of one attribute to independant searches.
        Search is case insensitive.
        
        List of searchable annotation types: self.saoids.keys()

        example usage:
        search("norms", platform='minichip', sample='abcC3-') 
            finds all ids of normalized entries where platform is minchip and 
            sample is abcC3-
        search("norms", platform='minichip', sample=[ 'abcC3-', 'abcG15-']) 
            finds all ids of normalized entries where platform is minichip and 
            sample is abcC3- or those where platform is minichip and sample
            is abcG15-
        """
        
        #search for all combinations of values - this is slow!

        l = []
        for k,v in kwargs.items():
            if not issequencens(v):
                v = [ v ]
            l.append((k,v))

        ares = []

        for r in mxrange([ len(v) for k,v in l ]):
            dico = {}
            for i,a in enumerate(r):
                dico[l[i][0]] = l[i][1][a]

            res,_ = self.sq("action=search&object_id=%s&query=%s" \
                % (self.obidt(type), self.pq(dico)), bufadd=self.address)

            ares += res

        return sorted(set(nth(ares, 0)))


    def chipN(self, id):
        return list(self.chipNs([id]))[0]

    def chipR(self, id):
        return list(self.chipRs([id]))[0]
  
    def chipNs(self, ids):
          
        def sel(res, legend):
            #Drop unwanted columns - for efficiency
            res = onlyColumns(res, legend, ["spot_id", 'M'])
            legend = onlyColumns( [ legend ], legend, ["spot_id", 'M'])[0]
            return res, legend

        antss = self.downloadMulti("action=get_normalized_data&ids=$MULTI$", ids, chunk=2, transformfn=sel)
        for a,legend in antss:
            yield a   

    def chipNsN(self, ids, annots):
        """
        Download chips using new shorter format.
        """
        chip_map_ids = zip(ids,[ dict(a)['chips.chip_map_id'] for a in annots ])

        def separateByChipsMaps(l):
            begin = 0
            cm = l[0][1]
            cp = 0
            for id,m in l[1:]:
                cp += 1
                if m != cm:
                    yield l[begin:cp]
                    cm = m
                    begin = cp
            yield l[begin:cp+1]
        
        sep = list(separateByChipsMaps(chip_map_ids))
      
        def sel(res, legend):
            #Drop unwanted columns - for efficiency
            res = onlyColumns(res, legend, ["spot_id", 'M'])
            legend = onlyColumns( [ legend ], legend, ["spot_id", 'M'])[0]
            return res, legend

        def separatefn(res):
            #each one is own rown
            #genes are in the first row
            genes = res[0][1:]
            cids = nth(res,0)[1:]

            antss = {}
            for i,cid in enumerate(cids):
                row = i+1
                vals = res[row][1:]
                antss[cid] = [ list(a) for a in zip(genes, vals) ]
            return ['spot_id', 'M'], antss

        for part in sep:
            pids = nth(part,0)
            antss = self.downloadMulti("action=get_normalized_data&mergeexperiments=1&ids=$MULTI$", pids, chunk=10, transformfn=sel, separatefn=separatefn)
            for a, legend in antss:
                yield a

    def chipRs(self, id):
        antss = self.downloadMulti("action=get_raw_data&ids=$MULTI$", ids, chunk=2)
        for a,legend in antss:
            yield a
  
    def spotId(self):
        res,legend = self.sq("action=spot_id_mapping")
        res2 = onlyColumns(res, legend, ["spot_id", 'ddb_g', 'genename'])
        return res2

    def annotationTypes(self):
        """
        Returns all annotation types.
        """
        return self.aoids.keys()

    def objects(self):
        """
        Returns all objects.
        """
        return self.obids.keys()

    def __annotationTypes(self):
        """
        Returns list of [ annotation_object_id, name ]
        """
        res, legend = self.sq("action=get_annotation_types")
        res2 = onlyColumns(res, legend, ["annotation_object_id", 'name'])
        return res2

    def __objects(self):
        res, legend = self.sq("action=get_objects")
        res2 = onlyColumns(res, legend, ["object_id", 'object_name'])
        return res2

    def spotMap(self):
        spotids = self.spotId()
        spotmap = [ (a[0],a[1]) for a in spotids ]

        spotmapd = {}

        for a,b in spotmap:
            if a in spotmapd:
                spotmapd[a] = spotmapd[a] + "-" + b
            else:
                spotmapd[a] = b

        return spotmapd

    def dictionarize(self, ids, fn, *args, **kwargs):
        """
        Creates a dictionary from id: function result.
        Callback for each part done.
        """
        callback = kwargs.pop("callback", None)
        odic = {}
        for a,b in izip(ids, fn(*args, **kwargs)):
            odic[a] = b
            if callback: callback()
        return odic
        #return dict(zip(ids, list(fn(*args, **kwargs))))

    def groupAnnotations(self, annotations, prefix="", join=[], separate=None):
        """
        Groups annotations by joining annotation options in 'join' to
        the same group while keeping annotations with different values
        of options in 'separate' apart. If join is None, set join to all 
        annotation types. If separate is None, use all annotation types not 
        in join as separate.

        Returns list of tuples. Each tuple describes a group and contains:
            1. [ (name, chipid) ] - a list of tuples containing group element 
                name and it's chip id.
            2. dictionary of joined annotations of group elements.

        """

        if verbose:
            print "Grouping annotations"

        # Separate and join are not a regular expressions (because of group sorting).

        ignore=['.*']

        if join == None:
            join = self.saoids.keys()
    
        #If separate is None, use all annotation types not in join as separate.
        if separate == None:
            separate = sorted(set(self.saoids.keys())-set(join))

        annotationsIn = dict(annotations)

        if len(annotationsIn) == 0:
            return None

        #join should be top most specification
        topMost = [ self.aoidtr(a) for a in self.aoids.keys() ]

        def forceTopMost(l):
            for a in l:
                if a not in topMost:
                    raise Exception(a + " not top most annotation type." )
        forceTopMost(join)

        #all possible annotations
        annots = allAnnotationVals(annotationsIn.values())

        #separate are not ignored. others are ignored!
        ignoreNot = separate + join

        def isIgnored(a):
            for exp in ignoreNot:
                if re.match(exp, a):
                    return False
            for exp in ignore:
                if re.match(exp, a):
                    return True
            return False
    
        def ok(x):
            a,b = x
            if len(b) <= 1:
                return False
            if isIgnored(a):
                return False
            return True

        diffannots = dict(filter(ok, annots.items()))
        #print "DIFF ANNOTS", diffannots

        """
        Group those which have everything the same except of join while
        ignoring annotation which match meaningless + ignore.
        """
    
        def validAnnot(x):
            try: 
                fkdllfkd = self.aoidt(x)
                return True
            except:
                return False

        allowDiffSet = set(filter(validAnnot, join))
        danames = set(diffannots.keys())
        reallyDiff = danames & allowDiffSet
        leaveKeys = danames - reallyDiff

        #remove all non-neccessary keys form the annotation dictionary
        #and form tem in groups afterwards

        def dictWithOnly(dic, keys):
            """
            Returns a shallow copy of the dictionary with
            key not in keys removed.
            """
            odic = {}
            for a,b in dic.items():
                if a in keys:
                     odic[a] = b
            return odic

        candidates = dict([ (id, dictWithOnly(dict(d),leaveKeys) ) \
            for id,d in annotationsIn.items() ])

        def groupDicts(ldics):
            """
            Group same dictionaries together. Returns indices
            of original dictionaries
            """
            ldics = dict([(id,frozenset(d.items())) for id,d in ldics.items() ])
            
            dd = {}
            for i,d in ldics.items():
                group = dd.get(d, [])
                group.append(i)
                dd[d] = group

            return dd.values()

        groups = groupDicts(candidates)

        def joinedAnnotation(group):
            return allAnnotationVals([ annotationsIn[v] for v in group ])

        def csorted(l):
            """
            Return sorted list of string , but choose sorting type first:
            if all are integer, sort in integers and then return string
            """
            try:
                l = [ int(a) for a in l ]
            except:
                pass
            return [ str(a) for a in sorted(l) ]


        #sort the groups by their annotations
        def compareGroups(g1, g2):
            """
            Compare groups by consequtive comparison of annotatition values.
            Consider those in "separate" first.
            """

            #separated need to be top most!
            allannot = sorted(set(g1.keys()) - set(separate))
            allannot = separate + allannot

            for compnow in allannot:
                #try sorting values in lists as integers
                c = cmp(csorted(g1[compnow]), csorted(g2[compnow]))
                if c != 0:
                    return c

            return 0

        groups = sorted(groups, key=joinedAnnotation, cmp=compareGroups)

        """
        At this point the groups going to the same example table
        are already made. Now we need to control value order, example
        order and table order for those groups.

        Important decisions:
        - example order: by spot id
        - attribute order: ([1,2,3],[1,2]) -> [1,1],[1,2],[2,1],...
        - table order - comparison of annotation valeues
        """

        #respect specified order, but put sample in the front
        differentialNames = [ a for a in join if a in reallyDiff ]
        if 'sample' in differentialNames: #if sample is present, put it to front
            differentialNames.remove('sample')
            differentialNames = [ "sample" ] + differentialNames

        differentialVals = [ csorted(diffannots[a]) for a in differentialNames ]

        def difV(li):
            """
            Returns differential values for tuples of indices (input list).
            (i,j) returns the i-th value for first differential name, 
            and j-th values for second differential name.
            """
            #print "difV", li, [ differentialVals[i][v] for i,v in enumerate(li) ]
            return [ differentialVals[i][v] for i,v in enumerate(li) ]

        def nameInds(values, inds):
            """
            Create attribute name for given values. If there are multiple 
            repetitions of the same value, enumerate them.
            """
            def aux(i,a):
                #leave out sample annotation
                name = ','.join([ n + '=' + v if n != 'sample' else v for n,v in values ])
                if len(inds) > 1:
                    name += "_" +  str(i)
                return name

            return [ aux(i,a) for i,a in enumerate(inds) ]
 
        def extractNeeded(annotations, needValues):
            """
            Returns indices (input list) of annotation with all values as 
            specified in needValues.
            """
            ok = []
            for i,an in enumerate(annotations):
                sl = sum([ an[nam] != val for nam,val in needValues ])
                if sl == 0: 
                    ok.append(i)
            return ok

        exampleTables = []

        for group in groups:

            lens = [ len(a) for a in differentialVals ]
            groupnames = []
            groupids = []

            for valsi in mxrange([ b for b in lens ]):

                gannots = [ ll2dic(annotationsIn[g]) for g in group ]

                needValues = zip(differentialNames, difV(valsi))
                inds = extractNeeded(gannots, needValues)

                groupids = groupids + [ group[a] for a in inds ]
                groupnames = groupnames + nameInds(needValues, inds)

            #Add every annotation. If there are multiple annotation values
            #for the same annotation, add all annotation values.
            
            ca = joinedAnnotation(group)
            #ca = dict( [(a,csorted(b)) for a,b in ca.items()] )

            exampleTables.append((zip(groupnames, groupids), ca))

        return exampleTables

    def exampleTables(self, chipsm, groups, spotmap={}, annotations=None, callback=None, chipidname=False, exclude_constant_labels=False, annots={}, newannotations=False):
        """
        Create example tables from chip readings, spot mappings and 
        group specifications.

        groups are input from "groupAnnotations" function. 
        spotmap is a dictionary of { spotid: gene }
        """

        if verbose:
            print "Creating example tables"

        if annotations != None:
            annotations = dict(annotations)

        exampleTables = []

        ids = flatten([ nth(g[0], 1) for g in groups ])

        #where to put n-th occurance of some spot id
        def outPositionMap(ids):
            """
            Where to put n-th occurance of some spot id.
            Input: list of lists of ids.
            Output: map of id,repetition -> position
            """
            
            amap = {}

            for lids in ids:

                repeats = {}
                for id in lids:
                    rep = repeats.get(id, 0)
                    repeats[id] = rep+1
                    key = (id, rep)
                    amap[key] = amap.get(key, 0) + 1

            return dict( [ (k,i) for i,k in enumerate(sorted(amap.keys())) ] )

        posMap = outPositionMap([ nth(c,0) for c in chipsm.values() ])
        npos = len(posMap.values())

        if callback: callback()

        ddb = [ "" for a in range(npos) ] 

        for k,ind in posMap.items():
            ddb[ind] = spotmap.get(k[0], "#"+k[0]) 

        for group in groups:

            if verbose:
                print  "joining group", group

            groupnames = []
            groupvals = []
            groupannots = []

            pairs = group[0]

            annotc = None

            if len(group) > 1:
                annotc = group[1]

            for name,chipid in pairs:

                vals = [ None ] * npos

                repeats = {}
                putinds = []
                for id,v in chipsm[chipid]:
                    rep = repeats.get(id, 0)
                    repeats[id] = rep+1
                    key = (id, rep)
                    putind = posMap[key]
                    putinds.append(putind)
                    vals[putind] = v

                #regarding chipid names
                if chipidname:
                    groupnames.append(chipid) 
                else:
                    groupnames.append(name)

                if annots:
                    #add chip id to annotations
                    groupannots.append(annots[chipid]+ [['chipid', str(chipid)]])

                groupvals.append(vals)

            et = createExampleTable(groupnames, groupvals, groupannots, ddb, exclude_constant_labels=exclude_constant_labels)

            if newannotations:
                annotc = allAnnotationVals( [annots[v] for v in nth(pairs,1) ] )

            annotc["chipids"] = nth(group[0], 1)

            et.setattr("annot", annotc)
            exampleTables.append(et)
            
            if callback: callback()
    
        return exampleTables


    def getData(self, type="norms", join=["time"],  exclude_constant_labels=False, separate=None, average=median, ids=None, callback=None, chipidname=False, **kwargs):
        """
        Returns a list of examples tables for a given search query and post-processing
        instructions.

        Parameters: 
            average: function used for combining multiple reading of the same spot on
                a chip. If None, no averaging is done. Fuction should take a list
                of floats and return an "averaged" float.
            join: a list of annotation types which can be different in a single example
                table. Chips are grouped in groups, which can contain chips which have
                same annotations.
            separate: annotation types by which we consider groups as separate ones.
                If blank, take all those not in join.
            ids: a list of chip ids. If present, use this ids instead of making
                a search.
            exclude_constant_labels: if a label has the same value in whole 
            example table, remove it

        Defaults: Median averaging. Join by time.
        """


        def optcb():
            if callback: callback()

        cbc = CallBack(1, optcb, callbacks=10)

        if not ids:
            #returns ids of elements that match the search function
            ids = self.search(type, **kwargs)

        cbc.end()

        #downloads annotations
        cbc = CallBack(len(ids), optcb, callbacks=10)

        readall = self.dictionarize(ids, self.annotations, type, ids, all=True, callback=cbc)

        read = {}
        for a,b in readall.items():
            read[a] = self.keepOnlyMeaningful(b)

        annotsinlist = [] #annotations in the same order
        for id in ids:
            annotsinlist.append(readall[id])

        if verbose:
            print zip(ids,[ dict(a)['chips.chip_map_id'] for a in annotsinlist ])

        cbc.end()

        cbc = CallBack(1, optcb, callbacks=10)
        #make annotation groups
        etsa = self.groupAnnotations(read, join=join, separate=separate )
        cbc.end()

        #here could user intervent. till now downloads were small

        import time
        tstart =  time.time()

        #here download actually happens

        cbc = CallBack(len(ids), optcb, callbacks=999-50)
        if type == "norms":
            chipd = self.dictionarize(ids, self.chipNsN, ids, annotsinlist, callback=cbc)
            #chipd = self.dictionarize(ids, self.chipNs, ids, callback=cbc)
        else:
            chipd = self.dictionarize(ids, self.chipRs, ids, callback=cbc)

        cbc.end()
        
        #create example tables grouped according to user's wishes
        if verbose:
            print "DOWNLOAD TIME", time.time() - tstart

        cbc = CallBack(len(etsa)+1, optcb, callbacks=10)
        ets = self.exampleTables(chipd, etsa, spotmap=self.spotMap(), callback=cbc, chipidname=chipidname, annots=read, exclude_constant_labels=exclude_constant_labels)
        cbc.end()

        cbc = CallBack(len(ets), optcb, callbacks=10)

        #if average function is given, use it to join same spotids
        if average != None:
            etsa = []
            for et in ets:
                eta = averageAttributes(et, fn=average)
                eta.setattr("annot", et.annot)
                etsa.append(eta)
                cbc()
            ets = etsa

        cbc.end()

        return ets

    def get_single_data(self, type="norms", average=median, ids=None, 
            callback=None, exclude_constant_labels=False, **kwargs):
        """
        Get data in a single example table with labels of individual attributes
        set to annotations.
        """
        mtables = self.getData(type=type, join=None, separate=[], average=average, ids=ids, callback=callback, chipidname=False, exclude_constant_labels=exclude_constant_labels, **kwargs)
        if len(mtables) != 1:
            mtablesShouldHaveOnlyOneElementError()
        return mtables[0]

def allAnnotationVals(annots):
    """
    All annotation valuess for given annotations
    in a dict of { name: set of possible values } pairs.
    """
    av = {}
    for a in annots:
        for name,val in a:
            cvals = av.get(name, set([]))
            cvals.add(val)
            av[name] = cvals
    return av



def createExampleTable(names, vals, annots, ddb, cname="DDB", \
        exclude_constant_labels=False):
    """
    Create an ExampleTable for this group. Attributes are those in
    names. 
    """
    attributes = [ orange.FloatVariable(n, numberOfDecimals=3) \
        for n in names ]

    #exclusion of names with constant values
    annotsvals = allAnnotationVals(annots)
    oknames = set(annotsvals.keys())
    if exclude_constant_labels:
        oknames = set(nth(filter(lambda x: len(x[1]) > 1, 
            annotsvals.items()), 0))

    #print oknames

    for a,an in zip(attributes, annots):
        a.setattr("attributes", dict([(name,val) for name,val in an if name in oknames]))

    domain = orange.Domain(attributes, False)
    ddbv = orange.StringVariable(cname)
    id = orange.newmetaid()
    domain.addmeta(id, ddbv)

    examples = []
    for v,d in izip(izip(*vals), ddb):
        ex = orange.Example(domain, [ floatOrUnknown(a) for a in v ])
        ex[cname] = d
        examples.append(ex)

    return orange.ExampleTable(domain,examples)


def averageAttributes(data, joinc="DDB", fn=median):
    """
    Averages attributes with the same "join" parameter using
    specified function. Builds a now dataset. Order of
    "join" parameter stays the same with regard to first
    appearance.
    """

    if verbose:
        print "Averaging attributes"

    valueso = []
    valuess = set(valueso)

    attributes = [ a for a in data.domain.attributes ]
    domain = data.domain

    #accumulate values
    valuesddb = dict( [ (at,{}) for at in attributes ])

    for ex in data:
        #join attribute - ddb
        j = str(ex[joinc])

        if j not in valuess:
            valueso.append(j)
            valuess.add(j)

        for a in attributes:
            val = ex[a]
            l = valuesddb[a].get(j, [])
            if not val.isSpecial():
                l.append(val.native())
            valuesddb[a][j] = l

    #print "len valueso", len(valueso)
    #print sorted(set([ str(ex[join]) for ex in data ]))

    #apply function fn to each attribute

    """
    for i,at in enumerate(data.domain.attributes):
        print valuesddb[at]["DDB_G0282817"], "CI" + data.annot["chipids"][i]
    """

    for a in attributes:
        for n,v in valuesddb[a].items():
            valuesddb[a][n] = floatOrUnknown(fn(v))
            #if n == "DDB_G0282817": print valuesddb[a][n]

    #create a new example table reusing the domain
    examples = []
    for v in valueso:
        example = orange.Example(domain, \
            [ valuesddb[a][v] for a in attributes ] )
        example[joinc] = v
        examples.append(example)

    return orange.ExampleTable(domain, examples)

def floatOrUnknown(a):
    """
    Converts an element to float if possible.
    If now, output "?".
    """
    try:
        return float(a)
    except:
        return "?"

class CallBack():
    """
    Converts "allparts" callbacks into by "callbacks"
    specified number of callbacks of function fn.
    """

    def __init__(self, allparts, fn, callbacks=100):
        self.allparts = allparts
        self.lastreport = 0.00001
        self.getparts = 0
        self.increase = 1.0/callbacks
        self.callbacks = callbacks
        self.fn = fn
        self.cbs = 0

    def __call__(self):
        self.getparts += 1
        done = float(self.getparts)/self.allparts
        while done > self.lastreport + self.increase:
            self.lastreport += self.increase
            self.fn()
            self.cbs += 1

    def end(self):
        while self.cbs < self.callbacks:
            self.fn()
            self.cbs += 1


class BufferSQLite(object):

    def __init__(self, filename, compress=True):
        self.compress = compress
        self.filename = filename
        self.conn = self.connect()

    def clear(self):
        """
        Removes all entries in the buffer
        """
        self.conn.close()
        os.remove(self.filename)
        self.conn = self.connect()


    def connect(self):
        import sqlite3
        conn = sqlite3.connect(self.filename)
        c = conn.cursor()
        c.execute('''create table if not exists buf
        (address text primary key, time text, con blob)''')
        c.close()
        conn.commit()
        return conn

    def contains(self, addr):
        c = self.conn.cursor()
        c.execute('select address from buf where address=?', (addr,))
        lc = len(list(c))
        c.close()
        if lc == 0:
            return False
        else:
            return True

    def add(self, addr, con):
        import cPickle, zlib, sqlite3
        if verbose:
            print "Adding", addr
        c = self.conn.cursor()
        if self.compress:
            bin = sqlite3.Binary(zlib.compress(cPickle.dumps(con)))
        else:
            bin = sqlite3.Binary(cPickle.dumps(con))
        c.execute('insert into buf values (?,?,?)', (addr, "0", bin))
        c.close()
        self.conn.commit()

    def get(self, addr):
        import cPickle, zlib
        if verbose:
            print "getting from buffer", addr
            t = time.time()
        c = self.conn.cursor()
        c.execute('select con from buf where address=?', (addr,))
        ls = list(c)
        first = ls[0][0]
        if verbose:
            print time.time() - t
        if self.compress:
            rc = cPickle.loads(zlib.decompress(first))
        else:
            rc = cPickle.loads(str(first))
        c.close()

        if verbose:
            print time.time() - t
        return rc



import urllib2
import orngServerFiles
import pickle

def download_url(url, repeat=2):
    def do():
        return urllib2.urlopen(url)

    if repeat <= 0:
        do()
    else:
        try:
            return do()
        except:
            return download_url(url, repeat=repeat-1)

def empty_none(s):
    if s:
        return s
    else:
        return None

class DictyBase(object):

    domain = "dictybase"
    filename = "information_mappings.pck"
    tags = [ "Dictyostelium discoideum", "gene", "essential", "dictyBase" ] 
 
    @classmethod
    def version(cls):
        orngServerFiles.localpath_download(cls.domain, cls.filename)
        return orngServerFiles.info(cls.domain, cls.filename)["datetime"]
    
    @classmethod
    def download_information(cls):
        """ 
        Downloads gene information and parses it. 
        Returns a dictionary {ID: (name, synonyms, products)}
        """
        s = download_url("http://www.dictybase.org/db/cgi-bin/dictyBase/download/download.pl?area=general&ID=gene_information.txt").read()
        out = []
        for l in txt2ll(s, separ='\t', lineSepar='\n')[1:]:
            if len(l) == 4:
                id = l[0]
                name = l[1]
                synonyms = filter(None, l[2].split(", "))
                products = l[3]
                out.append((id, name, synonyms, products))
        return dict((a,(b,c,d)) for a,b,c,d in out)

    @classmethod
    def download_mappings(cls):
        """ 
        Downloads DDB-GeneID-UniProt mappings and parses them. 
        Returns a list of (ddb, ddb_g, uniprot) triplets.
        
        2009/04/07: ddb's appear unique
        """
        s = download_url("http://www.dictybase.org/db/cgi-bin/dictyBase/download/download.pl?area=general&ID=DDB-GeneID-UniProt.txt").read()
        out = []
        for l in txt2ll(s, separ='\t', lineSepar='\n')[1:]:
            if len(l) == 3:
                ddb = empty_none(l[0])
                ddb_g = empty_none(l[1])
                uniprot = empty_none(l[2])
                out.append((ddb, ddb_g, uniprot))
        return out

    @classmethod
    def pickle_data(cls):
        info = cls.download_information()
        mappings = cls.download_mappings()
        return pickle.dumps((info,mappings), -1)

    def __init__(self):
        fn = orngServerFiles.localpath_download(self.domain, self.filename)
        self.info, self.mappings = pickle.load(open(fn, 'rb'))

if __name__=="__main__":
    verbose = 0

    a = DictyBase()
    print len(a.info)

    #dbc = DatabaseConnection("http://asterix.fri.uni-lj.si/microarray/api/index.php?", buffer=BufferSQLite("../tmpbuf1233"))
    dbc = DatabaseConnection("http://purple.bioch.bcm.tmc.edu/~anup/index.php?", buffer=BufferSQLite("../tmpbufnew"))

    print dbc.annotationOptions()

    count = 0
    def cb():
        global count
        count += 1
        #print "CBBB", count

    et = dbc.get_single_data(sample=[ "tagA-", "pkaC-"], callback=cb, exclude_constant_labels=True)
    print et.domain
    print et.domain[0].attributes
    
    """
    ets = dbc.getData(sample="tagA-", callback=cb)
    for i,et in enumerate(ets):
        et.save("%s/T%d.tab" % ("multid2", i))
        print et.annot
    """

