## problem je v tem, da moras steti gen veckrat, za vsak tip Evidence enkrat
import cPickle, math

### misc
# need to replace with a "standard" function
binoms = {}
def binom(n,k):
    global binoms

    r = binoms.get( str( (n, k)), None)
    if r:
        return r

    if n==0 or k==0 or n==k:
        binoms[ str( (n, k))] = 1
        return 1
    else: 
        r = binom(n-1,k-1) + binom(n-1,k)
        binoms[ str( (n,k))] = r
        return r
###

evidenceTypesOrdered = [
'IMP',
'IGI',
'IPI',
'ISS',
'IDA',
'IEP',
'IEA',
'TAS',
'NAS',
'ND',
'IC'
]

evidenceTypes = {
'IMP': 'inferred from mutant phenotype',
'IGI': 'inferred from genetic interaction', ## [with <database:gene_symbol[allele_symbol]>]',
'IPI': 'inferred from physical interaction', ## [with <database:protein_name>]',
'ISS': 'inferred from sequence similarity', ## [with <database:sequence_id>] ',
'IDA': 'inferred from direct assay',
'IEP': 'inferred from expression pattern',
'IEA': 'inferred from electronic annotation', ## [to <database:id>]',
'TAS': 'traceable author statement',
'NAS': 'non-traceable author statement',
'ND': 'no biological data available ',
'IC': 'inferred by curator'
}

## returns DAG (dictionary) of maxdepth
def DAGfilterForDepth(dag, node, cdepth):
    filteredDAG = {}
    if cdepth == 0:
        filteredDAG[node] = []

    if cdepth > 0:
        filteredDAG[node] = dag[node]
        for (childNode, rtype) in dag.get(node, []):
            tmpd = DAGfilterForDepth(dag, childNode, cdepth - 1)
            for (key, items) in tmpd.items():
                filteredDAG[key] = items
    return filteredDAG

## creates DAG from GO and two lists of GOIDs:
##   allGOIDs: GOIDS that can be included in DAG (because have genes assigned to them)
##   sigGOIDs: GOIDs that were marked as significant
##             DAG must show only the significant GOIDs and all other GOIDs needed to display a minimal and correct(full) GO
##             correct/full - without any missing paths from bottom of DAG to root
##
## creates a minimal GO displaying all significant GOIDs with genes annotated to them
##
def createGODAGtoDisplay(Lgo, allGOIDs, sigGOIDs):
    DAG = {}
    ## select significant GOIDs (sigGOIDs) and their children that are also in allGOIDs (all terms with genes from cluster set)
    for goid in allGOIDs:
        ## select significant nodes
        sig = goid in sigGOIDs
        ## or those that have a significant parent
        if not(sig):
            nodeparents = Lgo['rGO'].get(goid, []) ## all parents of the node with goid
            ## look if any of the parents significant
            for sigID in sigGOIDs:
                sig = sig or (sigID in nodeparents)
                if sig: break ## found one significant parent, we must include this goid into DAG
        
        if sig:
            ## select only those children that are in allGOIDs
            nodechildren = Lgo['GO'].get(goid, [])
            DAG[goid] = [(childgoid, rtype) for (childgoid, rtype) in nodechildren if childgoid in allGOIDs]

    ## now we must fill in the the rest of the DAG with all the connections root -> children worth displaying (the ones already in the DAG)
    worthDisp = DAG.keys()
    goidsToFillWith = []
    for goid in worthDisp: ## find all parents of those allGOIDs nodes worthDisplaying that are not yet present in DAG
        nodeparents = Lgo['rGO'].get(goid, []) ## all parents
        goidsToFillWith.extend( [parentgoid for parentgoid in nodeparents if (parentgoid in allGOIDs and parentgoid not in worthDisp and parentgoid not in goidsToFillWith)])

    ## fill DAG
    for goid in goidsToFillWith:
        nodechildren = Lgo['GO'].get(goid, [])
        assert( DAG.get(goid, None) == None) ## should not be already present
        DAG[goid] = [(childgoid, rtype) for (childgoid, rtype) in nodechildren if (childgoid in worthDisp + goidsToFillWith)]
    
    ## connect 'root' to all top nodes in DAG
    rootchildren = Lgo['GO'].get('root', [])
    DAG['root'] = [(goid, rtype) for (goid, rtype) in rootchildren if goid in DAG.keys()] ## make root->node connections only for those terms that have significant children

    return DAG

def printNode(dag, goid2term, node, level, allgs, sigs):
    print "\t"*level + str(int(node in allgs)) + " " + str(int(node in sigs)) + " " + str(node) + " " + str(goid2term.get(node, '?'))
    for (childNode, rtype) in dag.get(node, []):
        printNode(dag, goid2term, childNode, level + 1, allgs, sigs)

def printDAG(dag, goid2term, allgs, sigs):
    printNode(dag, goid2term, 'root', 0, allgs, sigs)

def nodeDepth(dag, node, level):
    depth = level
    for (childNode, rype) in dag.get(node, []):
        depth = max(depth, nodeDepth(dag, childNode, level + 1))
    return depth

def DAGdepth(dag):
    return nodeDepth(dag, 'root', 0)

### populate GO with genes
##
## for a given list of genes, annotation and GO it returns a list of all direct and indirect annotations to GO terms
## optional constraint: a list of GO term only for which to return results
##
## used: for cluster and reference frequencies of GO terms
##
def populateGO(geneList, Lann, Lgo, LonlyGOIDs=None, progressBar = None, progressStart = 0.0, progressPart = 100.0):
    if LonlyGOIDs and len(LonlyGOIDs) == 0:
        return {}, {}, {}

    genesGOIDdirect = {}
    genesGOIDindirect = {}
    genesGOIDboth = {}

    pcn = 0.0
    if LonlyGOIDs:
        ## go over only the GOIDs in list
        for daGOID in LonlyGOIDs:
            if progressBar:
                progressBar(progressStart + progressPart * pcn / len(LonlyGOIDs))
                pcn += 1.0
            geneAnn = Lann['GOID2gene'].get(daGOID, None)
            if not(geneAnn): continue

            for (daGene, daNOT, daEvidence, daAspect, daDB_Object_Type) in geneAnn:
                if daAspect <> Lgo['aspect']: continue ## skip annotations different from the loaded GO aspect
                if daGene not in geneList:
                    continue

                ## first include the direct annotation
                tmpl = genesGOIDdirect.get(daGOID, [])
                if (daGene, daEvidence) not in tmpl:
                    genesGOIDdirect[daGOID] = tmpl + [(daGene, daEvidence)] ## update only if GO term in list of terms to return, or list None

                ## update both
                tmpl = genesGOIDboth.get(daGOID, [])
                if (daGene, daEvidence) not in tmpl:
                    genesGOIDboth[daGOID] = tmpl + [(daGene, daEvidence)] ## update only if GO term in list of terms to return, or list None

                ## then all the indirect: by going over all parents of the daGOID, and make indirect annotations to those GO terms
                for GOID in Lgo['rGO'].get(daGOID, []): ##.get(daGOID, []): ## use the reverse GO info to go to all parents
    ##                if GOID == 'root': continue

                    if GOID not in LonlyGOIDs:
                        continue

                    tmpl = genesGOIDindirect.get(GOID, [])
                    if (daGene, daEvidence) not in tmpl:
                        genesGOIDindirect[GOID] = tmpl + [(daGene, daEvidence)]

                    tmpl = genesGOIDboth.get(GOID, [])
                    if (daGene, daEvidence) not in tmpl:
                        genesGOIDboth[GOID] = tmpl + [(daGene, daEvidence)]

    else:
        ## go over all genes and find the apropriate GOIDs
        for gene in geneList: ## go over genes
            if progressBar:
                progressBar(progressStart + progressPart * pcn / len(geneList))
                pcn += 1.0

            geneAnn = Lann['gene2GOID'].get(gene, None)
            if not(geneAnn): continue

            for (daGOID, daNOT, daEvidence, daAspect, daDB_Object_Type) in geneAnn:
                if daAspect <> Lgo['aspect']: continue ## skip annotations different from the loaded GO aspect
    ##            if daNOT <> '': continue ## should we skip those annotations that tell when a gene is not part of a specific GO term?

                ## first include the direct annotation
                tmpl = genesGOIDdirect.get(daGOID, [])
                if (gene, daEvidence) not in tmpl:
                    genesGOIDdirect[daGOID] = tmpl + [(gene, daEvidence)] ## update only if GO term in list of terms to return, or list None

                ## update both
                tmpl = genesGOIDboth.get(daGOID, [])
                if (gene, daEvidence) not in tmpl:
                    genesGOIDboth[daGOID] = tmpl + [(gene, daEvidence)] ## update only if GO term in list of terms to return, or list None

                ## then all the indirect: by going over all parents of the daGOID, and make indirect annotations to those GO terms
                for GOID in Lgo['rGO'].get(daGOID, []): ##.get(daGOID, []): ## use the reverse GO info to go to all parents
    ##                if GOID == 'root': continue

                    tmpl = genesGOIDindirect.get(GOID, [])
                    if (gene, daEvidence) not in tmpl:
                        genesGOIDindirect[GOID] = tmpl + [(gene, daEvidence)]

                    tmpl = genesGOIDboth.get(GOID, [])
                    if (gene, daEvidence) not in tmpl:
                        genesGOIDboth[GOID] = tmpl + [(gene, daEvidence)]

    return genesGOIDdirect, genesGOIDindirect, genesGOIDboth


### main function
## for a given set of genes and reference it recreates the results made with the GO Term Finder Web interface
##
## inputs:
##    clusteSet:    genes for which we want to find significant GO terms
##    referenceSet: list of genes that define our reference, if None or empty list then all genes in the genome are considered for reference
##    evidences:    list of evidence codes to consider from annotation when calculating the reference values for each GO term
##
##    annotation    must always be availabe to this function
##
##
## outputs:
##    sorted list: list of GOterms: (pval, number of genes in term, term GOID), sorted by the increasing p value
##                 used when filtering GO terms by pvalue and number of genes
##
##    GO:          DAG made out of all necessary GO terms needed to display all the genes annotated to significant GO terms
##                                       = they need to be (grand)children of a significant GO term and (grand)parents of a GO term
##                                         that has a direct annotation to a gene in a significant GO term
##
##    GOtermValues:  dictionary of GO terms, giving the calculated values:
##                    term description
##                    cluster frequency
##                    reference/genome frequency
##                    p value
##                    genes annotated to the term directly and indirectly
##                    genes annotated to the term directly
##
##    genesInClusterSet
##    genesInReferenceSet
##    (the last two are used to calculate the relative frequencies)
##
lastFindTermsReference = [0, None, None]
def findTerms(annotation, GO, clusterSet, referenceSet = None, evidences = None, progressBar = None, progressStart = 0.0, progressPart = 100.0):
    global lastFindTermsReference
    if evidences and len(evidences) == 0:
        evidences = None

    ## CLUSTER GENE FREQUENCIES
    ## count the direct and indirect annotations for cluster genes
    n = len(clusterSet) ## number of genes in cluster; used in the calculation of the p value
    clusterGenesGOIDdirect, clusterGenesGOIDindirect, clusterGenesGOID = populateGO(clusterSet, annotation, GO, None, progressBar, progressStart, progressPart / 5.0)
    ## when calculating the p value we use both, the direct and indirect count
    ## but when selecting a node or a subtree we can use the direct and/or indirect


    ## REFERENCE GO TERM (ALL GENES) FREQUENCIES
    ## calculate the reference values only in GO terms that belong to cluster genes only, last parameter in populateGO
    ## if we already made this calculation than reuse it, otherwise do it and remember the results

    ## referenceSet can be None: whole genome
    ##              can be a list of genes that need to be used for reference
    resID = str(id(clusterSet)) + str(id(referenceSet)) + str(id(annotation)) + str(id(GO))
    prevID, prevRefGenesGOID, prevReferenceSet = lastFindTermsReference
    if resID == prevID:
        refGenesGOID = prevRefGenesGOID ## same inputs, we can use old results
        referenceSet = prevReferenceSet
    else:
        ## input changed, new references need to be calculated
        if referenceSet == None:
            referenceSet = annotation['gene2GOID'].keys() ## use all genes in the annotation
        ## calculate frequencies for all GO terms
        ## so we don't have to do it next time, the clusterGenes set changes, because it takes a lot of time anyway
        refGenesGOIDdirect, refGenesGOIDindirect, refGenesGOID = populateGO(referenceSet, annotation, GO, clusterGenesGOID.keys(), progressBar, progressStart + progressPart / 2.0, 4.0 * progressPart / 5.0)
        lastFindTermsReference = [resID, refGenesGOID, referenceSet]

    N = len(referenceSet) ## number of genes in reference; used in the calculation of the p value
    ## the reference set for the whole genome should be calculated at the time of data loading (GO or annotation)
    ## and kept in memory, otherwise it will be too slow and unbearable to the user

    GOtermValues = {}
    sortedGOIDs = []
    if N > 0:
        for (GOID, genesInGOID) in clusterGenesGOID.items(): ## calculate values
            ## count the number of different genes in reference
            lst = refGenesGOID.get(GOID, [])
            genesInRef = []
            for (gene, daEvidence) in lst:
                if (not(evidences) or daEvidence in evidences) and gene not in genesInRef:
                    genesInRef.append( gene)

            ## count the number of different genes in cluster; direct and indirect annotation type
            genesInCluster = []
            for (gene, daEvidence) in genesInGOID:
                if (not(evidences) or daEvidence in evidences) and gene not in genesInCluster:
                    genesInCluster.append( gene)
            ##
            ## count the number of different genes in cluster; direct annotation only
            genesInClusterDirect = []
            genesInGOIDdirect = clusterGenesGOIDdirect.get(GOID, [])
            for (gene, daEvidence) in genesInGOIDdirect:
                if (not(evidences) or daEvidence in evidences) and gene not in genesInClusterDirect:
                    genesInClusterDirect.append( gene)

            G = len(genesInRef) ## reference frequency = all genes in reference for the selected GO term
            p = float(G) / float(N)
            x = len(genesInCluster) ## cluster frequency
            pval = sum([(binom(n, j) * math.pow(p, j) * math.pow(1.0-p, n-j)) for j in range(x, n+1, 1)])
            GOtermValues[GOID] = (GO['term'].get(GOID, '?')[0], x, G, pval, genesInCluster, genesInClusterDirect)
            sortedGOIDs.append( (pval, x, GOID))
        sortedGOIDs.sort()

    if progressBar:
        progressBar(progressStart + progressPart)
    return sortedGOIDs, GOtermValues, clusterSet, referenceSet


if __name__=="__main__":
    ### load the annotation and GO
    ### all of this is done somewhere inside the widget
    annotationFile = r"Annotations\Saccharomyces cerevisiae.annotation"
    GOfile = r"GO\200312-biological_process.go"

    ## read in the annotation
    testAnnotation = cPickle.load(open(annotationFile, 'r'))
    ## annotation is a dictionary: {'GOID2gene': GOID2gene, 'gene2GOID': gene2GOID, 'evidenceTypes': evidenceTypes}

    ## read in the GO info
    testGO = cPickle.load(open(GOfile, 'r'))
    ## GO = {'aspect': aspect, 'term': term, 'relationTypes': termRelTypes, 'GO': GOIDtoGOID[aspect], 'rGO': rGOIDtoGOID[aspect]}
    ###

    print "GO and annotation loaded"

    ## test inputs
    genes = ['YPD1', 'WHI4', 'SHS1', 'GCS1', 'HO', 'YDL228C', 'SSB1', 'PTP1', 'BRE4', 'OST4', 'YDL233W', 'GYP7']
    evidencesToConsider = [] ##['TAS']

    sortedGOIDs, GOtermValues, clusterSet, referenceSet = findTerms(testAnnotation, testGO, genes, None, evidencesToConsider)

    ## display
    alpha = 0.1
    n = len(clusterSet)
    N = len(referenceSet)

    sigGOIDs = [] ## significant GOID to display in GO
    for (p, x, GOID) in sortedGOIDs:
        if p > alpha: break ## end of significant GO terms reached
        if x <= 1:
            continue ## not worth mentioning
        sigGOIDs.append( GOID)
        GOterm, x, G, pval, genesInGOID, genesInGOIDdirect = GOtermValues[GOID]
        print GOID, "\t", GOterm, "\t", len(genesInGOID), "out of", n, "\t", G, "out of", N, "\t", pval, "\t", genesInGOID

    DAG = createGODAGtoDisplay(testGO, GOtermValues.keys(), sigGOIDs)
    printDAG(DAG, testGO['term'], GOtermValues.keys(), sigGOIDs)

    print
    print "filtered DAG"
    fDAG = DAGfilterForDepth(DAG, 'root', 1000)
    print "fDAG:", fDAG
    printDAG(fDAG, testGO['term'], GOtermValues.keys(), sigGOIDs)

