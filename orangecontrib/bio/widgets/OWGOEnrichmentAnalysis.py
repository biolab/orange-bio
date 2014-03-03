"""
<name>GO Browser</name>
<description>Enrichment analysis for Gene Ontology terms.</description>
<contact>Ales Erjavec</contact>
<icon>icons/GOBrowser.svg</icon>
<priority>2020</priority>
"""

from __future__ import absolute_import, with_statement

from collections import defaultdict
from functools import partial
import gc
import sys, os, tarfile, math
from os.path import join as p_join

from Orange.orng import orngServerFiles
from Orange.orng.orngDataCaching import data_hints
from Orange.OrangeWidgets import OWGUI
from Orange.OrangeWidgets.OWWidget import *
from Orange.OrangeWidgets.OWConcurrent import ThreadExecutor

from .. import obiGene, obiGO, obiProb, obiTaxonomy
from .utils.download import EnsureDownloaded

NAME = "GO Browser"
DESCRIPTION = "Enrichment analysis for Gene Ontology terms."
ICON = "icons/GOBrowser.svg"
PRIORITY = 2020

INPUTS = [("Cluster Examples", Orange.data.Table,
           "SetClusterDataset", Single + Default),
          ("Reference Examples", Orange.data.Table,
           "SetReferenceDataset")]

OUTPUTS = [("Selected Examples", Orange.data.Table),
           ("Unselected Examples", Orange.data.Table),
           ("Example With Unknown Genes", Orange.data.Table),
           ("Enrichment Report", Orange.data.Table)]

REPLACES = ["_bioinformatics.widgets.OWGOEnrichmentAnalysis.OWGOEnrichmentAnalysis"]


dataDir = orngServerFiles.localpath("GO")

def listAvailable():
    files = orngServerFiles.listfiles("GO")
    ret = {}
    for file in files:
        tags = orngServerFiles.info("GO", file)["tags"]
        td = dict([tuple(tag.split(":")) for tag in tags if tag.startswith("#") and ":" in tag])
        if "association" in file.lower():
            ret[td.get("#organism", file)] = file
    orgMap = {"352472":"44689"}
    essential = ["gene_association.%s.tar.gz" % obiGO.from_taxid(id) for id in obiTaxonomy.essential_taxids() if obiGO.from_taxid(id)]
    essentialNames = [obiTaxonomy.name(id) for id in obiTaxonomy.essential_taxids() if obiGO.from_taxid(id)]
    ret.update(zip(essentialNames, essential))
    return ret

class _disablegc(object):
    def __enter__(self):
        gc.disable()
    def __exit__(self, *args):
        gc.enable()

def getOrgFileName(org):
    from Orange.orng import orngServerFiles
    files = orngServerFiles.listfiles("go")
    return [f for f in files if org in f].pop()

class TreeNode(object):
    def __init__(self, tuple, children):
        self.tuple = tuple
        self.children = children

class GOTreeWidget(QTreeWidget):
    def contextMenuEvent(self, event):
        QTreeWidget.contextMenuEvent(self, event)
##        print event.x(), event.y()
        term = self.itemAt(event.pos()).term
        self._currMenu = QMenu()
        self._currAction = self._currMenu.addAction("View term on AmiGO website")
##        self.connect(self, SIGNAL("triggered(QAction*)"), partial(self.BrowserAction, term))
        self.connect(self._currAction, SIGNAL("triggered()"), lambda :self.BrowserAction(term))
        self._currMenu.popup(event.globalPos())

    def BrowserAction(self, term):
        import webbrowser
        if isinstance(term, obiGO.Term):
            term = term.id
        webbrowser.open("http://amigo.geneontology.org/cgi-bin/amigo/term-details.cgi?term="+term)
        
    def paintEvent(self, event):
        QTreeWidget.paintEvent(self, event)
        if getattr(self, "_userMessage", None):
            painter = QPainter(self.viewport())
            font = QFont(self.font())
            font.setPointSize(15)
            painter.setFont(font)
            painter.drawText(self.viewport().geometry(), Qt.AlignCenter, self._userMessage)
            painter.end()


class OWGOEnrichmentAnalysis(OWWidget):
    settingsList = ["annotationIndex", "useReferenceDataset", "aspectIndex",
                    "geneAttrIndex", "geneMatcherSettings",
                    "filterByNumOfInstances", "minNumOfInstances",
                    "filterByPValue", "maxPValue", "selectionDirectAnnotation",
                    "filterByPValue_nofdr", "maxPValue_nofdr",
                    "selectionDisjoint", "selectionType",
                    "selectionAddTermAsClass", "useAttrNames", "probFunc"
                    ]

    contextHandlers = {"": DomainContextHandler(
                                "",
                                ["geneAttrIndex", "useAttrNames",
                                 "annotationIndex", "geneMatcherSettings"],
                                matchValues=1)
                       }

    def __init__(self, parent=None, signalManager=None, name="GO Browser"):
        OWWidget.__init__(self, parent, signalManager, name)
        self.inputs = [("Cluster Examples", ExampleTable,
                        self.SetClusterDataset, Default),
                       ("Reference Examples", ExampleTable,
                        self.SetReferenceDataset, Single + NonDefault)]

        self.outputs = [("Selected Examples", ExampleTable, Default),
                        ("Unselected Examples", ExampleTable, Default),
                        ("Example With Unknown Genes", ExampleTable, Default),
                        ("Enrichment Report", ExampleTable)]

        self.annotationIndex = 0
        self.autoFindBestOrg = False
        self.useReferenceDataset = 0
        self.aspectIndex = 0
        self.geneAttrIndex = 0
        self.useAttrNames = False
        self.geneMatcherSettings = [True, False, False, False]
        self.filterByNumOfInstances = False
        self.minNumOfInstances = 1
        self.filterByPValue = True
        self.maxPValue = 0.1
        self.filterByPValue_nofdr = False
        self.maxPValue_nofdr = 0.1
        self.probFunc = 0
        self.selectionDirectAnnotation = 0
        self.selectionDisjoint = 0
        self.selectionAddTermAsClass = 0
        self.selectionChanging = 0
        
        # check usage of all evidences
        for etype in obiGO.evidenceTypesOrdered:
            varName = "useEvidence" + etype
            if varName not in self.settingsList: 
                self.settingsList.append(varName)
            code = compile("self.%s = True" % (varName), ".", "single")
            exec(code)
        self.annotationCodes = []
        
        self.loadSettings()
        
        #############
        ##GUI
        #############
        self.tabs = OWGUI.tabWidget(self.controlArea)
        ##Input tab
        self.inputTab = OWGUI.createTabPage(self.tabs, "Input")
        box = OWGUI.widgetBox(self.inputTab, "Info")
        self.infoLabel = OWGUI.widgetLabel(box, "No data on input\n")
        
        OWGUI.button(box, self, "Ontology/Annotation Info", callback=self.ShowInfo,
                     tooltip="Show information on loaded ontology and annotations",
                     debuggingEnabled=0)
        box = OWGUI.widgetBox(self.inputTab, "Organism", addSpace=True)
        self.annotationComboBox = OWGUI.comboBox(box, self, "annotationIndex",
                            items = self.annotationCodes, callback=self.Update,
                            tooltip="Select organism", debuggingEnabled=0)

        self.geneAttrIndexCombo = OWGUI.comboBox(self.inputTab, self, "geneAttrIndex",
                            box="Gene names", callback=self.Update,
                            tooltip="Use this attribute to extract gene names from input data")
        OWGUI.checkBox(self.geneAttrIndexCombo.box, self, "useAttrNames", "Use data attributes names",
                       disables=[(-1, self.geneAttrIndexCombo)], callback=self.Update, 
                       tooltip="Use attribute names for gene names")
        OWGUI.button(self.geneAttrIndexCombo.box, self, "Gene matcher settings", 
                     callback=self.UpdateGeneMatcher, 
                     tooltip="Open gene matching settings dialog", 
                     debuggingEnabled=0)
        
        self.referenceRadioBox = OWGUI.radioButtonsInBox(self.inputTab, self, "useReferenceDataset", 
                                                         ["Entire genome", "Reference set (input)"],
                                                         tooltips=["Use entire genome for reference",
                                                                   "Use genes from Referece Examples input signal as reference"],
                                                         box="Reference", callback=self.Update)
        self.referenceRadioBox.buttons[1].setDisabled(True)
        OWGUI.radioButtonsInBox(self.inputTab, self, "aspectIndex", ["Biological process",
                                                                     "Cellular component",
                                                                     "Molecular function"], 
                                box="Aspect", callback=self.Update)
        
        self.geneAttrIndexCombo.setDisabled(bool(self.useAttrNames))
        
        ##Filter tab
        self.filterTab = OWGUI.createTabPage(self.tabs, "Filter")
        box = OWGUI.widgetBox(self.filterTab, "Filter GO Term Nodes", addSpace=True)
        OWGUI.checkBox(box, self, "filterByNumOfInstances", "Genes",
                       callback=self.FilterAndDisplayGraph, 
                       tooltip="Filter by number of input genes mapped to a term")
        OWGUI.spin(OWGUI.indentedBox(box), self, 'minNumOfInstances', 1, 100, 
                   step=1, label='#:', labelWidth=15, 
                   callback=self.FilterAndDisplayGraph, 
                   callbackOnReturn=True, 
                   tooltip="Min. number of input genes mapped to a term")
        
        OWGUI.checkBox(box, self, "filterByPValue_nofdr", "p-value",
                       callback=self.FilterAndDisplayGraph, 
                       tooltip="Filter by term p-value")
        OWGUI.doubleSpin(OWGUI.indentedBox(box), self, 'maxPValue_nofdr', 1e-8, 1, 
                         step=1e-8,  label='p:', labelWidth=15, 
                         callback=self.FilterAndDisplayGraph, 
                         callbackOnReturn=True, 
                         tooltip="Max term p-value")

        #use filterByPValue for FDR, as it was the default in prior versions
        OWGUI.checkBox(box, self, "filterByPValue", "FDR",
                       callback=self.FilterAndDisplayGraph, 
                       tooltip="Filter by term FDR")
        OWGUI.doubleSpin(OWGUI.indentedBox(box), self, 'maxPValue', 1e-8, 1, 
                         step=1e-8,  label='p:', labelWidth=15, 
                         callback=self.FilterAndDisplayGraph, 
                         callbackOnReturn=True, 
                         tooltip="Max term p-value")

        box = OWGUI.widgetBox(box, "Significance test")

        OWGUI.radioButtonsInBox(box, self, "probFunc", ["Binomial", "Hypergeometric"], 
                                tooltips=["Use binomial distribution test", 
                                          "Use hypergeometric distribution test"], 
                                callback=self.Update)
        box = OWGUI.widgetBox(self.filterTab, "Evidence codes in annotation", 
                              addSpace=True)
        self.evidenceCheckBoxDict = {}
        for etype in obiGO.evidenceTypesOrdered:
            self.evidenceCheckBoxDict[etype] = OWGUI.checkBox(box, self, "useEvidence"+etype, etype,
                                            callback=self.Update, tooltip=obiGO.evidenceTypes[etype])
        
        ##Select tab
        self.selectTab = OWGUI.createTabPage(self.tabs, "Select")
        box = OWGUI.radioButtonsInBox(self.selectTab, self, "selectionDirectAnnotation", 
                                      ["Directly or Indirectly", "Directly"], 
                                      box="Annotated genes", 
                                      callback=self.ExampleSelection)
        
        box = OWGUI.widgetBox(self.selectTab, "Output", addSpace=True)
        OWGUI.radioButtonsInBox(box, self, "selectionDisjoint", 
                                btnLabels=["All selected genes", 
                                           "Term-specific genes", 
                                           "Common term genes"], 
                                tooltips=["Outputs genes annotated to all selected GO terms", 
                                          "Outputs genes that appear in only one of selected GO terms", 
                                          "Outputs genes common to all selected GO terms"], 
                                callback=[self.ExampleSelection,
                                          self.UpdateAddClassButton])
        self.addClassCB = OWGUI.checkBox(box, self, "selectionAddTermAsClass",
                                         "Add GO Term as class", 
                                         callback=self.ExampleSelection)

        # ListView for DAG, and table for significant GOIDs
        self.DAGcolumns = ['GO term', 'Cluster', 'Reference', 'p-value', 'FDR', 'Genes', 'Enrichment']
        
        self.splitter = QSplitter(Qt.Vertical, self.mainArea)
        self.mainArea.layout().addWidget(self.splitter)

        # list view
        self.listView = GOTreeWidget(self.splitter)
        self.listView.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.listView.setAllColumnsShowFocus(1)
        self.listView.setColumnCount(len(self.DAGcolumns))
        self.listView.setHeaderLabels(self.DAGcolumns)
        
        self.listView.header().setClickable(True)
        self.listView.header().setSortIndicatorShown(True)
        self.listView.setSortingEnabled(True)
        self.listView.setItemDelegateForColumn(6, EnrichmentColumnItemDelegate(self))
        self.listView.setRootIsDecorated(True)

        
        self.connect(self.listView, SIGNAL("itemSelectionChanged()"), self.ViewSelectionChanged)
        
        # table of significant GO terms
        self.sigTerms = QTreeWidget(self.splitter)
        self.sigTerms.setColumnCount(len(self.DAGcolumns))
        self.sigTerms.setHeaderLabels(self.DAGcolumns)
        self.sigTerms.setSortingEnabled(True)
        self.sigTerms.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.sigTerms.setItemDelegateForColumn(6, EnrichmentColumnItemDelegate(self))
        
        self.connect(self.sigTerms, SIGNAL("itemSelectionChanged()"), self.TableSelectionChanged)
        self.splitter.show()

        self.sigTableTermsSorted = []
        self.graph = {}
        
        self.loadedAnnotationCode = "---"
        
        self.inputTab.layout().addStretch(1)
        self.filterTab.layout().addStretch(1)
        self.selectTab.layout().addStretch(1)
        
        self.resize(1000, 800)

        self.clusterDataset = None
        self.referenceDataset = None
        self.ontology = None
        self.annotations = None
        self.treeStructRootKey = None
        self.probFunctions = [obiProb.Binomial(), obiProb.Hypergeometric()]
        self.selectedTerms = []

        self.connect(self, SIGNAL("widgetStateChanged(QString, int, QString)"), self.onStateChanged)

        self.setBlocking(True)
        self._executor = ThreadExecutor()
        self._init = EnsureDownloaded(
            [("Taxonomy", "ncbi_taxonomy.tar.gz"),
             ("GO", "taxonomy.pickle")]
        )
        self._init.finished.connect(self.UpdateOrganismComboBox)
        self._executor.submit(self._init)

    def UpdateOrganismComboBox(self):
        try:
            self.annotationFiles = listAvailable()
            self.annotationCodes = sorted(self.annotationFiles.keys())
            self.annotationComboBox.clear()
            self.annotationComboBox.addItems(self.annotationCodes)
            self.annotationComboBox.setCurrentIndex(self.annotationIndex)
        finally:
            self.setBlocking(False)

    def UpdateGeneMatcher(self):
        dialog = GeneMatcherDialog(self, defaults=self.geneMatcherSettings, modal=True)
        if dialog.exec_():
            self.geneMatcherSettings = [getattr(dialog, item[0]) for item in dialog.items]
            if self.annotations:
                self.SetGeneMatcher()
                if self.clusterDataset:
                    self.Update()
                
    def Update(self):
        if self.clusterDataset:
            pb = OWGUI.ProgressBar(self, 100)
            self.Load(pb=pb)
            self.FilterUnknownGenes()
            graph = self.Enrichment(pb=pb)
            self.SetGraph(graph)

    def UpdateGOAndAnnotation(self, tags=[]):
        from .OWUpdateGenomicsDatabases import OWUpdateGenomicsDatabases
        w = OWUpdateGenomicsDatabases(parent = self, searchString=" ".join(tags))
        w.setModal(True)
        w.show()
        self.UpdateAnnotationComboBox()

    def UpdateAnnotationComboBox(self):
        if self.annotationCodes:
            curr = self.annotationCodes[min(self.annotationIndex, len(self.annotationCodes)-1)]
        else:
            curr = None
        self.annotationFiles = listAvailable()
        self.annotationCodes = self.annotationFiles.keys()
        index = curr and self.annotationCodes.index(curr) or 0
        self.annotationComboBox.clear()
        self.annotationComboBox.addItems(self.annotationCodes)
        self.annotationComboBox.setCurrentIndex(index)
        if not self.annotationCodes:
            self.error(0, "No downloaded annotations!!\nClick the update button and update annotationa for at least one organism!")
        else:
            self.error(0)

    def SetGenesComboBox(self):
        self.candidateGeneAttrs = self.clusterDataset.domain.variables + self.clusterDataset.domain.getmetas().values()
        self.candidateGeneAttrs = filter(lambda v: v.varType in [orange.VarTypes.String,
                                                                 orange.VarTypes.Other,
                                                                 orange.VarTypes.Discrete], 
                                         self.candidateGeneAttrs)
        self.geneAttrIndexCombo.clear()
        self.geneAttrIndexCombo.addItems([a.name for a in  self.candidateGeneAttrs])

    def FindBestGeneAttrAndOrganism(self):
        if self.autoFindBestOrg:  
            organismGenes = dict([(o,set(go.getCachedGeneNames(o))) for o in self.annotationCodes])
        else:
            currCode = self.annotationCodes[min(self.annotationIndex, len(self.annotationCodes)-1)]
            filename = p_join(dataDir, self.annotationFiles[currCode])
            try:
                f = tarfile.open(filename)
                info = [info for info in f.getmembers() if info.name.startswith("gene_names")].pop()
                geneNames = cPickle.loads(f.extractfile(info).read().replace("\r\n", "\n"))
            except Exception, ex:
                geneNames = cPickle.loads(open(p_join(filename, "gene_names.pickle")).read().replace("\r\n", "\n"))
            organismGenes = {currCode: set(geneNames)}
        candidateGeneAttrs = self.clusterDataset.domain.attributes + self.clusterDataset.domain.getmetas().values()
        candidateGeneAttrs = filter(lambda v: v.varType in [orange.VarTypes.String, 
                                                            orange.VarTypes.Other, 
                                                            orange.VarTypes.Discrete], 
                                    candidateGeneAttrs)
        attrNames = [v.name for v in self.clusterDataset.domain.variables]
        cn = {}
        for attr in candidateGeneAttrs:
            vals = [str(e[attr]) for e in self.clusterDataset]
            if any("," in val for val in vals):
                vals = reduce(list.__add__, (val.split(",") for val in vals))
            for organism, s in organismGenes.items():
                l = filter(lambda a: a in s, vals)
                cn[(attr,organism)] = len(set(l))
        for organism, s in organismGenes.items():
            l = filter(lambda a: a in s, attrNames)
            cn[("_var_names_", organism)] = len(set(l))
            
        cn = cn.items()
        cn.sort(lambda a,b:-cmp(a[1],b[1]))
        ((bestAttr, organism), count) = cn[0]
        if bestAttr=="_var_names_" and count<=len(attrNames)/10.0 or \
           bestAttr!="_var_names_" and count<=len(self.clusterDataset)/10.0:
            return
        
        self.annotationIndex = self.annotationCodes.index(organism)
        if bestAttr=="_var_names_":
            self.useAttrNames = True
            self.geneAttrIndex = 0
        else:
            self.useAttrNames = False
            self.geneAttrIndex = candidateGeneAttrs.index(bestAttr)

    def SetClusterDataset(self, data=None):
        if not self.annotationCodes:
            QTimer.singleShot(200, lambda: self.SetClusterDataset(data))
            return
        self.closeContext()
        self.clusterDataset = data
        self.infoLabel.setText("\n")
        if data:
            self.SetGenesComboBox()
            try:
                taxid = data_hints.get_hint(data, "taxid", "")
                code = obiGO.from_taxid(taxid)
                filename = "gene_association.%s.tar.gz" % code
                if filename in self.annotationFiles.values():
                    self.annotationIndex = \
                            [i for i, name in enumerate(self.annotationCodes) \
                             if self.annotationFiles[name] == filename].pop()
            except Exception:
                pass
            self.useAttrNames = data_hints.get_hint(data, "genesinrows",
                                                    self.useAttrNames)
            self.openContext("", data)
            self.Update()
        else:
            self.infoLabel.setText("No data on input\n")
            self.warning(0)
            self.warning(1)
            self.openContext("", None)
            self.ClearGraph()
            self.send("Selected Examples", None)
            self.send("Unselected Examples", None)
            self.send("Example With Unknown Genes", None)
            self.send("Enrichment Report", None)

    def SetReferenceDataset(self, data=None):
        self.referenceDataset = data
        self.referenceRadioBox.buttons[1].setDisabled(not bool(data))
        self.referenceRadioBox.buttons[1].setText("Reference set")
        if self.clusterDataset and self.useReferenceDataset:
            self.useReferenceDataset = 0 if not data else 1
            graph = self.Enrichment()
            self.SetGraph(graph)
        elif self.clusterDataset:
            self.UpdateReferenceSetButton()

    def UpdateReferenceSetButton(self):
        allgenes, refgenes = None, None
        if self.referenceDataset:
            try:
                allgenes = self.GenesFromExampleTable(self.referenceDataset)
            except Exception:
                allgenes = []
            refgenes, unknown = self.FilterAnnotatedGenes(allgenes)
        self.referenceRadioBox.buttons[1].setDisabled(not bool(allgenes))
        self.referenceRadioBox.buttons[1].setText("Reference set " + ("(%i genes, %i matched)" % (len(allgenes), len(refgenes)) if allgenes and refgenes else ""))

    def GenesFromExampleTable(self, data):
        if self.useAttrNames:
            genes = [v.name for v in data.domain.variables]
        else:
            attr = self.candidateGeneAttrs[min(self.geneAttrIndex, len(self.candidateGeneAttrs) - 1)]
            genes = [str(ex[attr]) for ex in data if not ex[attr].isSpecial()]
            if any("," in gene for gene in genes):
                self.information(0, "Separators detected in gene names. Assuming multiple genes per example.")
                genes = reduce(list.__add__, (genes.split(",") for genes in genes))
        return genes
        
    def FilterAnnotatedGenes(self, genes):
        matchedgenes = self.annotations.GetGeneNamesTranslator(genes).values()
        return matchedgenes, [gene for gene in genes if gene not in matchedgenes]
        
    def FilterUnknownGenes(self):
        if not self.useAttrNames and self.candidateGeneAttrs:
            geneAttr = self.candidateGeneAttrs[min(self.geneAttrIndex, len(self.candidateGeneAttrs)-1)]
            examples = []
            for ex in self.clusterDataset:
                if not any(self.annotations.genematcher.match(n.strip()) for n in str(ex[geneAttr]).split(",")):
                    examples.append(ex)

            self.send("Example With Unknown Genes", examples and orange.ExampleTable(examples) or None)
        else:
            self.send("Example With Unknown Genes", None)

    def Load(self, pb=None):
        go_files, tax_files = orngServerFiles.listfiles("GO"), orngServerFiles.listfiles("Taxonomy")
        calls = []
        pb, finish = (OWGUI.ProgressBar(self, 0), True) if pb is None else (pb, False)
        count = 0
        if not tax_files:
            calls.append(("Taxonomy", "ncbi_taxnomy.tar.gz"))
            count += 1
        org = self.annotationCodes[min(self.annotationIndex, len(self.annotationCodes)-1)]
        if org != self.loadedAnnotationCode:
            count += 1
            if self.annotationFiles[org] not in go_files:
                calls.append(("GO", self.annotationFiles[org]))
                count += 1
                
        if "gene_ontology_edit.obo.tar.gz" not in go_files:
            calls.append(("GO", "gene_ontology_edit.obo.tar.gz"))
            count += 1
        if not self.ontology:
            count += 1
        pb.iter += count*100
        
        for i, args in enumerate(calls):
            orngServerFiles.localpath_download(*args, **dict(callback=pb.advance))
            
        i = len(calls)
        if not self.ontology:
            self.ontology = obiGO.Ontology(progressCallback=lambda value: pb.advance())
            i+=1
        if org != self.loadedAnnotationCode:
            self.annotations = None
            gc.collect() # Force run garbage collection.
            code = self.annotationFiles[org].split(".")[-3]
            self.annotations = obiGO.Annotations(code, genematcher=obiGene.GMDirect(), progressCallback=lambda value: pb.advance())
            i+=1
            self.loadedAnnotationCode = org
            count = defaultdict(int)
            geneSets = defaultdict(set)

            for anno in self.annotations.annotations:
                count[anno.evidence]+=1
                geneSets[anno.evidence].add(anno.geneName)
            for etype in obiGO.evidenceTypesOrdered:
                self.evidenceCheckBoxDict[etype].setEnabled(bool(count[etype]))
                self.evidenceCheckBoxDict[etype].setText(etype+": %i annots(%i genes)" % (count[etype], len(geneSets[etype])))
        if finish:
            pb.finish()
            
    def SetGeneMatcher(self):
        if self.annotations:
            taxid = self.annotations.taxid
            matchers = []
            for matcher, use in zip([obiGene.GMGO, obiGene.GMKEGG, obiGene.GMNCBI, obiGene.GMAffy], self.geneMatcherSettings):
                if use:
                    try:
                        if taxid == "352472":
                            matchers.extend([matcher(taxid), obiGene.GMDicty(),
                                            [matcher(taxid), obiGene.GMDicty()]])
                            # The reason machers are duplicated is that we want `matcher` or `GMDicty` to
                            # match genes by them self if possible. Only use the joint matcher if they fail.   
                        else:
                            matchers.append(matcher(taxid))
                    except Exception, ex:
                        print ex
            self.annotations.genematcher = obiGene.matcher(matchers)
            self.annotations.genematcher.set_targets(self.annotations.geneNames)
            
    def Enrichment(self, pb=None):
        pb = OWGUI.ProgressBar(self, 100) if pb is None else pb
        if not self.annotations.ontology:
            self.annotations.ontology = self.ontology
            
        if isinstance(self.annotations.genematcher, obiGene.GMDirect):
            self.SetGeneMatcher()
        self.error(1)
        self.warning([0, 1])
        try:    
            if self.useAttrNames:
                clusterGenes = [v.name for v in self.clusterDataset.domain.variables]
                self.information(0)
            else:
                geneAttr = self.candidateGeneAttrs[min(self.geneAttrIndex, len(self.candidateGeneAttrs)-1)]
                clusterGenes = [str(ex[geneAttr]) for ex in self.clusterDataset if not ex[geneAttr].isSpecial()]
                if any("," in gene for gene in clusterGenes):
                    self.information(0, "Separators detected in cluster gene names. Assuming multiple genes per example.")
                    clusterGenes = reduce(list.__add__, (genes.split(",") for genes in clusterGenes))
                else:
                    self.information(0)
        except Exception, ex:
            self.error(1, "Failed to extract gene names from input dataset! %s" % str(ex))
            return {}
        genesCount = len(clusterGenes)
        genesSetCount = len(set(clusterGenes))
        
        self.clusterGenes = clusterGenes = self.annotations.GetGeneNamesTranslator(clusterGenes).values()
        
#        self.clusterGenes = clusterGenes = filter(lambda g: g in self.annotations.aliasMapper or g in self.annotations.additionalAliases, clusterGenes)
        self.infoLabel.setText("%i unique genes on input\n%i (%.1f%%) genes with known annotations" % (genesSetCount, len(clusterGenes), 100.0*len(clusterGenes)/genesSetCount if genesSetCount else 0.0))
        
        referenceGenes = None
        if self.referenceDataset:
            try:
                if self.useAttrNames:
                    referenceGenes = [v.name for v in self.referenceDataset.domain.variables]
                    self.information(1)
                else:
                    referenceGenes = [str(ex[geneAttr]) for ex in self.referenceDataset if not ex[geneAttr].isSpecial()]
                    if any("," in gene for gene in clusterGenes):
                        self.information(1, "Separators detected in reference gene names. Assuming multiple genes per example.")
                        referenceGenes = reduce(list.__add__, (genes.split(",") for genes in referenceGenes))
                    else:
                        self.information(1)

                refc = len(referenceGenes)
#                referenceGenes = filter(lambda g: g in self.annotations.aliasMapper or g in self.annotations.additionalAliases, referenceGenes)
                referenceGenes = self.annotations.GetGeneNamesTranslator(referenceGenes).values()
                self.referenceRadioBox.buttons[1].setText("Reference set (%i genes, %i matched)" % (refc, len(referenceGenes)))
                self.referenceRadioBox.buttons[1].setDisabled(False)
                self.information(2)
            except Exception, er:
                if not self.referenceDataset:
                    self.information(2, "Unable to extract gene names from reference dataset. Using entire genome for reference")
                else:
                    self.referenceRadioBox.buttons[1].setText("Reference set")
                    self.referenceRadioBox.buttons[1].setDisabled(True)
                referenceGenes = self.annotations.geneNames
                self.useReferenceDataset = 0
        else:
            self.useReferenceDataset = 0
        if not self.useReferenceDataset:
            self.information(2)
            self.information(1)
            referenceGenes = self.annotations.geneNames
        self.referenceGenes = referenceGenes
        evidences = []
        for etype in obiGO.evidenceTypesOrdered:
            if getattr(self, "useEvidence"+etype):
                evidences.append(etype)
        aspect = ["P", "C", "F"][self.aspectIndex]
        
        if clusterGenes:
            self.terms = terms = self.annotations.GetEnrichedTerms(clusterGenes, referenceGenes, evidences, aspect=aspect,
                                                                   prob=self.probFunctions[self.probFunc], useFDR=False,
                                                                   progressCallback=lambda value:pb.advance() )
            ids = []
            pvals = []
            for i,d in self.terms.items():
                ids.append(i)
                pvals.append(d[1])
            for i,fdr in zip(ids, obiProb.FDR(pvals)): #save FDR as the last part of the tuple
                terms[i] = tuple(list(terms[i]) + [ fdr ])

        else:
            self.terms = terms = {}
        if not self.terms:
            self.warning(0, "No enriched terms found.")
        else:
            self.warning(0)
            
        pb.finish()
        self.treeStructDict = {}
        ids = self.terms.keys()
        
        self.treeStructRootKey = None
        
        parents = {}
        for id in ids:
            parents[id] = set([term for typeId, term in self.ontology[id].related])
            
        children = {}
        for term in self.terms:
            children[term] = set([id for id in ids if term in parents[id]])
            
        for term in self.terms:
            self.treeStructDict[term] = TreeNode(self.terms[term], children[term])
            if not self.ontology[term].related and not getattr(self.ontology[term], "is_obsolete", False):
                self.treeStructRootKey = term
        return terms
        
    def FilterGraph(self, graph):
        if self.filterByPValue_nofdr:
            graph = obiGO.filterByPValue(graph, self.maxPValue_nofdr)
        if self.filterByPValue: #FDR
            graph = dict(filter(lambda (k, e): e[3] <= self.maxPValue, graph.items()))
        if self.filterByNumOfInstances:
            graph = dict(filter(lambda (id,(genes, p, rc, fdr)):len(genes)>=self.minNumOfInstances, graph.items()))
        return graph

    def FilterAndDisplayGraph(self):
        if self.clusterDataset:
            self.graph = self.FilterGraph(self.originalGraph)
            if self.originalGraph and not self.graph:
                self.warning(1, "All found terms were filtered out.")
            else:
                self.warning(1)
            self.ClearGraph()
            self.DisplayGraph()

    def SetGraph(self, graph=None):
        self.originalGraph = graph
        if graph:
            self.FilterAndDisplayGraph()
        else:
            self.graph = {}
            self.ClearGraph()

    def ClearGraph(self):
        self.listView.clear()
        self.listViewItems=[]
        self.sigTerms.clear()

    def DisplayGraph(self):
        fromParentDict = {}
        self.termListViewItemDict = {}
        self.listViewItems = []
        enrichment = lambda t:float(len(t[0])) / t[2] * (float(len(self.referenceGenes))/len(self.clusterGenes))
        maxFoldEnrichment = max([enrichment(term) for term in self.graph.values()] or [1])


        def addNode(term, parent, parentDisplayNode):
            if (parent, term) in fromParentDict:
                return
            if term in self.graph:
                displayNode = GOTreeWidgetItem(self.ontology[term], self.graph[term], len(self.clusterGenes), len(self.referenceGenes), maxFoldEnrichment, parentDisplayNode)
                displayNode.goId = term
                self.listViewItems.append(displayNode)
                if term in self.termListViewItemDict:
                    self.termListViewItemDict[term].append(displayNode)
                else:
                    self.termListViewItemDict[term] = [displayNode]
                fromParentDict[(parent, term)] = True
                parent = term
            else:
                displayNode = parentDisplayNode

            for c in self.treeStructDict[term].children:
                addNode(c, parent, displayNode)

        if self.treeStructDict:
            addNode(self.treeStructRootKey, None, self.listView)

        terms = self.graph.items()
        terms = sorted(terms, key=lambda item: item[1][1])
        self.sigTableTermsSorted = [t[0] for t in terms]

        self.sigTerms.clear()
        for i, (t_id, (genes, p_value, refCount, fdr)) in enumerate(terms):
            item = GOTreeWidgetItem(self.ontology[t_id],
                                    (genes, p_value, refCount, fdr),
                                    len(self.clusterGenes),
                                    len(self.referenceGenes),
                                    maxFoldEnrichment,
                                    self.sigTerms)
            item.goId = t_id

        self.listView.expandAll()
        for i in range(5):
            self.listView.resizeColumnToContents(i)
            self.sigTerms.resizeColumnToContents(i)
        self.sigTerms.resizeColumnToContents(6)
        width = min(self.listView.columnWidth(0), 350)
        self.listView.setColumnWidth(0, width)
        self.sigTerms.setColumnWidth(0, width)

        # Create and send the enrichemnt report table.
        termsDomain = orange.Domain(
            [orange.StringVariable("GO Term Id"),
             orange.StringVariable("GO Term Name"),
             orange.FloatVariable("Cluster Frequency"),
             orange.FloatVariable("Reference Frequency"),
             orange.FloatVariable("p-value"),
             orange.FloatVariable("FDR"),
             orange.FloatVariable("Enrichment"),
             orange.StringVariable("Genes")
             ], None)

        terms = [[t_id,
                  self.ontology[t_id].name,
                  float(len(genes)) / len(self.clusterGenes),
                  float(r_count) / len(self.referenceGenes),
                  p_value,
                  fdr,
                  float(len(genes)) / len(self.clusterGenes) * \
                  float(len(self.referenceGenes)) / r_count,
                  ",".join(genes)
                  ]
                 for t_id, (genes, p_value, r_count, fdr) in terms]

        if terms:
            termsTable = orange.ExampleTable(termsDomain, terms)
        else:
            termsTable = orange.ExampleTable(termsDomain)
        self.send("Enrichment Report", termsTable)

    def ViewSelectionChanged(self):
        if self.selectionChanging:
            return

        self.selectionChanging = 1
        self.selectedTerms = []
        selected = self.listView.selectedItems()
        self.selectedTerms = list(set([lvi.term.id for lvi in selected]))
        self.ExampleSelection()
        self.selectionChanging = 0

    def TableSelectionChanged(self):
        if self.selectionChanging:
            return
        
        self.selectionChanging = 1
        self.selectedTerms = []
        selectedIds = set([self.sigTerms.itemFromIndex(index).goId for index in self.sigTerms.selectedIndexes()])
        
        for i in range(self.sigTerms.topLevelItemCount()):
            item = self.sigTerms.topLevelItem(i)
            selected = item.goId in selectedIds
            term = item.goId
            
            if selected:
                self.selectedTerms.append(term)
                
            for lvi in self.termListViewItemDict[term]:
                try:
                    lvi.setSelected(selected)
                    if selected: lvi.setExpanded(True)
                except RuntimeError:    ##Underlying C/C++ object deleted (why??)
                    pass
                
        self.ExampleSelection()
        self.selectionChanging = 0
            
    
    def UpdateAddClassButton(self):
        self.addClassCB.setEnabled(self.selectionDisjoint == 1)
        
    def ExampleSelection(self):
        selectedExamples = []
        unselectedExamples = []
        selectedGenes = []
        
        #change by Marko. don't do anything if there is no3 dataset 
        if not self.clusterDataset:
            return
        
        selectedGenes = reduce(set.union, [v[0] for id, v in self.graph.items() if id in self.selectedTerms], set())
        evidences = []
        for etype in obiGO.evidenceTypesOrdered:
            if getattr(self, "useEvidence"+etype):
                evidences.append(etype)
        allTerms = self.annotations.GetAnnotatedTerms(selectedGenes, 
                          directAnnotationOnly=self.selectionDirectAnnotation, 
                          evidenceCodes=evidences)
            
        if self.selectionDisjoint:
            count = defaultdict(int)
            for term in self.selectedTerms:
                for g in allTerms.get(term, []):
                    count[g]+=1
            ccount = 1 if self.selectionDisjoint==1 else len(self.selectedTerms)
            selectedGenes = [gene for gene, c in count.items() if c==ccount and gene in selectedGenes]
        else:
            selectedGenes = reduce(set.union, [allTerms.get(term, []) for term in self.selectedTerms], set())

        if self.useAttrNames:
            vars = [self.clusterDataset.domain[gene] for gene in set(selectedGenes)]
            newDomain = orange.Domain(vars, self.clusterDataset.domain.classVar)
            newdata = orange.ExampleTable(newDomain, self.clusterDataset)
            self.send("Selected Examples", newdata)
            self.send("Unselected Examples", None)
        elif self.candidateGeneAttrs:
            geneAttr = self.candidateGeneAttrs[min(self.geneAttrIndex, len(self.candidateGeneAttrs)-1)]
            if self.selectionDisjoint == 1:
                goVar = orange.EnumVariable("GO Term", values=list(self.selectedTerms))
                newDomain = orange.Domain(self.clusterDataset.domain.variables, goVar)
                newDomain.addmetas(self.clusterDataset.domain.getmetas())
#            else:
#                goVar = orange.StringVariable("GO Terms")
#                newDomain = orange.Domain(self.clusterDataset.domain)
#                newDomain.addmeta(orange.newmetaid(), goVar)
            
            
            for ex in self.clusterDataset:
                if not ex[geneAttr].isSpecial() and any(gene in selectedGenes for gene in str(ex[geneAttr]).split(",")):
                    if self.selectionDisjoint == 1 and self.selectionAddTermAsClass:
                        terms = filter(lambda term: any(gene in self.graph[term][0] for gene in str(ex[geneAttr]).split(",")) , self.selectedTerms)
                        term = sorted(terms)[0]
                        ex =  orange.Example(newDomain, ex)
                        ex[goVar] = goVar(term)
#                        ex.setclass(newClass(term))
                    selectedExamples.append(ex)
                else:
                    unselectedExamples.append(ex)
                    
            if selectedExamples:
                selectedExamples = orange.ExampleTable(selectedExamples)
            else:
                selectedExamples = None
                
            if unselectedExamples:
                unselectedExamples = orange.ExampleTable(unselectedExamples)
            else:
                unselectedExamples = None

            self.send("Selected Examples", selectedExamples)
            self.send("Unselected Examples", unselectedExamples)

    def ShowInfo(self):
        dialog = QDialog(self)
        dialog.setModal(False)
        dialog.setLayout(QVBoxLayout())
        label = QLabel(dialog)
        label.setText("Ontology:\n"+self.ontology.header if self.ontology else "Ontology not loaded!")
        dialog.layout().addWidget(label)

        label = QLabel(dialog)
        label.setText("Annotations:\n"+self.annotations.header.replace("!", "") if self.annotations else "Annotations not loaded!")
        dialog.layout().addWidget(label)
        dialog.show()
        
    def sendReport(self):
        self.reportSettings("Settings", [("Organism", self.annotationCodes[min(self.annotationIndex, len(self.annotationCodes) - 1)]),
                                         ("Significance test", ("Binomial" if self.probFunc == 0 else "Hypergeometric") )])
        self.reportSettings("Filter", ([("Min cluster size", self.minNumOfInstances)] if self.filterByNumOfInstances else []) + \
                                      ([("Max p-value", self.maxPValue_nofdr)] if self.filterByPValue_nofdr else []) + \
                                      ([("Max FDR", self.maxPValue)] if self.filterByPValue else []))

        def treeDepth(item):
            return 1 + max([treeDepth(item.child(i)) for i in range(item.childCount())] +[0])
        
        def printTree(item, level, treeDepth):
            text = '<tr>' + '<td width=16px></td>' * level
            text += '<td colspan="%i">%s: %s</td>' % (treeDepth - level, item.term.id, item.term.name)
            text += ''.join('<td>%s</td>' % item.text(i) for i in range(1, 5) + [6]) + '</tr>\n'
            for i in range(item.childCount()):
                text += printTree(item.child(i), level + 1, treeDepth)
            return text
        
        treeDepth = max([treeDepth(self.listView.topLevelItem(i)) for i in range(self.listView.topLevelItemCount())] + [0])
        
        tableText = '<table>\n<tr>' + ''.join('<th>%s</th>' % s for s in ["Term:", "List:", "Reference:", "p-value:", "FDR", "Enrichment:"]) + '</tr>'
        
        treeText = '<table>\n' +  '<th colspan="%i">%s</th>' % (treeDepth, "Term:") 
        treeText += ''.join('<th>%s</th>' % s for s in ["List:", "Reference:", "p-value:", "FDR", "Enrichment:"]) + '</tr>'
        
        for index in range(self.sigTerms.topLevelItemCount()):
            item = self.sigTerms.topLevelItem(index)
            tableText += printTree(item, 0, 1) 
        tableText += '</table>' 
        
        for index in range(self.listView.topLevelItemCount()):
            item = self.listView.topLevelItem(index)
            treeText += printTree(item, 0, treeDepth)
        
        self.reportSection("Enriched Terms")
        self.reportRaw(tableText)
        
        self.reportSection("Enriched Terms in the Ontology Tree")
        self.reportRaw(treeText)
        
    def onStateChanged(self, stateType, id, text):
        if stateType == "Warning":
            self.listView._userMessage = text
            self.listView.viewport().update()
            
    def onDeleteWidget(self):
        """ Called before the widget is removed from the canvas.
        """
        self.annotations = None
        self.ontology = None
        gc.collect() # Force collection
        

fmtp = lambda score: "%0.5f" % score if score > 10e-4 else "%0.1e" % score
fmtpdet = lambda score: "%0.9f" % score if score > 10e-4 else "%0.5e" % score

class GOTreeWidgetItem(QTreeWidgetItem):
    def __init__(self, term, enrichmentResult, nClusterGenes, nRefGenes, maxFoldEnrichment, parent):
        QTreeWidgetItem.__init__(self, parent)
        self.term = term
        self.enrichmentResult = enrichmentResult
        self.nClusterGenes = nClusterGenes
        self.nRefGenes = nRefGenes
        self.maxFoldEnrichment = maxFoldEnrichment
        self.enrichment = enrichment = lambda t:float(len(t[0])) / t[2] * (float(nRefGenes)/nClusterGenes)
        self.setText(0, term.name)
        fmt = "%" + str(-int(math.log(nClusterGenes))) + "i (%.2f%%)"
        self.setText(1, fmt % (len(enrichmentResult[0]), 100.0*len(self.enrichmentResult[0])/nClusterGenes))
        fmt = "%" + str(-int(math.log(nRefGenes))) + "i (%.2f%%)"
        self.setText(2, fmt % (enrichmentResult[2], 100.0*enrichmentResult[2]/nRefGenes))
        self.setText(3, fmtp(enrichmentResult[1]))
        self.setToolTip(3, fmtpdet(enrichmentResult[1]))
        self.setText(4, fmtp(enrichmentResult[3])) #FDR
        self.setToolTip(4, fmtpdet(enrichmentResult[3]))
        self.setText(5, ", ".join(enrichmentResult[0]))
        self.setText(6, "%.2f" % (enrichment(enrichmentResult)))
        self.setToolTip(6, "%.2f" % (enrichment(enrichmentResult)))
        self.setToolTip(0, "<p>" + term.__repr__()[6:].strip().replace("\n", "<br>"))
        self.sortByData = [term.name, len(self.enrichmentResult[0]), enrichmentResult[2], enrichmentResult[1], enrichmentResult[3], ", ".join(enrichmentResult[0]), enrichment(enrichmentResult)]

    def data(self, col, role):
        if role == Qt.UserRole:
            return QVariant(self.enrichment(self.enrichmentResult) / self.maxFoldEnrichment)
        else:
            return QTreeWidgetItem.data(self, col, role)

    def __lt__(self, other):
        col = self.treeWidget().sortColumn()
        return self.sortByData[col] < other.sortByData[col]
    
class EnrichmentColumnItemDelegate(QItemDelegate):
    def paint(self, painter, option, index):
        self.drawBackground(painter, option, index)
        value, ok = index.data(Qt.UserRole).toDouble()
        if ok:
            painter.save()
            painter.setBrush(QBrush(Qt.white, Qt.SolidPattern))
            painter.drawRect(option.rect)
            painter.setBrush(QBrush(Qt.blue, Qt.SolidPattern))
            painter.drawRect(option.rect.x(), option.rect.y(), value*(option.rect.width()-1), option.rect.height()-1)
            painter.restore()
        else:
            QItemDelegate.paint(self, painter, option, index)
        
        
class GeneMatcherDialog(OWWidget):
    items = [("useGO", "Use gene names from Gene Ontology annotations"),
             ("useKEGG", "Use gene names from KEGG Genes database"),
             ("useNCBI", "Use gene names from NCBI Gene info database"),
             ("useAffy", "Use Affymetrix platform reference ids")]
    settingsList = [item[0] for item in items]
    def __init__(self, parent=None, defaults=[True, False, False, False], enabled=[False, True, True, True], **kwargs):
        OWWidget.__init__(self, parent, **kwargs)
        for item, default in zip(self.items, defaults):
            setattr(self, item[0], default)
            
        self.loadSettings()
        for item, enable in zip(self.items, enabled):
            cb = OWGUI.checkBox(self, self, *item)
            cb.setEnabled(enable)
            
        box = OWGUI.widgetBox(self, orientation="horizontal")
        OWGUI.button(box, self, "OK", callback=self.accept)
        OWGUI.button(box, self, "Cancel", callback=self.reject)
        
        
if __name__=="__main__":
    import sys
    app = QApplication(sys.argv)
    w=OWGOEnrichmentAnalysis()
    data = orange.ExampleTable("brown-selected.tab")
    w.show()
    w.SetClusterDataset(data)
    app.exec_()
    w.saveSettings()
        
