"""
<name>KEGG Pathways</name>
<description>Browse KEGG pathways that include an input set of genes.</description>
<priority>2030</priority>
<icon>icons/KEGGPathways.svg</icon>
"""

from __future__ import absolute_import, with_statement

import sys
import gc
import webbrowser
from operator import add, itemgetter
from contextlib import contextmanager

from PyQt4.QtGui import (
    QTreeWidget, QTreeWidgetItem, QItemSelectionModel, QSplitter,
    QAction, QMenu, QGraphicsView, QGraphicsScene, QFont,
    QBrush, QColor, QPen, QTransform, QPainter, QPainterPath,
    QGraphicsItem, QGraphicsPathItem, QGraphicsPixmapItem, QPixmap,
    QMessageBox, QKeySequence
)

from PyQt4.QtCore import (
    Qt, QObject, QMetaObject, QTimer, Q_ARG, QRectF, SIGNAL, pyqtSlot
)

import Orange

from Orange.utils import serverfiles
from Orange.orng.orngDataCaching import data_hints

from Orange.OrangeWidgets import OWGUI
from Orange.OrangeWidgets.OWWidget import *
from Orange.OrangeWidgets.OWItemModels import VariableListModel
from Orange.OrangeWidgets.OWConcurrent import \
    ThreadExecutor, Task, methodinvoke


from .. import kegg
from .. import geneset


NAME = "KEGG Pathways"
DESCRIPTION = "Browse KEGG pathways that include an input set of genes."
ICON = "icons/KEGGPathways.svg"
PRIORITY = 2030

INPUTS = [("Examples", Orange.data.Table, "SetData", Default),
          ("Reference", Orange.data.Table, "SetRefData")]
OUTPUTS = [("Selected Examples", Orange.data.Table, Default),
           ("Unselected Examples", Orange.data.Table)]

REPLACES = ["_bioinformatics.widgets.OWKEGGPathwayBrowser.OWKEGGPathwayBrowser"]


def split_and_strip(string, sep=None):
    return [s.strip() for s in string.split(sep)]


def path_from_graphics(graphics):
    """
    Return a constructed `QPainterPath` for a KEGG pathway graphics element.
    """
    path = QPainterPath()
    x, y, w, h = [int(graphics.get(c, 0)) for c in
                  ["x", "y", "width", "height"]]
    type = graphics.get("type", "rectangle")
    if type == "rectangle":
        path.addRect(QRectF(x - w / 2, y - h / 2, w, h))
    elif type == "roundrectangle":
        path.addRoundedRect(QRectF(x - w / 2, y - h / 2, w, h), 10, 10)
    elif type == "circle":
        path.addEllipse(QRectF(x - w / 2, y - h / 2, w, h))
    else:
        ValueError("Unknown graphics type %r." % type)
    return path


class EntryGraphicsItem(QGraphicsPathItem):
    """
    An Graphics Item with actions for an overlay of a KEGG pathway image.
    """
    def __init__(self, graphics, *args):
        QGraphicsPathItem.__init__(self, *args)
        path = path_from_graphics(graphics)
        self.setPath(path)
        self.setAcceptHoverEvents(True)
        self._actions = []
        self.link = None

    def hoverEnterEvent(self, event):
        self.setBrush(QBrush(QColor(0, 100, 0, 100)))

    def hoverLeaveEvent(self, event):
        self.setBrush(QBrush(Qt.NoBrush))

    def contextMenuEvent(self, event):
        if self._actions:
            self._menu = menu = QMenu()
            for action in self._actions:
                menu.addAction(action)
            menu.popup(event.screenPos())

    def itemChange(self, change, value):
        if change == QGraphicsItem.ItemSelectedHasChanged:
            self.setPen(QPen(Qt.red if self.isSelected() else Qt.blue, 2))

        return QGraphicsPathItem.itemChange(self, change, value)


class GraphicsPathwayItem(QGraphicsPixmapItem):
    """
    A Graphics Item displaying a KEGG Pathway image with optional
    marked objects.

    """
    def __init__(self, pathway, objects, *args, **kwargs):
        QGraphicsPixmapItem.__init__(self, *args)
        self.setTransformationMode(Qt.SmoothTransformation)
        self.setPathway(pathway)
        self.setMarkedObjects(objects,
                              name_mapper=kwargs.get("name_mapper", {}))

    def setPathway(self, pathway):
        """
        Set pathway
        """
        self.pathway = pathway
        if pathway:
            image_filename = pathway.get_image()
            self._pixmap = QPixmap(image_filename)
        else:
            self._pixmap = QPixmap()
        self.setPixmap(self._pixmap)

    def setMarkedObjects(self, objects, name_mapper={}):
        for entry in self.pathway.entries() if self.pathway else []:
            if entry.type == "group":
                continue
            graphics = entry.graphics
            contained_objects = [obj for obj in objects if obj in entry.name]
            item = EntryGraphicsItem(graphics, self, self.scene())
            item.setToolTip(self.tooltip(entry, contained_objects,
                                         name_mapper))
            item._actions = self.actions(entry, contained_objects)
            item.marked_objects = contained_objects
            if contained_objects:
                item.setPen(QPen(Qt.blue, 2))
                item.setFlag(QGraphicsItem.ItemIsSelectable, True)

    def actions(self, entry, marked_objects=[]):
        actions = []
        type = entry.type
        if marked_objects:
            action = QAction("View genes on kegg website", None)
            org = set([s.split(":")[0] for s in marked_objects]).pop()
            genes = [s.split(":")[-1] for s in marked_objects]
            address = ("http://www.genome.jp/dbget-bin/www_bget?" +
                       "+".join([org] + genes))
            action.connect(action,
                           SIGNAL("triggered()"),
                           lambda toggled=False, address=address:
                               webbrowser.open(address))
            actions.append(action)
        elif hasattr(entry, "link"):
            action = QAction("View %s on KEGG website" % str(type), None)
            action.connect(action,
                           SIGNAL("triggered()"),
                           lambda toggled=False, address=entry.link: \
                               webbrowser.open(address))
            actions.append(action)
        return actions

    def tooltip(self, entry, objects, name_mapper={}):
        names = [obj for obj in objects if obj in entry.name]
        names = [name_mapper.get(name, name) for name in names]
        text = entry.name[:16] + " ..." if len(entry.name) > 20 else entry.name
        text = "<p>%s</p>" % text
        if names:
            text += "<br>".join(names)
        return text

    def contextMenuEvent(self, event):
        self._menu = menu = QMenu()
        action = menu.addAction("View this pathway on KEGG website")
        address = ("http://www.kegg.jp/kegg-bin/show_pathway?%s%s" %
                   (self.pathway.org, self.pathway.number))
        action.connect(action, SIGNAL("triggered()"),
                       lambda: webbrowser.open(address))
        menu.popup(event.screenPos())


class PathwayView(QGraphicsView):
    def __init__(self, master, *args):
        QGraphicsView.__init__(self, *args)
        self.master = master

        self.setHorizontalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        self.setVerticalScrollBarPolicy(Qt.ScrollBarAsNeeded)

        self.setRenderHints(QPainter.Antialiasing)
        scene = QGraphicsScene(self)
        self.pixmapGraphicsItem = QGraphicsPixmapItem(None, scene)
        self.setScene(scene)

        self.setMouseTracking(True)
        self.viewport().setMouseTracking(True)

        self.setFocusPolicy(Qt.WheelFocus)

    def SetPathway(self, pathway=None, objects=[]):
        self.scene().clear()
        self.pathway = pathway
        self.objects = objects
        self.pathwayItem = GraphicsPathwayItem(
            pathway, objects, None,
            name_mapper=getattr(self.master, "uniqueGenesDict", {})
        )

        self.scene().addItem(self.pathwayItem)
        self.scene().setSceneRect(self.pathwayItem.boundingRect())
        self.updateTransform()

    def resizeEvent(self, event):
        self.updateTransform()
        return QGraphicsView.resizeEvent(self, event)

    def updateTransform(self):
        if self.master.autoResize:
            self.fitInView(self.scene().sceneRect().adjusted(-1, -1, 1, 1),
                           Qt.KeepAspectRatio)
        else:
            self.setTransform(QTransform())

    def paintEvent(self, event):
        QGraphicsView.paintEvent(self, event)
        if getattr(self, "_userMessage", None):
            painter = QPainter(self.viewport())
            font = QFont(self.font())
            font.setPointSize(15)
            painter.setFont(font)
            painter.drawText(self.viewport().geometry(), Qt.AlignCenter,
                             self._userMessage)
            painter.end()


class OWKEGGPathwayBrowser(OWWidget):
    settingsList = ["organismIndex", "geneAttrIndex", "autoCommit",
                    "autoResize", "useReference", "useAttrNames",
                    "showOrthology"]

    contextHandlers = {
        "": DomainContextHandler(
            "",
            [ContextField("organismIndex",
                          DomainContextHandler.Required +
                          DomainContextHandler.IncludeMetaAttributes),
             ContextField("geneAttrIndex",
                          DomainContextHandler.Required +
                          DomainContextHandler.IncludeMetaAttributes),
             ContextField("useAttrNames",
                          DomainContextHandler.Required +
                          DomainContextHandler.IncludeMetaAttributes)]
        )
    }

    def __init__(self, parent=None, signalManager=None, name="KEGG Pathways"):
        OWWidget.__init__(self, parent, signalManager, name, wantGraph=True)
        self.inputs = [("Examples", Orange.data.Table, self.SetData),
                       ("Reference", Orange.data.Table, self.SetRefData)]
        self.outputs = [("Selected Examples", Orange.data.Table),
                        ("Unselected Examples", Orange.data.Table)]
        self.organismIndex = 0
        self.geneAttrIndex = 0
        self.autoCommit = False
        self.autoResize = True
        self.useReference = False
        self.useAttrNames = 0
        self.showOrthology = True

        self.loadSettings()

        self.organismCodes = []
        self._changedFlag = False

        self.controlArea.setMaximumWidth(250)
        box = OWGUI.widgetBox(self.controlArea, "Info")
        self.infoLabel = OWGUI.widgetLabel(box, "No data on input\n")

        # Organism selection.
        box = OWGUI.widgetBox(self.controlArea, "Organism")
        self.organismComboBox = OWGUI.comboBox(
            box, self, "organismIndex",
            items=[],
            callback=self.Update,
            addSpace=True,
            debuggingEnabled=0,
            tooltip="Select the organism of the input genes")

        # Selection of genes attribute
        box = OWGUI.widgetBox(self.controlArea, "Gene attribute")
        self.geneAttrCandidates = VariableListModel(parent=self)
        self.geneAttrCombo = OWGUI.comboBox(
            box, self, "geneAttrIndex", callback=self.Update)
        self.geneAttrCombo.setModel(self.geneAttrCandidates)

        OWGUI.checkBox(box, self, "useAttrNames",
                       "Use variable names",
                       disables=[(-1, self.geneAttrCombo)],
                       callback=self.Update)

        self.geneAttrCombo.setDisabled(bool(self.useAttrNames))

        OWGUI.separator(self.controlArea)

        OWGUI.checkBox(self.controlArea, self, "useReference",
                       "From signal",
                       box="Reference",
                       callback=self.Update)

        OWGUI.separator(self.controlArea)

        OWGUI.checkBox(self.controlArea, self, "showOrthology",
                       "Show pathways in full orthology",
                       box="Orthology",
                       callback=self.UpdateListView)

        OWGUI.checkBox(self.controlArea, self, "autoResize",
                       "Resize to fit",
                       box="Image",
                       callback=self.UpdatePathwayViewTransform)

        box = OWGUI.widgetBox(self.controlArea, "Cache Control")

        OWGUI.button(box, self, "Clear cache",
                     callback=self.ClearCache,
                     tooltip="Clear all locally cached KEGG data.")

        OWGUI.separator(self.controlArea)

        box = OWGUI.widgetBox(self.controlArea, "Selection")
        cb = OWGUI.checkBox(box, self, "autoCommit", "Commit on update")
        button = OWGUI.button(box, self, "Commit", callback=self.Commit,
                              default=True)
        OWGUI.setStopper(self, button, cb, "_changedFlag", self.Commit)

        OWGUI.rubber(self.controlArea)

        spliter = QSplitter(Qt.Vertical, self.mainArea)
        self.pathwayView = PathwayView(self, spliter)
        self.pathwayView.scene().selectionChanged.connect(
            self._onSelectionChanged
        )
        self.mainArea.layout().addWidget(spliter)

        self.listView = QTreeWidget(spliter)
        spliter.addWidget(self.listView)

        self.listView.setAllColumnsShowFocus(1)
        self.listView.setColumnCount(4)
        self.listView.setHeaderLabels(["Pathway", "P value",
                                       "Genes", "Reference"])

        self.listView.setSelectionMode(QTreeWidget.SingleSelection)

        self.listView.setSortingEnabled(True)

        self.listView.setMaximumHeight(200)

        self.connect(self.listView,
                     SIGNAL("itemSelectionChanged()"),
                     self.UpdatePathwayView)

        self.connect(self.graphButton,
                     SIGNAL("clicked()"),
                     self.saveGraph)

        select = QAction(
            "Select All", self,
            shortcut=QKeySequence.SelectAll
        )
        select.triggered.connect(self.selectAll)
        self.addAction(select)

        self.data = None
        self.refData = None

        self.resize(800, 600)

        self.connect(self,
                     SIGNAL("widgetStateChanged(QString, int, QString)"),
                     self.onStateChange)

        self.has_new_data = False
        self.has_new_reference_set = False

        self._executor = ThreadExecutor()
        self.setEnabled(False)
        self.setBlocking(True)
        QTimer.singleShot(0, self._initialize)
        self.infoLabel.setText("Fetching organism definitions\n")

    def _initialize(self):
        # First try to import slumber to see if we can even use the
        # kegg module.
        try:
            import slumber
        except ImportError:
            QMessageBox.warning(self,
                "'slumber' library required.",
                '<p>Please install '
                '<a href="http://pypi.python.org/pypi/slumber">slumber</a> '
                'library to use KEGG Pathways widget.</p>'
            )
            self.infoLabel.setText(
                '<p>Please install '
                '<a href="http://pypi.python.org/pypi/slumber">slumber</a> '
                'library to use KEGG Pathways widget.</p>'
            )
            self.error(0, "Missing slumber/requests library")
            return

        progress = methodinvoke(self, "setProgress", (float,))

        def get_genome():
            """Return a KEGGGenome with the common org entries precached."""
            genome = kegg.KEGGGenome()

            essential = genome.essential_organisms()
            common = genome.common_organisms()
            # Remove duplicates of essential from common.
            # (essential + common list as defined here will be used in the
            # GUI.)
            common = [c for c in common if c not in essential]

            # TODO: Add option to specify additional organisms not
            # in the common list.

            keys = map(genome.org_code_to_entry_key, essential + common)

            genome.pre_cache(keys, progress_callback=progress)
            return (keys, genome)

        self._genomeTask = task = Task(function=get_genome)
        task.finished.connect(self._initializeOrganisms)

        self.progressBarInit()
        self._executor.submit(task)

    def _initializeOrganisms(self):
        self.progressBarFinished()
        try:
            keys, genome = self._genomeTask.result()
        except Exception as err:
            self.error(0, str(err))
            return

        entries = [genome[key] for key in keys]
        items = [entry.definition for entry in entries]
        codes = [entry.organism_code for entry in entries]

        self.organismCodes = codes
        self.organismComboBox.clear()
        self.organismComboBox.addItems(items)
        self.organismComboBox.setCurrentIndex(self.organismIndex)

        self.setEnabled(True)
        self.setBlocking(False)
        self.infoLabel.setText("No data on input\n")

    def Clear(self):
        """
        Clear the widget state.
        """
        self.queryGenes = []
        self.referenceGenes = []
        self.genes = {}
        self.uniqueGenesDict = {}
        self.revUniqueGenesDict = {}
        self.pathways = {}
        self.org = None
        self.geneAttrCandidates[:] = []

        self.infoLabel.setText("No data on input\n")
        self.listView.clear()
        self.pathwayView.SetPathway(None)

        self.send("Selected Examples", None)
        self.send("Unselected Examples", None)

    def SetData(self, data=None):
        self.closeContext()
        self.data = data
        self.warning(0)
        self.error(0)
        self.information(0)

        if data is not None:
            vars = data.domain.variables + data.domain.getmetas().values()
            vars = [var for var in vars
                    if isinstance(var, (Orange.feature.String,
                                        Orange.feature.Discrete))]
            self.geneAttrCandidates[:] = vars

            # Try to guess the gene name variable
            names_lower = [v.name.lower() for v in vars]
            scores = [(name == "gene", "gene" in name)
                      for name in names_lower]
            imax, _ = max(enumerate(scores), key=itemgetter(1))
            self.geneAttrIndex = imax

            taxid = data_hints.get_hint(data, "taxid", None)
            if taxid:
                try:
                    code = kegg.from_taxid(taxid)
                    self.organismIndex = self.organismCodes.index(code)
                except Exception as ex:
                    print ex, taxid

            self.useAttrNames = data_hints.get_hint(data, "genesinrows",
                                                    self.useAttrNames)

            self.openContext("", data)
        else:
            self.Clear()

        self.has_new_data = True

    def SetRefData(self, data=None):
        self.refData = data
        self.information(1)

        self.has_new_reference_set = True

    def handleNewSignals(self):
        if self.has_new_data or (self.has_new_reference_set and \
                                 self.useReference):
            self.Update()

            self.has_new_data = False
            self.has_new_reference_set = False

    def UpdateListView(self):
        self.bestPValueItem = None
        self.listView.clear()
        if not self.data:
            return

        allPathways = self.org.pathways()
        allRefPathways = kegg.pathways("map")

        items = []
        kegg_pathways = kegg.KEGGPathways()

        org_code = self.organismCodes[min(self.organismIndex,
                                          len(self.organismCodes) - 1)]

        if self.showOrthology:
            self.koOrthology = kegg.KEGGBrite("ko00001")
            self.listView.setRootIsDecorated(True)
            path_ids = set([s[-5:] for s in self.pathways.keys()])

            def _walkCollect(koEntry):
                num = koEntry.title[:5] if koEntry.title else None
                if num in path_ids:
                    return ([koEntry] +
                            reduce(lambda li, c: li + _walkCollect(c),
                                   [child for child in koEntry.entries],
                                   []))
                else:
                    c = reduce(lambda li, c: li + _walkCollect(c),
                               [child for child in koEntry.entries],
                               [])
                    return c + (c and [koEntry] or [])

            allClasses = reduce(lambda li1, li2: li1 + li2,
                                [_walkCollect(c) for c in self.koOrthology],
                                [])

            def _walkCreate(koEntry, lvItem):
                item = QTreeWidgetItem(lvItem)
                id = "path:" + org_code + koEntry.title[:5]

                if koEntry.title[:5] in path_ids:
                    p = kegg_pathways.get_entry(id)
                    if p is None:
                        # In case the genesets still have obsolete entries
                        name = koEntry.title
                    else:
                        name = p.name
                    genes, p_value, ref = self.pathways[id]
                    item.setText(0, name)
                    item.setText(1, "%.5f" % p_value)
                    item.setText(2, "%i of %i" % (len(genes), len(self.genes)))
                    item.setText(3, "%i of %i" % (ref, len(self.referenceGenes)))
                    item.pathway_id = id if p is not None else None
                else:
                    if id in allPathways:
                        text = kegg_pathways.get_entry(id).name
                    else:
                        text = koEntry.title
                    item.setText(0, text)

                    if id in allPathways:
                        item.pathway_id = id
                    elif "path:map" + koEntry.title[:5] in allRefPathways:
                        item.pathway_id = "path:map" + koEntry.title[:5]
                    else:
                        item.pathway_id = None

                for child in koEntry.entries:
                    if child in allClasses:
                        _walkCreate(child, item)

            for koEntry in self.koOrthology:
                if koEntry in allClasses:
                    _walkCreate(koEntry, self.listView)

            self.listView.update()
        else:
            self.listView.setRootIsDecorated(False)
            pathways = self.pathways.items()
            pathways.sort(lambda a, b: cmp(a[1][1], b[1][1]))

            for id, (genes, p_value, ref) in pathways:
                item = QTreeWidgetItem(self.listView)
                item.setText(0, kegg_pathways.get_entry(id).name)
                item.setText(1, "%.5f" % p_value)
                item.setText(2, "%i of %i" % (len(genes), len(self.genes)))
                item.setText(3, "%i of %i" % (ref, len(self.referenceGenes)))
                item.pathway_id = id
                items.append(item)

        self.bestPValueItem = items and items[0] or None
        self.listView.expandAll()
        for i in range(4):
            self.listView.resizeColumnToContents(i)

        if self.bestPValueItem:
            index = self.listView.indexFromItem(self.bestPValueItem)
            self.listView.selectionModel().select(
                index, QItemSelectionModel.ClearAndSelect
            )

    def UpdatePathwayView(self):
        items = self.listView.selectedItems()

        if len(items) > 0:
            item = items[0]
        else:
            item = None

        self.Commit()
        item = item or self.bestPValueItem
        if not item or not item.pathway_id:
            self.pathwayView.SetPathway(None)
            return

        def get_kgml_and_image(pathway_id):
            """Return an initialized KEGGPathway with pre-cached data"""
            p = kegg.KEGGPathway(pathway_id)
            p._get_kgml()  # makes sure the kgml file is downloaded
            p._get_image_filename()  # makes sure the image is downloaded
            return (pathway_id, p)

        self.setEnabled(False)
        self._pathwayTask = Task(
            function=lambda: get_kgml_and_image(item.pathway_id)
        )
        self._pathwayTask.finished.connect(self._onPathwayTaskFinshed)
        self._executor.submit(self._pathwayTask)

    def _onPathwayTaskFinshed(self):
        self.setEnabled(True)
        pathway_id, self.pathway = self._pathwayTask.result()
        self.pathwayView.SetPathway(
            self.pathway,
            self.pathways.get(pathway_id, [[]])[0]
        )

    def UpdatePathwayViewTransform(self):
        self.pathwayView.updateTransform()

    def Update(self):
        """
        Update (recompute enriched pathways) the widget state.
        """
        if not self.data:
            return

        self.error(0)
        self.information(0)

        # XXX: Check data in setData, do not even alow this to be executed if
        # data has no genes
        try:
            genes = self.GeneNamesFromData(self.data)
        except ValueError:
            self.error(0, "Cannot extract gene names from input.")
            genes = []

        if not self.useAttrNames and any("," in gene for gene in genes):
            genes = reduce(add, (split_and_strip(gene, ",")
                                 for gene in genes),
                           [])
            self.information(0,
                             "Separators detected in input gene names. "
                             "Assuming multiple genes per instance.")
        self.queryGenes = genes

        self.information(1)
        reference = None
        if self.useReference and self.refData:
            reference = self.GeneNamesFromData(self.refData)
            if not self.useAttrNames \
                    and any("," in gene for gene in reference):
                reference = reduce(add, (split_and_strip(gene, ",")
                                         for gene in reference),
                                   [])
                self.information(1,
                                 "Separators detected in reference gene "
                                 "names. Assuming multiple genes per "
                                 "instance.")

        org_code = self.SelectedOrganismCode()

        def run_enrichment(org_code, genes, reference=None, progress=None):
            org = kegg.KEGGOrganism(org_code)
            if reference is None:
                reference = org.get_genes()

            # Map 'genes' and 'reference' sets to unique KEGG identifiers
            unique_genes, _, _ = org.get_unique_gene_ids(set(genes))
            unique_ref_genes, _, _ = org.get_unique_gene_ids(set(reference))

            taxid = kegg.to_taxid(org.org_code)
            # Map the taxid back to standard 'common' taxids
            # (as used by 'geneset') if applicable
            r_tax_map = dict((v, k) for k, v in
                             kegg.KEGGGenome.TAXID_MAP.items())
            if taxid in r_tax_map:
                taxid = r_tax_map[taxid]

            # We use the kegg pathway gene sets provided by 'geneset' for
            # the enrichment calculation.

            # Ensure we are using the latest genesets
            # TODO: ?? Is updating the index enough?
            serverfiles.update(geneset.sfdomain, "index.pck")
            kegg_gs_collections = geneset.collections(
                (("KEGG", "pathways"), taxid)
            )

            pathways = pathway_enrichment(
                kegg_gs_collections, unique_genes.keys(),
                unique_ref_genes.keys(),
                callback=progress
            )
            # Ensure that pathway entries are pre-cached for later use in the
            # list/tree view
            kegg_pathways = kegg.KEGGPathways()
            kegg_pathways.pre_cache(
                pathways.keys(), progress_callback=progress
            )

            return pathways, org, unique_genes, unique_ref_genes

        self.progressBarInit()
        self.setEnabled(False)
        self.infoLabel.setText("Retrieving...\n")

        progress = methodinvoke(self, "setProgress", (float,))
        self._enrichTask = Task(
            function=lambda:
                run_enrichment(org_code, genes, reference, progress)
        )
        self._enrichTask.finished.connect(self._onEnrichTaskFinished)
        self._executor.submit(self._enrichTask)

    def _onEnrichTaskFinished(self):
        self.setEnabled(True)
        self.setBlocking(False)

        self.progressBarFinished()
        try:
            pathways, org, unique_genes, unique_ref_genes = \
                self._enrichTask.result()
        except Exception:
            raise

        self.org = org
        self.genes = unique_genes.keys()
        self.uniqueGenesDict = unique_genes
        self.revUniqueGenesDict = dict([(val, key) for key, val in
                                        self.uniqueGenesDict.items()])
        self.referenceGenes = unique_ref_genes.keys()
        self.pathways = pathways

        if not self.pathways:
            self.warning(0, "No enriched pathways found.")
        else:
            self.warning(0)

        count = len(set(self.queryGenes))
        self.infoLabel.setText(
            "%i unique gene names on input\n"
            "%i (%.1f%%) genes names matched" %
            (count, len(unique_genes),
             100.0 * len(unique_genes) / count if count else 0.0)
        )

        self.UpdateListView()

    @pyqtSlot(float)
    def setProgress(self, value):
        self.progressBarValue = value

    def GeneNamesFromData(self, data):
        """
        Extract and return gene names from `data`.
        """
        if self.useAttrNames:
            genes = [str(v.name).strip() for v in data.domain.attributes]
        elif self.geneAttrCandidates:
            index = min(self.geneAttrIndex, len(self.geneAttrCandidates) - 1)
            geneAttr = self.geneAttrCandidates[index]
            genes = [str(e[geneAttr]) for e in data
                     if not e[geneAttr].isSpecial()]
        else:
            raise ValueError("No gene names in data.")
        return genes

    def SelectedOrganismCode(self):
        """
        Return the selected organism code.
        """
        return self.organismCodes[min(self.organismIndex,
                                      len(self.organismCodes) - 1)]

    def selectAll(self):
        """
        Select all items in the pathway view.
        """
        changed = False
        scene = self.pathwayView.scene()
        with disconnected(scene.selectionChanged, self._onSelectionChanged):
            for item in scene.items():
                if item.flags() & QGraphicsItem.ItemIsSelectable and \
                        not item.isSelected():
                    item.setSelected(True)
                    changed = True
        if changed:
            self._onSelectionChanged()

    def _onSelectionChanged(self):
        # Item selection in the pathwayView/scene has changed
        if self.autoCommit:
            self.Commit()
        else:
            self._changedFlag = True

    def Commit(self):
        if self.data:
            selectedItems = self.pathwayView.scene().selectedItems()
            selectedGenes = reduce(set.union, [item.marked_objects
                                               for item in selectedItems],
                                   set())

            if self.useAttrNames:
                selectedVars = [self.data.domain[self.uniqueGenesDict[gene]]
                                for gene in selectedGenes]
                newDomain = Orange.data.Domain(selectedVars, 0)
                data = Orange.data.Table(newDomain, self.data)
                self.send("Selected Examples", data)
            elif self.geneAttrCandidates:
                geneAttr = self.geneAttrCandidates[min(self.geneAttrIndex,
                                                       len(self.geneAttrCandidates) - 1)]
                selectedExamples = []
                otherExamples = []
                for ex in self.data:
                    names = [self.revUniqueGenesDict.get(name, None)
                             for name in split_and_strip(str(ex[geneAttr]), ",")]
                    if any(name and name in selectedGenes for name in names):
                        selectedExamples.append(ex)
                    else:
                        otherExamples.append(ex)

                if selectedExamples:
                    selectedExamples = Orange.data.Table(selectedExamples)
                else:
                    selectedExamples = None

                if otherExamples:
                    otherExamples = Orange.data.Table(otherExamples)
                else:
                    otherExamples = None

                self.send("Selected Examples", selectedExamples)
                self.send("Unselected Examples", otherExamples)
        else:
            self.send("Selected Examples", None)
            self.send("Unselected Examples", None)

    def ClearCache(self):
        from ..kegg import caching
        try:
            caching.clear_cache()
        except Exception, ex:
            QMessageBox.warning(self, "Cache clear", ex.args[0])

    def onStateChange(self, stateType, id, text):
        if stateType == "Warning":
            self.pathwayView._userMessage = text
            self.pathwayView.viewport().update()

    def saveGraph(self):
        from Orange.OrangeWidgets.OWDlgs import OWChooseImageSizeDlg
        sizeDlg = OWChooseImageSizeDlg(self.pathwayView.scene(), parent=self)
        sizeDlg.exec_()

    def onDeleteWidget(self):
        """
        Called before the widget is removed from the canvas.
        """
        self.org = None
        gc.collect()  # Force collection

        self._executor.shutdown(wait=False)


from ..utils import stats


def pathway_enrichment(genesets, genes, reference, prob=None, callback=None):
    result_sets = []
    p_values = []
    if prob is None:
        prob = stats.Hypergeometric()

    for i, gs in enumerate(genesets):
        cluster = gs.genes.intersection(genes)
        ref = gs.genes.intersection(reference)
        k = len(cluster)
        N = len(reference)
        m = len(ref)
        n = len(genes)
        if k:
            p_val = prob.p_value(k, N, m, n)
            result_sets.append((gs.id, cluster, ref))
            p_values.append(p_val)
        if callback is not None:
            callback(100.0 * i / len(genesets))

    # FDR correction
    p_values = stats.FDR(p_values)

    return dict([(id, (genes, p_val, len(ref)))
                 for (id, genes, ref), p_val in zip(result_sets, p_values)])


@contextmanager
def disconnected(signal, slot):
    """
    A context manager disconnecting a slot from a signal.
    ::

        with disconnected(scene.selectionChanged, self.onSelectionChanged):
            # Can change item selection in a scene without
            # onSelectionChanged being invoked.
            do_something()

    """
    signal.disconnect(slot)
    try:
        yield
    finally:
        signal.connect(slot)


if __name__ == "__main__":
    app = QApplication(sys.argv)
    w = OWKEGGPathwayBrowser()
    w.show()

    data = Orange.data.Table("brown-selected.tab")

    QTimer.singleShot(3000, lambda: w.SetData(data))
    QTimer.singleShot(3500, w.handleNewSignals)

    app.exec_()
    w.saveSettings()
