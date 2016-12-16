"""
Set Enrichment
--------------

"""
import sys
import math
import operator
import itertools
import traceback
import types
import concurrent.futures
from collections import defaultdict
from functools import reduce, partial

import numpy as np

from AnyQt.QtWidgets import (
    QTreeWidget, QTreeWidgetItem, QTreeView, QLineEdit, QCompleter,
    QHBoxLayout, QStyle, QStyledItemDelegate, QApplication
)
from AnyQt.QtGui import (
    QBrush, QColor, QFont, QStandardItemModel, QStandardItem
)
from AnyQt.QtCore import (
    Qt, QRect, QSize, QModelIndex, QStringListModel, QThread, QThreadPool,
    Slot
)

import Orange
from Orange.widgets.utils.datacaching import data_hints
from Orange.widgets import widget, gui, settings
from Orange.widgets.utils.concurrent import ThreadExecutor, Task, methodinvoke

from orangecontrib.bio import gene, geneset, taxonomy, utils

from ..widgets.utils.download import EnsureDownloaded


def gsname(geneset):
    return geneset.name if geneset.name else geneset.id

fmtp = lambda score: "%0.5f" % score if score > 10e-4 else "%0.1e" % score
fmtpdet = lambda score: "%0.9f" % score if score > 10e-4 else "%0.5e" % score

# A translation table mapping punctuation, ... to spaces
_TR_TABLE = dict((ord(c), ord(" ")) for c in ".,!?()[]{}:;'\"<>")


def word_split(string):
    """
    Split a string into a list of words.
    """
    return string.translate(_TR_TABLE).split()


class BarItemDelegate(QStyledItemDelegate):
    def __init__(self, parent, brush=QBrush(QColor(255, 170, 127)),
                 scale=(0.0, 1.0)):
        super().__init__(parent)
        self.brush = brush
        self.scale = scale

    def paint(self, painter, option, index):
        if option.widget is not None:
            style = option.widget.style()
        else:
            style = QApplication.instance().style()

        style.drawPrimitive(
            QStyle.PE_PanelItemViewRow, option, painter, option.widget)
        style.drawPrimitive(
            QStyle.PE_PanelItemViewItem, option, painter, option.widget)
        rect = option.rect
        val = index.data(Qt.DisplayRole)

        if isinstance(val, float):
            if np.isfinite(val):
                minv, maxv = self.scale
                val = (val - minv) / (maxv - minv)
                rect = rect.adjusted(
                    1, 1, - int(rect.width() * (1.0 - val)) - 2, -2)
            else:
                rect = QRect()

            painter.save()
            if option.state & QStyle.State_Selected:
                painter.setOpacity(0.75)
            painter.setBrush(self.brush)
            painter.drawRect(rect)
            painter.restore()


def name_or_none(taxid):
    """Return organism name for ncbi taxid or None if not found.
    """
    try:
        return taxonomy.name(taxid)
    except taxonomy.UnknownSpeciesIdentifier:
        return None


def fulfill(value):
    """Return a pre-fulfilled Future yielding `value`."""
    f = concurrent.futures.Future()
    f.set_result(value)
    return f


def withtraceback(func):
    def f(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except Exception as ex:
            ex._traceback = traceback.format_exc()
            raise
    return f


def memoize(func):
    cache = {}

    def memoized(arg):
        try:
            return cache[arg]
        except KeyError:
            cache[arg] = value = func(arg)
            return value

    return memoized


class OWSetEnrichment(widget.OWWidget):
    name = "Set Enrichment"
    description = ""
    icon = "../widgets/icons/GeneSetEnrichment.svg"
    priority = 5000

    inputs = [("Data", Orange.data.Table, "setData", widget.Default),
              ("Reference", Orange.data.Table, "setReference")]
    outputs = [("Data subset", Orange.data.Table)]

    settingsHandler = settings.DomainContextHandler()

    taxid = settings.ContextSetting(None)
    speciesIndex = settings.ContextSetting(0)
    genesinrows = settings.ContextSetting(False)
    geneattr = settings.ContextSetting(0)
    categoriesCheckState = settings.ContextSetting({})

    useReferenceData = settings.Setting(False)
    useMinCountFilter = settings.Setting(True)
    useMaxPValFilter = settings.Setting(True)
    useMaxFDRFilter = settings.Setting(True)
    minClusterCount = settings.Setting(3)
    maxPValue = settings.Setting(0.01)
    maxFDR = settings.Setting(0.01)
    autocommit = settings.Setting(False)

    Ready, Initializing, Loading, RunningEnrichment = 0, 1, 2, 4

    def __init__(self, parent=None):
        super().__init__(parent)

        self.geneMatcherSettings = [False, False, True, False]

        self.data = None
        self.referenceData = None
        self.taxid_list = []

        self.__genematcher = (None, fulfill(gene.matcher([])))
        self.__invalidated = False

        self.currentAnnotatedCategories = []
        self.state = None
        self.__state = OWSetEnrichment.Initializing

        box = gui.widgetBox(self.controlArea, "Info")
        self.infoBox = gui.widgetLabel(box, "Info")
        self.infoBox.setText("No data on input.\n")

        self.speciesComboBox = gui.comboBox(
            self.controlArea, self,
            "speciesIndex", "Species",
            callback=self.__on_speciesIndexChanged)

        box = gui.widgetBox(self.controlArea, "Entity names")
        self.geneAttrComboBox = gui.comboBox(
            box, self, "geneattr", "Entity feature", sendSelectedValue=0,
            callback=self.updateAnnotations)

        cb = gui.checkBox(
            box, self, "genesinrows", "Use feature names",
            callback=self.updateAnnotations,
            disables=[(-1, self.geneAttrComboBox)])
        cb.makeConsistent()

#         gui.button(box, self, "Gene matcher settings",
#                    callback=self.updateGeneMatcherSettings,
#                    tooltip="Open gene matching settings dialog")

        self.referenceRadioBox = gui.radioButtonsInBox(
            self.controlArea,
            self, "useReferenceData",
            ["All entities", "Reference set (input)"],
            tooltips=["Use entire genome (for gene set enrichment) or all " +
                      "available entities for reference",
                      "Use entities from Reference Examples input signal " +
                      "as reference"],
            box="Reference", callback=self.updateAnnotations)

        box = gui.widgetBox(self.controlArea, "Entity Sets")
        self.groupsWidget = QTreeWidget(self)
        self.groupsWidget.setHeaderLabels(["Category"])
        box.layout().addWidget(self.groupsWidget)

        hLayout = QHBoxLayout()
        hLayout.setSpacing(10)
        hWidget = gui.widgetBox(self.mainArea, orientation=hLayout)
        gui.spin(hWidget, self, "minClusterCount",
                 0, 100, label="Entities",
                 tooltip="Minimum entity count",
                 callback=self.filterAnnotationsChartView,
                 callbackOnReturn=True,
                 checked="useMinCountFilter",
                 checkCallback=self.filterAnnotationsChartView)

        pvalfilterbox = gui.widgetBox(hWidget, orientation="horizontal")
        cb = gui.checkBox(
            pvalfilterbox, self, "useMaxPValFilter", "p-value",
            callback=self.filterAnnotationsChartView)

        sp = gui.doubleSpin(
            pvalfilterbox, self, "maxPValue", 0.0, 1.0, 0.0001,
            tooltip="Maximum p-value",
            callback=self.filterAnnotationsChartView,
            callbackOnReturn=True,
        )
        sp.setEnabled(self.useMaxFDRFilter)
        cb.toggled[bool].connect(sp.setEnabled)

        pvalfilterbox.layout().setAlignment(cb, Qt.AlignRight)
        pvalfilterbox.layout().setAlignment(sp, Qt.AlignLeft)

        fdrfilterbox = gui.widgetBox(hWidget, orientation="horizontal")
        cb = gui.checkBox(
            fdrfilterbox, self, "useMaxFDRFilter", "FDR",
            callback=self.filterAnnotationsChartView)

        sp = gui.doubleSpin(
            fdrfilterbox, self, "maxFDR", 0.0, 1.0, 0.0001,
            tooltip="Maximum False discovery rate",
            callback=self.filterAnnotationsChartView,
            callbackOnReturn=True,
        )
        sp.setEnabled(self.useMaxFDRFilter)
        cb.toggled[bool].connect(sp.setEnabled)

        fdrfilterbox.layout().setAlignment(cb, Qt.AlignRight)
        fdrfilterbox.layout().setAlignment(sp, Qt.AlignLeft)

        self.filterLineEdit = QLineEdit(
            self, placeholderText="Filter ...")

        self.filterCompleter = QCompleter(self.filterLineEdit)
        self.filterCompleter.setCaseSensitivity(Qt.CaseInsensitive)
        self.filterLineEdit.setCompleter(self.filterCompleter)

        hLayout.addWidget(self.filterLineEdit)
        self.mainArea.layout().addWidget(hWidget)

        self.filterLineEdit.textChanged.connect(
            self.filterAnnotationsChartView)

        self.annotationsChartView = QTreeView(
            alternatingRowColors=True,
            sortingEnabled=True,
            selectionMode=QTreeView.ExtendedSelection,
            rootIsDecorated=False,
            editTriggers=QTreeView.NoEditTriggers,
        )
        self.annotationsChartView.viewport().setMouseTracking(True)
        self.mainArea.layout().addWidget(self.annotationsChartView)

        contextEventFilter = gui.VisibleHeaderSectionContextEventFilter(
            self.annotationsChartView)
        self.annotationsChartView.header().installEventFilter(contextEventFilter)

        self.groupsWidget.itemClicked.connect(self.subsetSelectionChanged)
        gui.auto_commit(self.controlArea, self, "autocommit", "Commit")

        self.setBlocking(True)

        task = EnsureDownloaded(
            [(taxonomy.Taxonomy.DOMAIN, taxonomy.Taxonomy.FILENAME),
             (geneset.sfdomain, "index.pck")]
        )

        task.finished.connect(self.__initialize_finish)
        self.setStatusMessage("Initializing")
        self._executor = ThreadExecutor(
            parent=self, threadPool=QThreadPool(self))
        self._executor.submit(task)

    def sizeHint(self):
        return QSize(1024, 600)

    def __initialize_finish(self):
        # Finalize the the widget's initialization (preferably after
        # ensuring all required databases have been downloaded.

        sets = geneset.list_all()
        taxids = set(taxonomy.common_taxids() +
                     list(filter(None, [tid for _, tid, _ in sets])))
        organisms = [(tid, name_or_none(tid)) for tid in taxids]
        organisms = [(tid, name) for tid, name in organisms
                     if name is not None]

        organisms = [(None, "None")] + sorted(organisms)
        taxids = [tid for tid, _ in organisms]
        names = [name for _, name in organisms]
        self.taxid_list = taxids

        self.speciesComboBox.clear()
        self.speciesComboBox.addItems(names)
        self.genesets = sets

        if self.taxid in self.taxid_list:
            taxid = self.taxid
        else:
            taxid = self.taxid_list[0]

        self.taxid = None
        self.setCurrentOrganism(taxid)
        self.setBlocking(False)
        self.__state = OWSetEnrichment.Ready
        self.setStatusMessage("")

    def setCurrentOrganism(self, taxid):
        """Set the current organism `taxid`."""
        if taxid not in self.taxid_list:
            taxid = self.taxid_list[min(self.speciesIndex,
                                        len(self.taxid_list) - 1)]
        if self.taxid != taxid:
            self.taxid = taxid
            self.speciesIndex = self.taxid_list.index(taxid)
            self.refreshHierarchy()
            self._invalidateGeneMatcher()
            self._invalidate()

    def currentOrganism(self):
        """Return the current organism taxid"""
        return self.taxid

    def __on_speciesIndexChanged(self):
        taxid = self.taxid_list[self.speciesIndex]
        self.taxid = "< Do not look >"
        self.setCurrentOrganism(taxid)
        if self.__invalidated and self.data is not None:
            self.updateAnnotations()

    def clear(self):
        """Clear/reset the widget state."""
        self._cancelPending()
        self.state = None

        self.__state = self.__state & ~OWSetEnrichment.RunningEnrichment

        self._clearView()

        if self.annotationsChartView.model() is not None:
            self.annotationsChartView.model().clear()

        self.geneAttrComboBox.clear()
        self.geneAttrs = []
        self._updatesummary()

    def _cancelPending(self):
        """Cancel pending tasks."""
        if self.state is not None:
            self.state.results.cancel()
            self.state.namematcher.cancel()
            self.state.cancelled = True

    def _clearView(self):
        """Clear the enrichment report view (main area)."""
        if self.annotationsChartView.model() is not None:
            self.annotationsChartView.model().clear()

    def setData(self, data=None):
        """Set the input dataset with query gene names"""
        if self.__state & OWSetEnrichment.Initializing:
            self.__initialize_finish()

        self.error(0)
        self.closeContext()
        self.clear()

        self.groupsWidget.clear()
        self.data = data

        if data is not None:
            varlist = [var for var in data.domain.variables + data.domain.metas
                       if isinstance(var, Orange.data.StringVariable)]

            self.geneAttrs = varlist
            for var in varlist:
                self.geneAttrComboBox.addItem(*gui.attributeItem(var))

            oldtaxid = self.taxid
            self.geneattr = min(self.geneattr, len(self.geneAttrs) - 1)

            taxid = data_hints.get_hint(data, "taxid", "")
            if taxid in self.taxid_list:
                self.speciesIndex = self.taxid_list.index(taxid)
                self.taxid = taxid

            self.genesinrows = data_hints.get_hint(
                data, "genesinrows", self.genesinrows)

            self.openContext(data)
            if oldtaxid != self.taxid:
                self.taxid = "< Do not look >"
                self.setCurrentOrganism(taxid)

            self.refreshHierarchy()
            self._invalidate()

    def setReference(self, data=None):
        """Set the (optional) input dataset with reference gene names."""
        self.referenceData = data
        self.referenceRadioBox.setEnabled(bool(data))
        if self.useReferenceData:
            self._invalidate()

    def handleNewSignals(self):
        if self.__invalidated:
            self.updateAnnotations()

    def _invalidateGeneMatcher(self):
        _, f = self.__genematcher
        f.cancel()
        self.__genematcher = (None, fulfill(gene.matcher([])))

    def _invalidate(self):
        self.__invalidated = True

    def genesFromTable(self, table):
        if self.genesinrows:
            genes = [attr.name for attr in table.domain.attributes]
        else:
            geneattr = self.geneAttrs[self.geneattr]
            genes = [str(ex[geneattr]) for ex in table]
        return genes

    def getHierarchy(self, taxid):
        def recursive_dict():
            return defaultdict(recursive_dict)
        collection = recursive_dict()

        def collect(col, hier):
            if hier:
                collect(col[hier[0]], hier[1:])

        for hierarchy, t_id, _ in self.genesets:
            collect(collection[t_id], hierarchy)

        return (taxid, collection[taxid]), (None, collection[None])

    def setHierarchy(self, hierarchy, hierarchy_noorg):
        self.groupsWidgetItems = {}

        def fill(col, parent, full=(), org=""):
            for key, value in sorted(col.items()):
                full_cat = full + (key,)
                item = QTreeWidgetItem(parent, [key])
                item.setFlags(item.flags() | Qt.ItemIsUserCheckable |
                              Qt.ItemIsSelectable | Qt.ItemIsEnabled)
                if value:
                    item.setFlags(item.flags() | Qt.ItemIsTristate)

                checked = self.categoriesCheckState.get(
                    (full_cat, org), Qt.Checked)
                item.setData(0, Qt.CheckStateRole, checked)
                item.setExpanded(True)
                item.category = full_cat
                item.organism = org
                self.groupsWidgetItems[full_cat] = item
                fill(value, item, full_cat, org=org)

        self.groupsWidget.clear()
        fill(hierarchy[1], self.groupsWidget, org=hierarchy[0])
        fill(hierarchy_noorg[1], self.groupsWidget, org=hierarchy_noorg[0])

    def refreshHierarchy(self):
        self.setHierarchy(*self.getHierarchy(taxid=self.taxid_list[self.speciesIndex]))

    def selectedCategories(self):
        """
        Return a list of currently selected hierarchy keys.

        A key is a tuple of identifiers from the root to the leaf of
        the hierarchy tree.
        """
        return [key for key, check in self.getHierarchyCheckState().items()
                if check == Qt.Checked]

    def getHierarchyCheckState(self):
        def collect(item, full=()):
            checked = item.checkState(0)
            name = str(item.data(0, Qt.DisplayRole))
            full_cat = full + (name,)
            result = [((full_cat, item.organism), checked)]
            for i in range(item.childCount()):
                result.extend(collect(item.child(i), full_cat))
            return result

        items = [self.groupsWidget.topLevelItem(i)
                 for i in range(self.groupsWidget.topLevelItemCount())]
        states = itertools.chain(*(collect(item) for item in items))
        return dict(states)

    def subsetSelectionChanged(self, item, column):
        # The selected geneset (hierarchy) subset has been changed by the
        # user. Update the displayed results.
        # Update the stored state (persistent settings)
        self.categoriesCheckState = self.getHierarchyCheckState()
        categories = self.selectedCategories()

        if self.data is not None:
            if self._nogenematching() or \
                    not set(categories) <= set(self.currentAnnotatedCategories):
                self.updateAnnotations()
            else:
                self.filterAnnotationsChartView()

    def updateGeneMatcherSettings(self):
        raise NotImplementedError

        from .OWGOEnrichmentAnalysis import GeneMatcherDialog
        dialog = GeneMatcherDialog(self, defaults=self.geneMatcherSettings, enabled=[True] * 4, modal=True)
        if dialog.exec_():
            self.geneMatcherSettings = [getattr(dialog, item[0]) for item in dialog.items]
            self._invalidateGeneMatcher()
            if self.data is not None:
                self.updateAnnotations()

    def _genematcher(self):
        """
        Return a Future[gene.SequenceMatcher]
        """
        taxid = self.taxid_list[self.speciesIndex]

        current, matcher_f = self.__genematcher

        if taxid == current and \
                not matcher_f.cancelled():
            return matcher_f

        self._invalidateGeneMatcher()

        if taxid is None:
            self.__genematcher = (None, fulfill(gene.matcher([])))
            return self.__genematcher[1]

        matchers = [gene.GMGO, gene.GMKEGG, gene.GMNCBI, gene.GMAffy]
        matchers = [m for m, use in zip(matchers, self.geneMatcherSettings)
                    if use]

        def create():
            return gene.matcher([m(taxid) for m in matchers])

        matcher_f = self._executor.submit(create)
        self.__genematcher = (taxid, matcher_f)
        return self.__genematcher[1]

    def _nogenematching(self):
        return self.taxid is None or not any(self.geneMatcherSettings)

    def updateAnnotations(self):
        if self.data is None:
            return

        assert not self.__state & OWSetEnrichment.Initializing
        self._cancelPending()
        self._clearView()

        self.information(0)
        self.warning(0)
        self.error(0)

        if not self.genesinrows and len(self.geneAttrs) == 0:
            self.error(0, "Input data contains no columns with gene names")
            return

        self.__state = OWSetEnrichment.RunningEnrichment

        taxid = self.taxid_list[self.speciesIndex]
        self.taxid = taxid

        categories = self.selectedCategories()

        clusterGenes = self.genesFromTable(self.data)

        if self.referenceData is not None and self.useReferenceData:
            referenceGenes = self.genesFromTable(self.referenceData)
        else:
            referenceGenes = None

        self.currentAnnotatedCategories = categories

        genematcher = self._genematcher()

        self.progressBarInit()

        ## Load collections in a worker thread
        # TODO: Use cached collections if already loaded and
        # use ensure_genesetsdownloaded with progress report (OWSelectGenes)
        collections = self._executor.submit(geneset.collections, *categories)

        def refset_null():
            """Return the default background reference set"""
            col = collections.result()
            return reduce(operator.ior, (set(g.genes) for g in col), set())

        def refset_ncbi():
            """Return all NCBI gene names"""
            geneinfo = gene.NCBIGeneInfo(taxid)
            return set(geneinfo.keys())

        def namematcher():
            matcher = genematcher.result()
            match = matcher.set_targets(ref_set.result())
            match.umatch = memoize(match.umatch)
            return match

        def map_unames():
            matcher = namematcher.result()
            query = list(filter(None, map(matcher.umatch, querynames)))
            reference = list(filter(None, map(matcher.umatch, ref_set.result())))
            return query, reference

        if self._nogenematching():
            if referenceGenes is None:
                ref_set = self._executor.submit(refset_null)
            else:
                ref_set = fulfill(referenceGenes)
        else:
            if referenceGenes == None:
                ref_set = self._executor.submit(refset_ncbi)
            else:
                ref_set = fulfill(referenceGenes)

        namematcher = self._executor.submit(namematcher)
        querynames = clusterGenes

        state = types.SimpleNamespace()
        state.query_set = clusterGenes
        state.reference_set = referenceGenes
        state.namematcher = namematcher
        state.query_count = len(set(clusterGenes))
        state.reference_count = (len(set(referenceGenes))
                                 if referenceGenes is not None else None)

        state.cancelled = False

        progress = methodinvoke(self, "_setProgress", (float,))
        info = methodinvoke(self, "_setRunInfo", (str,))

        @withtraceback
        def run():
            info("Loading data")
            match = namematcher.result()
            query, reference = map_unames()
            gscollections = collections.result()

            results = []
            info("Running enrichment")
            p = 0
            for i, gset in enumerate(gscollections):
                genes = set(filter(None, map(match.umatch, gset.genes)))
                enr = set_enrichment(genes, reference, query)
                results.append((gset, enr))

                if state.cancelled:
                    raise UserInteruptException

                pnew = int(100 * i / len(gscollections))
                if pnew != p:
                    progress(pnew)
                    p = pnew
            progress(100)
            info("")
            return query, reference, results

        task = Task(function=run)
        task.resultReady.connect(self.__on_enrichment_finished)
        task.exceptionReady.connect(self.__on_enrichment_failed)
        result = self._executor.submit(task)
        state.results = result

        self.state = state
        self._updatesummary()

    def __on_enrichment_failed(self, exception):
        if not isinstance(exception, UserInteruptException):
            print("ERROR:", exception, file=sys.stderr)
            print(exception._traceback, file=sys.stderr)

        self.progressBarFinished()
        self.setStatusMessage("")
        self.__state &= ~OWSetEnrichment.RunningEnrichment

    def __on_enrichment_finished(self, results):
        assert QThread.currentThread() is self.thread()
        self.__state &= ~OWSetEnrichment.RunningEnrichment

        query, reference, results = results

        if self.annotationsChartView.model():
            self.annotationsChartView.model().clear()

        nquery = len(query)
        nref = len(reference)
        maxcount = max((len(e.query_mapped) for _, e in results),
                       default=1)
        maxrefcount = max((len(e.reference_mapped) for _, e in results),
                          default=1)
        nspaces = int(math.ceil(math.log10(maxcount or 1)))
        refspaces = int(math.ceil(math.log(maxrefcount or 1)))
        query_fmt = "%" + str(nspaces) + "s  (%.2f%%)"
        ref_fmt = "%" + str(refspaces) + "s  (%.2f%%)"

        def fmt_count(fmt, count, total):
            return fmt % (count, 100.0 * count / (total or 1))

        fmt_query_count = partial(fmt_count, query_fmt)
        fmt_ref_count = partial(fmt_count, ref_fmt)

        linkFont = QFont(self.annotationsChartView.viewOptions().font)
        linkFont.setUnderline(True)

        def item(value=None, tooltip=None, user=None):
            si = QStandardItem()
            if value is not None:
                si.setData(value, Qt.DisplayRole)
            if tooltip is not None:
                si.setData(tooltip, Qt.ToolTipRole)
            if user is not None:
                si.setData(user, Qt.UserRole)
            else:
                si.setData(value, Qt.UserRole)
            return si

        model = QStandardItemModel()
        model.setSortRole(Qt.UserRole)
        model.setHorizontalHeaderLabels(
            ["Category", "Term", "Count", "Reference count", "p-value",
             "FDR", "Enrichment"])
        for i, (gset, enrich) in enumerate(results):
            if len(enrich.query_mapped) == 0:
                continue
            nquery_mapped = len(enrich.query_mapped)
            nref_mapped = len(enrich.reference_mapped)

            row = [
                item(", ".join(gset.hierarchy)),
                item(gsname(gset), tooltip=gset.link),
                item(fmt_query_count(nquery_mapped, nquery),
                     tooltip=nquery_mapped, user=nquery_mapped),
                item(fmt_ref_count(nref_mapped, nref),
                     tooltip=nref_mapped, user=nref_mapped),
                item(fmtp(enrich.p_value), user=enrich.p_value),
                item(),  # column 5, FDR, is computed in filterAnnotationsChartView
                item(enrich.enrichment_score,
                     tooltip="%.3f" % enrich.enrichment_score,
                     user=enrich.enrichment_score)
            ]
            row[0].geneset = gset
            row[0].enrichment = enrich
            row[1].setData(gset.link, gui.LinkRole)
            row[1].setFont(linkFont)
            row[1].setForeground(QColor(Qt.blue))

            model.appendRow(row)

        self.annotationsChartView.setModel(model)
        self.annotationsChartView.selectionModel().selectionChanged.connect(
            self.commit
        )

        if not model.rowCount():
            self.warning(0, "No enriched sets found.")
        else:
            self.warning(0)

        allnames = set(gsname(geneset)
                       for geneset, (count, _, _, _) in results if count)

        allnames |= reduce(operator.ior,
                           (set(word_split(name)) for name in allnames),
                           set())

        self.filterCompleter.setModel(None)
        self.completerModel = QStringListModel(sorted(allnames))
        self.filterCompleter.setModel(self.completerModel)

        if results:
            max_score = max((e.enrichment_score for _, e in results
                             if np.isfinite(e.enrichment_score)),
                            default=1)

            self.annotationsChartView.setItemDelegateForColumn(
                6, BarItemDelegate(self, scale=(0.0, max_score))
            )

        self.annotationsChartView.setItemDelegateForColumn(
            1, gui.LinkStyledItemDelegate(self.annotationsChartView)
        )

        header = self.annotationsChartView.header()
        for i in range(model.columnCount()):
            sh = self.annotationsChartView.sizeHintForColumn(i)
            sh = max(sh, header.sectionSizeHint(i))
            self.annotationsChartView.setColumnWidth(i, max(min(sh, 300), 30))
#             self.annotationsChartView.resizeColumnToContents(i)

        self.filterAnnotationsChartView()

        self.progressBarFinished()
        self.setStatusMessage("")

    def _updatesummary(self):
        state = self.state
        if state is None:
            self.error(0,)
            self.warning(0)
            self.infoBox.setText("No data on input.\n")
            return

        text = "{.query_count} unique names on input\n".format(state)

        if state.results.done() and not state.results.exception():
            mapped, _, _ = state.results.result()
            ratio_mapped = (len(mapped) / state.query_count
                            if state.query_count else 0)
            text += ("%i (%.1f%%) gene names matched" %
                     (len(mapped), 100.0 * ratio_mapped))
        elif not state.results.done():
            text += "..."
        else:
            text += "<Error {}>".format(str(state.results.exception()))
        self.infoBox.setText(text)

        # TODO: warn on no enriched sets found (i.e no query genes
        # mapped to any set)

    def filterAnnotationsChartView(self, filterString=""):
        if self.__state & OWSetEnrichment.RunningEnrichment:
            return

        # TODO: Move filtering to a filter proxy model
        # TODO: Re-enable string search

        categories = set(", ".join(cat)
                         for cat, _ in self.selectedCategories())

#         filterString = str(self.filterLineEdit.text()).lower()

        model = self.annotationsChartView.model()

        def ishidden(index):
            # Is item at index (row) hidden
            item = model.item(index)
            item_cat = item.data(Qt.DisplayRole)
            return item_cat not in categories

        hidemask = [ishidden(i) for i in range(model.rowCount())]

        # compute FDR according the selected categories
        pvals = [model.item(i, 4).data(Qt.UserRole)
                 for i, hidden in enumerate(hidemask) if not hidden]
        fdrs = utils.stats.FDR(pvals)

        # update FDR for the selected collections and apply filtering rules
        itemsHidden = []
        fdriter = iter(fdrs)
        for index, hidden in enumerate(hidemask):
            if not hidden:
                fdr = next(fdriter)
                pval = model.index(index, 4).data(Qt.UserRole)
                count = model.index(index, 2).data(Qt.ToolTipRole)

                hidden = (self.useMinCountFilter and count < self.minClusterCount) or \
                         (self.useMaxPValFilter and pval > self.maxPValue) or \
                         (self.useMaxFDRFilter and fdr > self.maxFDR)

                if not hidden:
                    fdr_item = model.item(index, 5)
                    fdr_item.setData(fmtpdet(fdr), Qt.ToolTipRole)
                    fdr_item.setData(fmtp(fdr), Qt.DisplayRole)
                    fdr_item.setData(fdr, Qt.UserRole)

            self.annotationsChartView.setRowHidden(
                index, QModelIndex(), hidden)

            itemsHidden.append(hidden)

        if model.rowCount() and all(itemsHidden):
            self.information(0, "All sets were filtered out.")
        else:
            self.information(0)

        self._updatesummary()

    @Slot(float)
    def _setProgress(self, value):
        assert QThread.currentThread() is self.thread()
        self.progressBarSet(value, processEvents=None)

    @Slot(str)
    def _setRunInfo(self, text):
        self.setStatusMessage(text)

    def commit(self):
        if self.data is None or \
                self.__state & OWSetEnrichment.RunningEnrichment:
            return

        model = self.annotationsChartView.model()
        rows = self.annotationsChartView.selectionModel().selectedRows(0)
        selected = [model.item(index.row(), 0) for index in rows]
        mapped = reduce(operator.ior,
                        (set(item.enrichment.query_mapped)
                         for item in selected),
                        set())
        assert self.state.namematcher.done()
        matcher = self.state.namematcher.result()

        axis = 1 if self.genesinrows else 0
        if axis == 1:
            mapped = [attr for attr in self.data.domain.attributes
                      if matcher.umatch(attr.name) in mapped]

            newdomain = Orange.data.Domain(
                mapped, self.data.domain.class_vars, self.data.domain.metas)
            data = self.data.from_table(newdomain, self.data)
        else:
            geneattr = self.geneAttrs[self.geneattr]
            selected = [i for i, ex in enumerate(self.data)
                        if matcher.umatch(str(ex[geneattr])) in mapped]
            data = self.data[selected]
        self.send("Data subset", data)

    def onDeleteWidget(self):
        if self.state is not None:
            self._cancelPending()
            self.state = None
        self._executor.shutdown(wait=False)


class UserInteruptException(Exception):
    pass


from collections import namedtuple

enrichment_res = namedtuple(
    "enrichment_result",
    ["query_mapped",      #:: list
     "reference_mapped",  #:: list
     "p_value",           #:: float
     "enrichment_score"   #:: float
     ]
)


def set_enrichment(target, reference, query,
                   prob=utils.stats.Hypergeometric()):
    """
    :param set query: query set
    :param set target: target set
    :param set reference: the reference set

    """
    assert len(reference) > 0
    query_mapped = target.intersection(query)
    reference_mapped = target.intersection(reference)

    query_p = len(query_mapped) / len(query) if query else np.nan
    ref_p = len(reference_mapped) / len(reference) if reference else np.nan
    enrichment = query_p / ref_p if ref_p else np.nan

    return enrichment_res(
        list(query_mapped), list(reference_mapped),
        prob.p_value(len(query_mapped), len(reference),
                     len(reference_mapped), len(query)),
        enrichment
    )


if __name__ == "__main__":
    app = QApplication(sys.argv)
    w = OWSetEnrichment()
#     data = Orange.data.Table("yeast-class-RPR.tab")
    data = Orange.data.Table("brown-selected")
    w.setData(data)
#     w.setReference(data)
    w.handleNewSignals()
    w.show()
    app.exec_()
    w.saveSettings()
    w.onDeleteWidget()
