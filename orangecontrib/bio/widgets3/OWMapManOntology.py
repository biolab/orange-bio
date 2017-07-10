"""
Enrichment analysis on MapMan ontologies
----------------------------------------
"""
import sys
import os
import csv
import warnings
import enum
import io
import gzip

import tempfile
import numbers

import pickle
import logging
import itertools

import urllib.parse

from collections import defaultdict, deque, namedtuple, OrderedDict
from itertools import chain
from xml.sax.saxutils import escape

from types import SimpleNamespace as namespace

from typing import (
    Dict, List, Tuple, Optional, Union, Any, Iterable, Sequence, Callable,
    Generic, TypeVar, Hashable
)

import numpy

from AnyQt.QtCore import (
    Qt, QSize, QModelIndex, QAbstractProxyModel,
    QSortFilterProxyModel, QItemSelection, QItemSelectionRange,
    QItemSelectionModel, QEvent
)
from AnyQt.QtGui import QStandardItemModel, QStandardItem, QPalette
from AnyQt.QtWidgets import (
    QWidget, QTreeView, QComboBox, QApplication, QSplitter,
    QRadioButton, QButtonGroup, QProgressBar, QStackedLayout, QFormLayout,
    QStyledItemDelegate, QHBoxLayout, QSizePolicy,
    QStyle, QStylePainter, QStyleOptionComboBox, QToolTip
)
from AnyQt.QtCore import pyqtSlot as Slot

import Orange.data
from Orange.widgets import widget, gui, settings

from orangecontrib.bio import ontology
from orangecontrib.bio import gene
from orangecontrib.bio.utils import stats, environ

from orangecontrib.bio.widgets3.OWGEODatasets import retrieve_url

from orangecontrib.bio.widgets3.utils import disconnected, group_ranges
from orangecontrib.bio.widgets3.utils.data import append_columns
from orangecontrib.bio.widgets3.OWFeatureSelection import copy_variable

def download_file(url, targetpath, progress=None):
    # type: (str, str, Optional[Callable[[int, int], None]]) -> None
    """
    Download the contents specified by `url` into a directory `dirpath`.

    The contents will be writen into a unique temporary file within the
    same `os.path.dirname(targetpath)` and moved into place once the download
    completes.

    Parameters
    ----------
    url : str
        Url naming a resource to retrieve
    targetpath : str
        A local filesystem path naming a file where the contents will be
        stored.
    progress : Optional[Callable[[int, int], None]]
        Optional progress report callback. Will be called periodically
        with `(transferred, total)` bytes count. `total` can be `-1` if
        the total contents size cannot be determined beforehand.
    """
    dirpath, basename = os.path.split(targetpath)
    temp = tempfile.NamedTemporaryFile(
       prefix=basename, suffix=".part", dir=dirpath, delete=False)
    os.chmod(temp.name, 0o644)
    try:
        retrieve_url(url, temp, progress=progress)
    except BaseException:
        # cleanup temp file in case of an error
        try:
            temp.close()
            os.remove(temp.name)
        except OSError:
            # suppress 'expected' errors during cleanup, other errors should
            # propagate replacing err
            pass
        raise  # re-raise the error
    else:
        # Move the temporary into place (atomically replacing an existing file
        # if present)
        temp.close()
        os.replace(temp.name, targetpath)


def default_datadir():
    return os.path.join(environ.buffer_dir, "orangecontrib", "bio", "data")


class Resource:
    # Home page of the resource origin
    HOME = None
    # Resource version/release number/date
    VERSION = None
    # Final filename (basename) on local disk
    FILENAME = None
    # Url from which the resource can be downloaded
    URL = None
    # License under which the resource is provided
    LICENSE = None
    # Local cache directory relative to `datadir` (parameter)
    CACHEPATH = None

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    @classmethod
    def resolve(cls, **kwargs):
        home = cls.HOME.format(**kwargs)
        version = cls.VERSION.format(**kwargs)
        filename = cls.FILENAME.format(VERSION=version, **kwargs)
        url = cls.URL.format(FILENAME=urllib.parse.quote(filename),
                             VERSION=version, **kwargs)
        license = cls.LICENSE.format(**kwargs)
        cachepath = cls.CACHEPATH.format(FILENAME=filename, VERSION=version,
                                         **kwargs)
        return namespace(home=home, version=version, filename=filename,
                         url=url, license=license, cachepath=cachepath)

    def fetch(self, datadir=None, progress=None):
        if datadir is None:
            datadir = default_datadir()
        ns = self.resolve(**self.__dict__)
        cachepath = os.path.join(datadir, ns.cachepath)
        os.makedirs(cachepath, exist_ok=True)
        url = ns.url
        localpath = self.localpath(datadir=datadir)
        if not os.path.exists(localpath):
            download_file(url, localpath, progress=progress)

    def localpath(self, datadir=None):
        if datadir is None:
            datadir = default_datadir()
        ns = self.resolve(**self.__dict__)
        return os.path.join(datadir, ns.cachepath, ns.filename)

    def exists(self, datadir=None):
        return os.path.exists(self.localpath(datadir))


class GeneOntologyResource(Resource):
    HOME = "https://geneontology.org"
    VERSION = "2016-12-01"  # monthly export date
    FILENAME = "gene_ontology_edit.obo.gz"
    URL = "ftp://ftp.geneontology.org/go/ontology-archive/" \
          "gene_ontology_edit.obo.{VERSION}.gz"
    # VERSION = "4.2302"  # CVS revision
    # FILENAME = "gene_ontology_edit.obo"
    # URL = ("http://cvsweb.geneontology.org/cgi-bin/cvsweb.cgi/~checkout~/go/"
    #        "ontology/gene_ontology.obo?rev={VERSION};content-type=text%2Fplain")
    LICENSE = "CC BY 3.0"
    CACHEPATH = "geneontology/{VERSION}/"


class OntologyResource(Resource):
    HOME = "http://www.gomapman.org"
    VERSION = "2017-03-14"  # export date
    FILENAME = "ontology.obo"
    URL = "http://www.gomapman.org/export/{VERSION}/OBO/{FILENAME}"
    LICENSE = "CC BY-NC-SA 3.0"
    CACHEPATH = "gomapman/{VERSION}/OBO/"


class MappingResource(Resource):
    HOME = "http://www.gomapman.org"
    VERSION = "2017-03-14"  # export date
    FILENAME = "{organism.orgcode}_{mapping.identifiertag}_{VERSION}_mapping.txt.gz"
    URL = "http://www.gomapman.org/export/{VERSION}/mapman/{FILENAME}"
    LICENSE = "CC BY-NC-SA 3.0"
    CACHEPATH = "gomapman/{VERSION}/mapman/"

    Organism = namedtuple(
        "Organism", [
            "name",
            "ncbi_taxid",
            "common_name",
            "orgcode",
        ]
    )

    Mapping = namedtuple(
        "Mapping", [
            "organism",
            "identifiertag",
            "identifiername",
            "desc",
        ]
    )

    Arabidopsis = Organism(
        "Arabidopsis thaliana", "3702", "Thale cress", "ath")
    Potato = Organism(
        "Solanum tuberosum", "4113", "Potato", "stu")
    Rice = Organism(
        "Oryza sativa", "4530", "Rice", "osa")
    Tomato = Organism(
        "Solanum lycopersicum", "4081", "Tomato", "sly")
    Tobacco = Organism(
        "Nicotiana tabacum", "4097", "Common tobacco", "nta")
    Beet = Organism(
        "Beta vulgaris subsp. vulgaris", "3555", "Beet", "bvu")
    CacaoTree = Organism(
        "Theobroma cacao", "3641", "Cacao", "tca")

    Mappings = [
        # Arabidopsis thaliana
        Mapping(Arabidopsis, "Araport11", "Locus ID (AGI)",
                "Arabidopsis gene loci number (AGI)"),
        # Potato
        Mapping(Potato, "stNIB-v1", "stNIB-v1",
                "Combining PGSC Gene Model v3.4, ITAG Potato Gene Model v.1, "
                "POCI and StGIv13"),
        Mapping(Potato, "PGSC_gene", "PGSC gene",
                "Potato Genome Sequencing Consortium"),
        Mapping(Potato, "PGSC_protein", "PGSC protein",
                "Potato Genome Sequencing Consortium"),
        Mapping(Potato, "PGSC_transcript", "PGSC transcript",
                "Potato Genome Sequencing Consortium"),
        Mapping(Potato, "ITAG", "ITAG", ""),
        Mapping(Potato, "Agilent_4x44k", "Agilent_4x44k", ""),
        # Rice
        Mapping(Rice, "MSU_v7", "Locus ID (MSU)",
                "MSU Rice Genome Annotation Project"),
        # Tomato
        Mapping(Tomato, "SL2.40_ITAG2.3",  "Locus ID (SGN)",
                "International Tomato Annotation Group"),
        # Tobacco
        Mapping(Tobacco, "ntaUG17", "NCBI UniGene", ""),
        Mapping(Tobacco, "Affymetrix_A-AFFY-135", "", ""),
        Mapping(Tobacco, "Agilent_4x44k", "", ""),
        # Beet
        Mapping(Beet, "RefBeet1.1", "RefBeet", ""),
        Mapping(Beet, "Refbeet1.1_transcript", "RefBeet Transcript", ""),
        # Cacao Tree
        Mapping(CacaoTree, "Phytozome9.1", "Phytozome 9.1", ""),
        Mapping(CacaoTree, "Phytozome9.1_gene-PLAZA", "PLAZA", ""),
        Mapping(CacaoTree, "Phytozome9.1_transcript", "Phytozome Transcript", ""),
    ]

    def __init__(self, version=VERSION, mapping=Mappings[0], **kwargs):
        super().__init__(**kwargs)
        self.version = version
        self.mapping = mapping
        self.organism = mapping.organism
        self.identifiertag = mapping.identifiertag

    def resolve(self, **kwargs):
        kwargs["organism"] = self.mapping.organism
        kwargs["mapping"] = self.mapping
        return super().resolve(**kwargs)

    def open(self, datadir=None):
        if datadir is None:
            datadir = default_datadir()
        path = self.localpath(datadir=datadir)
        basename = os.path.basename(path)
        return gzip.open(self.localpath(datadir=datadir))


def cache_load_ontology(path):
    """
    Load/parse an .obo file into an OBOOntology instance caching the results
    in a pickle file next to path for reuse (faster future loads).
    """
    cacheversion = 0x1
    with open(path, "rb") as f:
        stat = os.fstat(f.fileno())
        # try to find/load a pickled cached file in the same dir
        try:
            with open(path + ".cache.pickle", "rb") as fcache:
                version = pickle.load(fcache)
                if version == cacheversion:
                    fingerprint = pickle.load(fcache)
                    if fingerprint == \
                            (stat.st_mtime_ns, stat.st_size, stat.st_ino):
                        return pickle.load(fcache)
        except (OSError, pickle.UnpicklingError):
            pass
        except Exception:
            log = logging.getLogger(__name__)
            log.exception("cache_load_ontology: Could not load cached %s",
                          path, exc_info=True)

        if path.endswith(".gz"):
            f = gzip.GzipFile(fileobj=f)

        ftext = io.TextIOWrapper(f, encoding="utf-8")
        ont = ontology.OBOOntology(ftext)
        ftext.detach()
        # TODO: atomic write/replace
        with open(path + ".cache.pickle", "wb") as fcache:
            pickle.dump(cacheversion, fcache)
            pickle.dump((stat.st_mtime_ns, stat.st_size, stat.st_ino), fcache)
            pickle.dump(ont, fcache)
        return ont


class RealDelegate(QStyledItemDelegate):
    """
    An Item delegate for displaying numerical columns
    """
    def __init__(self, parent=None, precision=4, **kwargs):
        super().__init__(parent, **kwargs)
        self.precision = precision

    def displayText(self, value, locale):
        if isinstance(value, numbers.Integral):
            return locale.toString(int(value))
        elif isinstance(value, numbers.Real):
            return locale.toString(float(value), "g", self.precision)
        else:
            return super().displayText(value, locale)

    def initStyleOption(self, option, index):
        super().initStyleOption(option, index)
        align = index.data(Qt.TextAlignmentRole)
        data = index.data(Qt.DisplayRole)
        if align is None and isinstance(data, numbers.Real):
            # Right align if the model does not specify otherwise
            option.displayAlignment = Qt.AlignRight | Qt.AlignVCenter


if not hasattr(QComboBox, "currentData"):
    class QComboBox(QComboBox):
        # new in Qt 5.2
        def currentData(self, role=Qt.UserRole):
            return self.itemData(self.currentIndex(), role)


class ComboBox(QComboBox):
    """
    Combo box, slightly more suitable for (faked) hierarchical display.
    """
    # Override for DisplayRole as displayed in the combo box for the current
    # item.
    CurrentDisplayRole = Qt.UserRole + 666

    def paintEvent(self, event):
        painter = QStylePainter(self)
        painter.setPen(self.palette().color(QPalette.Text))
        option = QStyleOptionComboBox()
        self.initStyleOption(option)
        painter.drawComplexControl(QStyle.CC_ComboBox, option)
        painter.drawControl(QStyle.CE_ComboBoxLabel, option)

    def initStyleOption(self, option):
        super().initStyleOption(option)
        text = self.currentData(ComboBox.CurrentDisplayRole)
        if text is not None:
            text = str(text)
        option.currentText = text or self.currentText()

    def event(self, event):
        if event.type() == QEvent.ToolTip and self.toolTip() == "":
            text = self.currentData(ComboBox.CurrentDisplayRole)
            if text is not None:
                text = str(text)
                QToolTip.showText(event.globalPos(), text, self)
                return True
        return super().event(event)


class DomainRole(enum.IntEnum):
    """
    An enum describing the position in the domain where columns should be
    placed.
    """
    Attribute, Class, Meta = 0, 1, 2

DomainRole.DisplayNames = {
    DomainRole.Attribute: "Attribute",
    DomainRole.Class: "Class variable",
    DomainRole.Meta: "Meta variable",
}


class OWMapManEnrichment(widget.OWWidget):
    name = "GoMapMan Ontology"
    icon = "icons/GoMapMan.svg"

    priority = 2020

    inputs = [
        widget.InputSignal(
            "Query", Orange.data.Table, "set_query_data", widget.Default,
            id="query-data",
            doc="Input data with query identifiers"),
        widget.InputSignal(
            "Background", Orange.data.Table, "set_reference_data",
            id="reference-data",
            doc="Optional data with reference ids (defines the background/"
                "prior distribution)"),
    ]
    outputs = [
        widget.OutputSignal(
            "Selected data", Orange.data.Table, widget.Default,
            id="selected-data"),
        widget.OutputSignal(
            "Data", Orange.data.Table,
            id="annotated-data",
            doc="The input query data with (selected) term membership "
                "indicator columns"
        )
    ]
    settingsHandler = settings.DomainContextHandler()

    selected_mapping_key = settings.Setting(None)    # type: Optional[str]
    selected_query_source = settings.ContextSetting(None, exclude_metas=False) \
        # type: Optional[Union[str, Orange.data.Variable]]
    use_reference_input = settings.Setting(True)     # type: bool

    count_filter_value = settings.Setting(3)         # type: int
    count_filter_enabled = settings.Setting(True)    # type: bool
    pvalue_filter_value = settings.Setting(0.05)     # type: float
    pvalue_filter_enabled = settings.Setting(False)  # type: bool
    fdr_filter_value = settings.Setting(0.05)        # type: float
    fdr_filter_enabled = settings.Setting(True)      # type: bool

    splitter_state = settings.Setting(b'')           # type: bytes
    header_states = settings.Setting((b'', b''))     # type: Tuple[bytes, bytes]
    append_term_indicators = settings.Setting(True)  # type: bool
    term_indicator_role = settings.Setting(
        int(DomainRole.Meta))                        # type: int
    auto_commit = settings.Setting(False)            # type: bool

    def __init__(self):
        super().__init__()
        #: Loaded ontology
        self.ontology = None         # type: ontology.OBOOntology
        #: The input query data set
        self.query_data = None       # type: Orange.data.Table
        #: The input reference data set
        self.reference_data = None   # type: Orange.data.Table

        #: The current list of query ids (derived from query_data)
        self.query_list = None       # type: List[str]
        #: The current list of reference ids (derived from reference_data
        #: or default full mapping id set)
        self.reference_list = None   # type: List[str]

        self.__data = namespace(
            query_source=None,
            query_source_vector=None,
            queryids=None,
            refids=None,
            ontology=None,
            mapping_key=None,
            mapping=None,
            matcher=None,
            queryids_resolved=None,
            refids_resolved=None,
            matcher_results=None,
            results=None,
        )
        # GUI
        box = gui.widgetBox(self.controlArea, "Info", addSpace=False)
        self.infolabel = gui.widgetLabel(box, "No data on input\n")
        self.controlArea.layout().setSpacing(-1)
        pbstack = QStackedLayout()
        pb = QProgressBar()
        pbstack.addWidget(QWidget(box))  # empty placeholder
        pbstack.addWidget(pb)
        pbstack.setCurrentIndex(0)
        box.layout().addLayout(pbstack)
        self.pb = namespace(stack=pbstack, widget=pb)

        # Drop down hierarchical menu for selection of id -> ontology mapping
        box = gui.widgetBox(self.controlArea, "Mapping", addSpace=False)

        self.mapping_combo = combo = ComboBox(
            objectName="mapping-combo-box",
            minimumContentsLength=16,
            sizeAdjustPolicy=QComboBox.AdjustToMinimumContentsLength,
        )
        model = combo.model()
        for org, orgmapping in itertools.groupby(
                MappingResource.Mappings,
                key=lambda m: m.organism):
            combo.model()
            item = QStandardItem("{0.name} ({0.common_name})".format(org))
            item.setToolTip("{0.name} (NCBI-TAXID:{0.ncbi_taxid}".format(org))
            item.setEnabled(False)
            item.setSelectable(False)
            model.appendRow([item])

            for mapping in orgmapping:  # type: MappingResource.Mapping
                m = MappingResource(mapping=mapping)
                ns = m.resolve(**m.__dict__)
                item = QStandardItem("  " + mapping.identifiertag)
                item.setToolTip(mapping.desc)  # need better tooltip
                # Organism | entity id db (human readable name) + version (tag)
                item.setData(m, Qt.UserRole)
                key = str(ns.filename)
                item.setData(key, Qt.UserRole + 1)
                item.setData("{0.name} / {1.identifiername}"
                             .format(org, mapping),
                             ComboBox.CurrentDisplayRole)
                model.appendRow(item)

        index = combo.findData(self.selected_mapping_key, Qt.UserRole + 1,
                               Qt.MatchExactly)
        if index != -1 and model.item(index).isSelectable():
            combo.setCurrentIndex(index)
        else:
            combo.setCurrentIndex(-1)
        combo.activated.connect(self.__mapping_activated)

        if combo.currentIndex() == -1:
            self.selected_mapping_key = ""
        self.__data.mapping_key = self.selected_mapping_key or None

        box.layout().addWidget(combo)
        box = gui.widgetBox(self.controlArea, "Entity identifiers",
                            addSpace=False)
        self.query_combo = combo = QComboBox(
            objectName="entity-identifier-combo",
            minimumContentsLength=12,
            sizeAdjustPolicy=QComboBox.AdjustToMinimumContentsLength,
            toolTip="Select the identifier source column/row.",
        )
        self.query_combo.activated.connect(self.__query_id_activated)
        box.layout().addWidget(combo)

        box = gui.widgetBox(self.controlArea, "Reference", addSpace=False)

        bbox = QButtonGroup(box, exclusive=True)
        b1 = QRadioButton("All entities", box)
        b2 = QRadioButton("Reference set (none provided)", box)
        bbox.addButton(b1, 1)
        bbox.addButton(b2, 2)
        if self.use_reference_input:
            b2.setChecked(True)
        else:
            b1.setChecked(True)

        @bbox.buttonClicked.connect
        def on_ref_toggled():
            self.use_reference_input = b2.isChecked()
            self.__data.refids = None
            self.__data.queryids_resolved = None
            self.__update()

        self.refbuttons = namespace(
            allbutton=b1,
            inputbutton=b2
        )

        box.layout().addWidget(b1)
        box.layout().addWidget(b2)
        self.controlArea.layout().addStretch(10)

        box = gui.widgetBox(self.controlArea, "Output", addSpace=False)
        append_check = gui.checkBox(box, self, "append_term_indicators",
                                    "Append term indicator columns")
        append_check.setAttribute(Qt.WA_LayoutUsesWidgetRect, True)

        form = QFormLayout(
            formAlignment=Qt.AlignLeft,
            labelAlignment=Qt.AlignLeft,
            fieldGrowthPolicy=QFormLayout.AllNonFixedFieldsGrow,
        )
        box.layout().addLayout(form)
        append_check.ensurePolished()
        form.setContentsMargins(
            gui.checkButtonOffsetHint(append_check), 0, 6, 6)
        position_cb = gui.comboBox(box, self, "term_indicator_role")
        position_cb.setEnabled(bool(self.append_term_indicators))

        for role in DomainRole:
            position_cb.addItem(DomainRole.DisplayNames[role])

        position_cb.activated.connect(self._invalidate_output)
        form.addRow("Place as:", position_cb)
        append_check.toggled[bool].connect(position_cb.setEnabled)

        gui.auto_commit(box, self, "auto_commit", "Commit", box=False,
                        commit=self.commit)

        # Main area views
        #################

        box = gui.widgetBox(self.mainArea, orientation=QHBoxLayout())
        gui.widgetLabel(box, "Filters:")
        box.layout().addSpacing(8)

        gui.spin(box, self, "count_filter_value", label="Entities",
                 minv=1, maxv=100, checked="count_filter_enabled",
                 tooltip="Minimum mapped ID count",
                 callbackOnReturn=True,
                 callback=self.__invalidate_filter,
                 checkCallback=self.__invalidate_filter)

        gui.spin(box, self, "pvalue_filter_value", label="P-Value",
                 minv=0.0001, maxv=1.0, spinType=float, decimals=4, step=0.0001,
                 checked="pvalue_filter_enabled",
                 tooltip="Maximum (unadjusted) p-value",
                 callbackOnReturn=True,
                 callback=self.__invalidate_filter,
                 checkCallback=self.__invalidate_filter
                 )
        gui.spin(box, self, "fdr_filter_value", label="FDR",
                 minv=0.0001, maxv=1.0, spinType=float, decimals=4, step=0.0001,
                 checked="fdr_filter_enabled",
                 tooltip="Maximum FDR adjusted p-value",
                 callback=self.__invalidate_filter,
                 callbackOnReturn=True,
                 checkCallback=self.__invalidate_filter
                 )

        box.layout().addStretch(10)
        box.setSizePolicy(QSizePolicy.MinimumExpanding,
                          QSizePolicy.Fixed)

        self.splitter = QSplitter(orientation=Qt.Vertical)
        self.mainArea.layout().addWidget(self.splitter)

        # tree display of the ontology
        self.treeview = QTreeView(
            sortingEnabled=True,
            alternatingRowColors=True,
            selectionMode=QTreeView.ExtendedSelection,
            editTriggers=QTreeView.NoEditTriggers,
        )

        self.treeview.setModel(OntologyEnrichmentResultsModel(parent=self))
        self.treeview.selectionModel().selectionChanged.connect(
            self.__tree_selection_changed, Qt.UniqueConnection)

        h = self.treeview.header()
        h.resizeSection(0, 300)
        d = gui.VisibleHeaderSectionContextEventFilter(self)
        h.installEventFilter(d)

        # Non-hierarchical table view of the results
        self.tableview = QTreeView(
            sortingEnabled=True,
            alternatingRowColors=True,
            rootIsDecorated=False,
            selectionMode=QTreeView.ExtendedSelection,
            editTriggers=QTreeView.NoEditTriggers
        )
        # Setup a default empty model to supply the header
        proxy = FilterProxyModel(parent=self)
        proxy.setSourceModel(EnrichmentTableModel(parent=proxy))
        self.tableview.setModel(proxy)
        # The table view is connected to two slots. The first syncs the
        # selection to the treeview, the second dispatches to commit method
        # TODO: Use one slot only (disconnected will change invoke order)
        self.tableview.selectionModel().selectionChanged.connect(
            self.__table_selection_changed)
        self.tableview.selectionModel().selectionChanged.connect(
            self.__selection_changed)

        h = self.tableview.header()
        h.resizeSection(0, 300)
        d = gui.VisibleHeaderSectionContextEventFilter(self)
        h.installEventFilter(d)

        self.splitter.addWidget(self.treeview)
        self.splitter.addWidget(self.tableview)
        self.splitter.setSizes([250, 100])

        # TODO: remove call here, ensure filter is updated when results change
        self.__invalidate_filter()  # __update_table_filter

        columns = EnrichmentTableModel.Column
        numericalcols = [columns.Count, columns.Reference, columns.FDR,
                         columns.Pval]

        tabledelegate = RealDelegate(parent=self)
        treedelegate = RealDelegate(parent=self)

        for col in numericalcols:
            self.tableview.setItemDelegateForColumn(int(col), tabledelegate)
            self.treeview.setItemDelegateForColumn(int(col), treedelegate)

        # Initial column visible state
        for header in [self.tableview.header(), self.treeview.header()]:
            header.setSectionHidden(int(columns.TermID), True)
            header.setSectionHidden(int(columns.LogEnrichment), True)
            header.setSortIndicator(int(columns.Pval), Qt.AscendingOrder)

        headers = [self.treeview.header(), self.tableview.header()]
        for header, state in zip(headers, self.header_states):
            if state:
                header.restoreState(state)

        self.splitter.restoreState(self.splitter_state)

    def clear_query_data(self):
        """
        Clear the widget state of the query_data state.
        """
        self.query_data = None
        self.__data.query_source = None
        self.__data.query_source_vector = None
        self.__data.queryids = None
        self.__data.queryids_resolved = None

        self.query_combo.clear()
        self.query_list = []
        self.reference_list = []
        with disconnected(self.tableview.selectionModel().selectionChanged,
                          self.__table_selection_changed):
            model = self.tableview.model().sourceModel()
            model.setSourceResults([])

        with disconnected(self.treeview.selectionModel().selectionChanged,
                          self.__tree_selection_changed):
            model = self.treeview.model()
            model.setGraph({})
            model.setResults([])

        self.unconditional_commit()

    def set_query_data(self, data):
        """
        Set the input query data set

        Parameters
        ----------
        data : Orange.data.Table
        """
        self.closeContext()
        self.clear_query_data()

        genevars = []
        geneannot = []
        iderror = ""
        if data is not None:
            genevars = [var for var in data.domain.metas
                        if isinstance(var, Orange.data.StringVariable)]
            geneannot = ["Use attribute names"] if data.domain.attributes else []

            if not (genevars or geneannot):
                iderror = "No string identifiers"
                data = None

        self.error(iderror)
        if genevars or geneannot:
            for item in genevars:
                self.query_combo.addItem(*gui.attributeItem(item),
                                         userData=item)
            if genevars and geneannot:
                # must not insert separator in the first position
                self.query_combo.insertSeparator(self.query_combo.count())

            for item in geneannot:
                self.query_combo.addItem(item, userData=item)

        self.query_data = data
        self.__data.results = None
        self.__data.queryids = None
        self.__data.queryids_resolved = None
        self.__data.query_source_vector = None
        self.__data.query_source = None

        def find(iterable, pred):
            # type: (Iterable[T], Callable[[T], bool]) -> int
            for i, el in enumerate(iterable):
                if pred(el):
                    return i
            return -1

        if data is not None:
            self.openContext(data)
            # restore selected query id column
            model = self.query_combo.model()
            index = find((model.item(i).data(Qt.UserRole)
                          for i in range(model.rowCount())),
                         lambda v: v == self.selected_query_source)
            if index != -1 and model.item(index).isSelectable():
                self.query_combo.setCurrentIndex(index)
            assert model.item(self.query_combo.currentIndex()).isSelectable()

        self.selected_query_source = self.query_combo.currentData()

    def set_reference_data(self, data):
        self.reference_data = data
        if self.use_reference_input:
            self.__data.results = None
            self.__data.refids = None
            self.__data.refids_resolved = None
        if data is None:
            self.refbuttons.inputbutton.setText("Reference set (none provided)")
        else:
            self.refbuttons.inputbutton.setText("Reference set")

    def handleNewSignals(self):
        if self.__data.results is None:
            self.__update()

        self.__update_infotext()
        self.unconditional_commit()

    def query_list_source(self):
        return self.query_combo.currentData()

    def __update_mapping(self):
        # update the loaded mapping
        mapping = self.mapping_combo.currentData()
        if mapping is None:
            self.__data.matcher = None
            self.__data.matcher_results = None
            self.__data.refids_resolved = None
            self.__data.queryids_resolved = None
            self.__data.mapping = None
            self.__data.mapping_key = None
            return
        assert isinstance(mapping, MappingResource)
        organism = mapping.organism

        res = MappingResource(mapping=mapping)
        res.fetch()

        with gzip.open(res.localpath()) as mcontents:
            f = io.TextIOWrapper(mcontents, encoding="utf-8")
            nodes, mapping = mapman_mapping(f)
            f.close()
            mapping = [(r.identifier, "GMM:" + r.bincode) for r in mapping]

        # id match
        if organism == MappingResource.Arabidopsis:
            matcher = gene.GMNCBI(organism.ncbi_taxid, ignore_case=False)
            matcher.set_targets({sid for sid, _ in mapping})
        else:
            matcher = None

        self.__data.matcher = matcher
        self.__data.queryids_resolved = None
        self.__data.refids_resolved = None
        self.__data.matcher_results = None
        self.__data.mapping = mapping
        self.__data.mapping_key = self.selected_mapping_key

    def __update_queryids(self):
        self.__data.queryids_resolved = None
        self.__data.queryids = None
        self.__data.queryids = None

        if self.query_data is None:
            return

        value = self.query_list_source()
        if value is None:
            return
        elif isinstance(value, Orange.data.StringVariable):
            query_vector = self.query_data[:, value].metas.ravel()
            queryids = query_vector[query_vector != value.Unknown]
            self.query_list = list(unique(queryids))
        elif isinstance(value, str):
            query_vector = numpy.array(
                [v.name for v in self.query_data.domain.attributes],
                dtype=object)
            self.query_list = list(unique(query_vector))
        else:
            assert False
        self.__data.query_source = value
        self.__data.query_source_vector = query_vector
        self.__data.queryids = self.query_list

    def __update_refids(self):
        self.__data.refids = None
        self.__data.refids_resolved = None

        if self.reference_data is None or self.query_data is None or \
                not self.use_reference_input:
            self.warning()
            return
        value = self.query_list_source()
        if isinstance(value, Orange.data.StringVariable) and \
                value not in self.reference_data.domain.metas:
            self.warning("{} is not in the reference dataset".format(value))
            return
        else:
            self.warning()

        if isinstance(value, Orange.data.StringVariable):
            refids = self.reference_data[:, value].metas
            refids = list(refids[refids != value.Unknown])
            self.reference_list = list(unique(refids))
        else:
            assert isinstance(value, str)
            refids = [v.name for v in self.reference_data.domain.attributes]
            self.reference_list = list(unique(refids))
        self.__data.refids = refids

    def __invalidate_filter(self):
        # Invalidate and update the effective enrichment results table filters
        # TODO: separate the tree/tree view filter update
        #       call __update_table_filter in __init__
        # TODO: restore selection in tree view after update
        proxy = self.tableview.model()
        model = proxy.sourceModel()
        assert isinstance(proxy, FilterProxyModel)
        assert isinstance(model, EnrichmentTableModel)
        filters = []
        subset_res = self.__data.results or []
        if self.count_filter_enabled:
            min_count = self.count_filter_value
            subset_res = (r for r in subset_res
                          if len(r.query_mapped) >= min_count)
            filters.append(
                FilterProxyModel.Filter(
                    column=int(model.Column.Count),
                    role=Qt.DisplayRole,
                    predicate=lambda val: val >= min_count)
            )
        if self.pvalue_filter_enabled:
            max_pval = self.pvalue_filter_value
            subset_res = (r for r in subset_res if r.p_value <= max_pval)
            filters.append(
                FilterProxyModel.Filter(
                    column=int(model.Column.Pval),
                    role=Qt.DisplayRole,
                    predicate=lambda val: val <= max_pval)
            )
        if self.fdr_filter_enabled:
            max_fdr = self.fdr_filter_value
            subset_res = (r for r in subset_res if r.fdr_value <= max_fdr)
            filters.append(
                FilterProxyModel.Filter(
                    column=int(model.Column.FDR),
                    role=Qt.DisplayRole,
                    predicate=lambda val: val <= max_fdr)
            )
        proxy.setFilters(filters)
        subset_res = list(subset_res)  # type: List[Result]
        self.__update_tree_view_filter()
        # Take the max 20 most significant terms and expand them in the
        # tree view
        nsig_expand = 20
        subset_res = sorted(subset_res, key=lambda r: r.p_value)
        subset = subset_res[:nsig_expand]
        subset = [r.node for r in subset]
        self.__treeview_ensure_expanded(subset)

    def __update_tree_view_filter(self):
        model = self.treeview.model()
        if self.__data.ontology is None:
            model.setGraph({})
            return

        subset_res = self.__data.results or []  # type: List[Result]
        # subset_res = self.__data.results_filterd or []  # type: List[Result]
        assert isinstance(model, OntologyEnrichmentResultsModel)
        min_count = self.count_filter_value if self.count_filter_enabled else 0
        max_pval = self.pvalue_filter_value if self.pvalue_filter_enabled else 1
        max_fdr =self.fdr_filter_value if self.fdr_filter_enabled else 1
        subset_res = [r for r in subset_res
                      if len(r.query_mapped) >= min_count and
                          r.p_value <= max_pval and
                          r.fdr_value <= max_fdr]
        graph = simple_graph(self.__data.ontology, reltypes=["is_a"])
        graph = {self.__data.ontology.term(n):
                     [self.__data.ontology.term(s) for s in succ]
                 for n, succ in graph.items()}
        termsubset = {r.node for r in subset_res}
        displayedgraph = dag_subgraph(graph, termsubset)
        selected = self.treeview.selectionModel().selectedRows(0)
        terms = [idx.data(Qt.UserRole) for idx in selected]
        model.setGraph(displayedgraph)
        terms = {t for t in terms if t in termsubset}
        with disconnected(self.treeview.selectionModel().selectionChanged,
                          self.__tree_selection_changed):
            self.__treeview_expand_and_select(list(terms), list(terms))

    def __update_enrichment(self):
        data = self.__data
        if data.mapping is None or data.ontology is None:
            self.__set_enrichment_results(None)
            return

        assert data.ontology is not None
        assert data.mapping is not None
        assert data.mapping_key is self.selected_mapping_key
        assert data.queryids is not None

        graph = simple_graph(data.ontology, reltypes=["is_a"])

        term_alt_id_mapping = {}
        for t in data.ontology.terms():
            term_alt_id_mapping.update(dict.fromkeys(t.alt_id, t.id))

        # map the ids in mapping through any term alt_id
        mapping = [(nameid, term_alt_id_mapping.get(termid, termid))
                   for nameid, termid in data.mapping]
        mapping = relation_list_to_multimap(mapping)

        queryids = data.queryids
        if data.refids is not None and self.use_reference_input:
            refids = data.refids
        else:
            refids = set(mapping.keys())

        idmapper = self.__data.matcher
        if idmapper is not None:
            def mapid(qid):
                # type: (str) -> Tuple[str, Optional[str]]
                if qid in mapping:
                    return qid, qid
                else:
                    return qid, idmapper.umatch(qid)
        else:
            def mapid(qid):
                # type: (str) -> Tuple[str, Optional[str]]
                if qid in mapping:
                    return qid, qid
                else:
                    return qid, None

        # mapping from (input) query/reference ids to (resolved) annotation ids
        inputid_to_annotationid = {
            qid: mapid(qid) for qid in set(queryids) | set(refids)
        }

        resolve_map = {inputid: mid
                       for inputid, mid in inputid_to_annotationid.values()
                       if mid is not None}
        # resolve query/reference ids to annotation ids
        queryids_resolved = {resolve_map[q]
                             for q in queryids if q in resolve_map}
        refids_resolved = {resolve_map[r]
                           for r in refids if r in resolve_map}
        # map query and reference ids onto the ontology terms
        querymap = {q: mapping.get(resolve_map.get(q, q), []) for q in queryids}
        refmap = {r: mapping.get(resolve_map.get(r, r), []) for r in refids}

        res = dag_enrichment(graph, querymap, refmap)
        res = [Result(data.ontology.term(r.node), *r[1:])
               for r in res]

        self.__data.matcher_results = inputid_to_annotationid
        self.__data.queryids_resolved = queryids_resolved
        self.__data.refids_resolved = refids_resolved
        self.__set_enrichment_results(res)

    def __set_enrichment_results(self, results):
        # type: (Optional[List[Result[ontology.Term, str]]]) -> None
        # set/clear the enrichment analysis results
        self.__data.results = results
        self.__update_enrichment_view()

    def __update_enrichment_view(self):
        results = self.__data.results
        if results is None:
            # Clear the enrichment views
            if self.tableview.model():
                self.tableview.model().sourceModel().setSourceResults([])
                assert not self.tableview.selectionModel().selectedRows(0)
            if self.treeview.model():
                model = self.treeview.model()
                model.setResults([])
            return

        model = self.tableview.model().sourceModel()
        assert isinstance(model, EnrichmentTableModel)
        model.setSourceResults(results)

        model = self.treeview.model()
        assert isinstance(model, OntologyEnrichmentResultsModel)
        model.setResults(results)
        if model.canFetchMore(QModelIndex()):
            model.fetchMore(QModelIndex())
        self.treeview.setExpanded(model.index(0, 0), True)

        self.__invalidate_filter()

    def __update_infotext(self):
        """
        Update the summary information text display
        """
        tmplt = (
            "<table>\n"
            "<tr><td>Ids</td><td>Mapped</td></tr>\n"
            "<tr><td>Query:</td><td>{} out of {}</td>\n"
            "<tr><td>Reference:</td><td>{} out of {}</td>\n"
            "</table>"
        )
        if self.query_data is not None and self.__data.mapping is not None:
            data = self.__data
            nquery = len(data.queryids or [])
            nquerymapped = len(data.queryids_resolved or [])

            if data.refids is None or not self.use_reference_input:
                nref = len({id for id, _ in data.mapping})
            else:
                nref = len(data.refids or [])
            nrefmapped = len(data.refids_resolved or [])
            text = tmplt.format(nquerymapped, nquery, nrefmapped, nref)
        elif self.query_data is not None:
            # no selected mapping
            text = "No mapping selected<br/><br/>"
        else:
            text = "No data on input<br/><br/>"

        self.infolabel.setTextFormat(Qt.RichText)
        self.infolabel.setText(text)

    def __update(self):
        if self.__data.mapping is None:
            self.__update_mapping()
        if self.__data.queryids is None:
            self.__update_queryids()
        if self.__data.refids is None:
            self.__update_refids()
        if self.__data.ontology is None:
            res = OntologyResource()
            res.fetch()
            self.__data.ontology = cache_load_ontology(res.localpath())
            self.ontology = self.__data.ontology
        if self.query_data is not None:
            self.__update_enrichment()
        self.__update_infotext()

    def __mapping_activated(self, index):
        # load the selected id -> term mapping
        key = self.mapping_combo.itemData(index, Qt.UserRole + 1)
        if self.__data.mapping_key == key:
            return

        assert key is not None
        self.selected_mapping_key = key
        # invalidate partial results and run update
        self.__data.mapping = None
        self.__data.mapping_key = self.selected_mapping_key
        self.__data.queryids = None
        self.__data.queryids_resolved = None
        self.__data.refids = None
        self.__data.refids_resolved = None
        self.__data.results = None
        self.__update()

    def __query_id_activated(self, index):
        # The selected query id column/row was changed by user.
        # Update the query ids.
        self.selected_query_source = self.query_combo.currentData()
        self.__data.queryids = None
        self.__data.queryids_resolved = None
        if self.use_reference_input:
            self.__data.refids = None
            self.__data.refids_resolved = None
        self.__data.results = None
        self.__update()

    def onDeleteWidget(self):
        self.clear_query_data()
        self.ontology = None
        self.__data = None

    @Slot()
    def __tree_selection_changed(self):
        # item selection in the tree view has changed. Sync with the
        # tableview's selection model
        selected_tree = self.treeview.selectionModel().selectedRows(0)
        terms = [index.data(Qt.UserRole) for index in selected_tree]
        assert all(isinstance(t, ontology.Term) for t in terms)
        selected = set(terms)
        table = self.tableview.model()

        rows = []
        for r in range(table.rowCount()):
            tidx = table.index(r, 0)
            tt = tidx.data(Qt.UserRole)
            if tt in selected:
                rows.append(r)

        row_ranges = group_ranges(rows)

        selection = QItemSelection()
        for start, stop in row_ranges:
            selection.append(
                QItemSelectionRange(table.index(start, 0),
                                    table.index(stop - 1, 0)))

        tableselection = self.tableview.selectionModel()
        with disconnected(tableselection.selectionChanged,
                          self.__table_selection_changed):
            tableselection.select(
                selection, (QItemSelectionModel.ClearAndSelect |
                            QItemSelectionModel.Rows))

    @Slot()
    def __table_selection_changed(self):
        selected_table = self.tableview.selectionModel().selectedRows(0)
        selected = [idx.data(Qt.UserRole) for idx in selected_table]
        assert all(isinstance(t, ontology.Term) for t in selected)
        treemodel = self.treeview.model()
        assert isinstance(treemodel, OntologyEnrichmentResultsModel)

        with disconnected(self.treeview.selectionModel().selectionChanged,
                          self.__tree_selection_changed):
            self.__treeview_expand_and_select(
                selected, selected,
                command=(QItemSelectionModel.ClearAndSelect |
                         QItemSelectionModel.Rows)
            )

    def __treeview_ensure_expanded(self, subset):
        # type: (List[ontology.Term]) -> None
        self.__treeview_expand_and_select(
            subset, [], QItemSelectionModel.NoUpdate)

    def __treeview_expand_and_select(self, expanded, selected,
                                     command=QItemSelectionModel.ClearAndSelect):
        # type: (List[ontology.Term], List[ontology.Term], QItemSelectionModel.Command) -> None
        view = self.treeview
        model = view.model()  # type: OntologyEnrichmentResultsModel
        graph = model.graph()
        invgraph = graph_invert(graph)

        expanded_subgraph = graph_select_closure(invgraph, expanded)
        expanded_terms = set(expanded_subgraph)
        selected_terms = set(selected)
        selected_indices = []  # type: List[QModelIndex]
        expanded_indices = []  # type: List[QModelIndex]

        def should_expand(index):
            assert isinstance(index.data(Qt.UserRole), ontology.Term)
            return index.data(Qt.UserRole) in expanded_terms

        def should_select(index):
            assert isinstance(index.data(Qt.UserRole), ontology.Term)
            return index.data(Qt.UserRole) in selected_terms

        def expand_and_select(view, index):
            # type: (QTreeView, QModelIndex) -> None
            model = view.model()
            if should_expand(index):
                expanded_indices.append(index)
                if model.canFetchMore(index):
                    model.fetchMore(index)

            if should_select(index):
                selected_indices.append(index)

            for i in range(model.rowCount(index)):
                expand_and_select(view, model.index(i, 0, index))

        if model.canFetchMore(QModelIndex()):
            model.fetchMore(QModelIndex())

        for row in range(model.rowCount()):
            expand_and_select(
                self.treeview, model.index(row, 0, QModelIndex()))

        for index in expanded_indices:
            view.setExpanded(index, True)

        selection = QItemSelection()
        for selected in selected_indices:
            selection.select(selected, selected)

        view.selectionModel().select(selection, command)

    @Slot()
    def __selection_changed(self):
        self._invalidate_output()

    def _invalidate_output(self):
        self.commit()

    def commit(self):
        def selected_rows(view, column=0):
            model = view.model()
            selection = view.selectionModel()
            rows = selection.selectedRows(column)
            if isinstance(model, QAbstractProxyModel):
                # map to source rows
                rows = [model.mapToSource(row) for row in rows]
            return rows

        if self.query_data is not None and self.__data.results is not None:
            results = self.__data.results  # type: List[Result]
            rows = selected_rows(self.tableview, column=0)
            terms = [idx.data(Qt.UserRole) for idx in rows]
            assert all(isinstance(term, ontology.OBOObject) for term in terms)
            terms = list(unique(terms))
            selected_termids = [t.id for t in terms]
            query_vector = self.__data.query_source_vector
            query_source = self.__data.query_source
            selected_termids_set = set(selected_termids)

            # selected term id to query ids that map to that set
            tid_2_query = {r.node.id: set(r.query_mapped) for r in results
                           if r.node.id in selected_termids_set}
            sel_qmapped = [tid_2_query[tid] for tid in selected_termids]
            indicator = [[qid in qm for qm in sel_qmapped]
                         for qid in query_vector]
            indicator = numpy.array(indicator, dtype=bool)
            assert indicator.ndim == 2
            selected_mask = numpy.logical_or.reduce(indicator, axis=1)
            indices = numpy.flatnonzero(selected_mask)

            if isinstance(query_source, Orange.data.StringVariable):
                assert len(query_vector) == len(self.query_data)
                selected_data = self.query_data[indices]
                if self.append_term_indicators:
                    term_vars = [Orange.data.DiscreteVariable(
                                     term.name, values=["No", "Yes"])
                                 for term in terms]
                    term_cols = list(zip(term_vars, indicator.T))
                else:
                    term_vars, term_cols = [], []
                select_var = Orange.data.DiscreteVariable(
                    "Selected", ["No", "Yes"])
                selected_cols = [(select_var, selected_mask)]
                # 'Selected' Column is always in metas, the 'terms'
                # positions are specified by the user.
                if not self.append_term_indicators:
                    args = {"metas": selected_cols}
                elif self.term_indicator_role == DomainRole.Attribute:
                    args = {"attributes": term_cols, "metas": selected_cols}
                elif self.term_indicator_role == DomainRole.Class:
                    args = {"class_vars": term_cols, "metas": selected_cols}
                elif self.term_indicator_role == DomainRole.Meta:
                    args = {"metas": term_cols + selected_cols}
                else:
                    assert False
                annotated_data = append_columns(self.query_data, **args)
            else:
                assert len(query_vector) == len(self.query_data.domain.attributes)
                selected_data = self.query_data.from_table(
                    Orange.data.Domain(
                        [self.query_data.domain[i] for i in indices],
                        self.query_data.domain.class_vars,
                        self.query_data.domain.metas),
                    self.query_data
                )

                if self.append_term_indicators:
                    newattrs = [copy_variable(v)
                                for v in self.query_data.domain.attributes]
                    assert len(newattrs) == len(query_vector)
                    for v, indicatorrow in zip(newattrs, indicator):
                        for term, ismember in zip(terms, indicatorrow):
                            v.attributes[term.name] = ismember  # str or not?
                    annotated_data = self.query_data.from_table(
                        Orange.data.Domain(
                            newattrs, self.query_data.domain.class_vars,
                            self.query_data.domain.metas),
                        self.query_data
                    )
                    annotated_data.X[:] = self.query_data.X[:]
                else:
                    annotated_data = self.query_data

        else:
            selected_data = None
            annotated_data = None

        self.send("Selected data", selected_data)
        self.send("Data", annotated_data)

    def sizeHint(self):
        sh = super().sizeHint()
        return sh.expandedTo(QSize(1000, 700))

    def closeEvent(self, event):
        """
        Reimplemented from OWWidget
        """
        # save view header states
        self.header_states = (
            bytes(self.treeview.header().saveState()),
            bytes(self.tableview.header().saveState()),
        )
        self.splitter_state = bytes(self.splitter.saveState())
        super().closeEvent(event)


def Term_html(term):
    # type: (ontology.Term) -> str
    """
    Summarize a term as a html fragment.
    """
    parts = [
        ("Id", term.id),
        ("Name", term.name),
        ("Definition", term.def_),
        ("Namespace", term.namespace),
        ("Synonyms", ", ".join(term.synonyms)),
        ("Comment", term.comment),
    ]

    parts = [(t, d) for t, d in parts if d]
    dttemplate = '<dt><b>{}</b></dt><dd>{}</dd>'

    def format_dt(t, d):
        return dttemplate.format(escape(t), escape(d))

    html = ["<dl>"] + list(format_dt(*p) for p in parts) + ["</dl>"]
    return "\n".join(html)


def join_elided(sep, maxlen, values, elidetemplate=" ..."):
    """
    Parameters
    ----------
    sep : str
        Join string separator
    maxlen : int
        Maximum (character) string length
    values : List[str]
        List of strings to join
    elidetemplate : str
        Elide template; string appended at the end to indicate elision took
        place.

    Returns
    -------
    string : str

    Example
    -------
    >>> join_elided(", ", 9, ["1", "2", "3", "4"], ", ...")
    '1, 2, ...'
    >>> join_elided(", ", 24, ["Alpha", "Beta", "Charlie", "Delta", "Echo"],
    ...             " (and {} more)")
    'Alpha, Beta (and 3 more)'
    >>> join_elided(", ", 0, [], "...")
    ''
    """
    def generate(sep, elidetemplate, values):
        count = len(values)
        length = 0
        parts = []
        for i, val in enumerate(values):
            if count - i > 1:
                elidesuffix = elidetemplate.format(count - i - 1)
            else:
                elidesuffix = ""
            length += len(val) + (len(sep) if parts else 0)
            parts.append(val)
            yield i, itertools.islice(parts, i + 1), length, elidesuffix

    best = None
    for i, parts, length, elide in generate(sep, elidetemplate, values):
        if length > maxlen:
            if best is None:
                best = sep.join(parts) + elide
            return best
        fulllen = length + len(elide)
        if fulllen <= maxlen or best is None:
            best = sep.join(parts) + elide
    else:
        return best or ""


# for memory conservation purposes (long rarely accessed tooltips)
class _Item(QStandardItem):
    """
    Parameters
    ----------
    data : Dict[Qt.ItemDataRole, Callable[[], Any]]
        A callback dictionary dispatching queries for particular roles.
    """
    def __init__(self, *args, data={}, **kwargs):
        super().__init__(*args, **kwargs)
        self.__data = data  # type: Dict[int, Callable[[], Any]]

    def data(self, role=Qt.DisplayRole):
        if role in self.__data:
            # TODO: What to do in case of errors?
            return self.__data[role]()
        else:
            return super().data(role)


class EnrichmentTableModel(QStandardItemModel):
    class Column(enum.IntEnum):
        Name, TermID, Count, Reference, FDR, Pval, Items, \
            Enrichment, LogEnrichment = range(9)

    HeaderLabels = [
        (Column.Name, "Name"),
        (Column.TermID, "Term ID"),
        (Column.Count, "Count"),
        (Column.Reference, "Reference"),
        (Column.FDR, "FDR(P value)"),
        (Column.Pval, "P value"),
        (Column.Items, "Genes"),
        (Column.Enrichment, "Enrichment"),
        (Column.LogEnrichment, "Log2(Fold Enrichment)")
    ]

    def __init__(self, parent=None, **kwargs):
        super().__init__(parent, **kwargs)
        for col, text in EnrichmentTableModel.HeaderLabels:
            item = QStandardItem(text)
            self.setHorizontalHeaderItem(int(col), item)

    def setSourceResults(self, results):
        """
        Parameters
        ----------
        results : List[Result]
        """
        self.setRowCount(0)  # reset/clear the model
        for res in results:  # type: Result[ontology.Term, str]
            term = res.node
            html = Term_html(term)

            termitem = QStandardItem(term.name)
            termitem.setData(term, Qt.UserRole)
            termitem.setToolTip(html)
            termid = QStandardItem(term.id)
            termid.setToolTip(html)
            count = QStandardItem()
            count.setData(len(res.query_mapped), Qt.DisplayRole)
            ref = QStandardItem()
            ref.setData(len(res.reference_mapped), Qt.DisplayRole)
            pval = QStandardItem()
            pval.setData(res.p_value, Qt.DisplayRole)
            fdr = QStandardItem()
            fdr.setData(res.fdr_value, Qt.DisplayRole)
            querylist = list(res.query_mapped)

            def tooltip(elts):
                return "<div>" + ", ".join(map(escape, elts)) + "</div>"

            genelist = _Item(
                data={
                    Qt.DisplayRole: lambda ql=querylist:
                        join_elided(", ", 80, ql, " (and {} more)"),
                    Qt.ToolTipRole: lambda ql=querylist: tooltip(ql),
                    Qt.UserRole: lambda ql=querylist: ql,
                }
            )
            enrich = QStandardItem()
            enrich.setData(res.enrichment_score, Qt.DisplayRole)
            logenrich = QStandardItem()
            logenrich.setData(float(numpy.log2(res.enrichment_score)),
                              Qt.DisplayRole)
            row = [termitem, termid, count, ref, fdr, pval, genelist,
                   enrich, logenrich]

            self.appendRow(row)


class LazyItemModel(QStandardItemModel):
    class LazyItem(QStandardItem):
        pass

    def __init__(self, parent=None, **kwargs):
        super().__init__(parent, **kwargs)

    def hasChildren(self, parent=QModelIndex()):
        if not parent.isValid():
            return self.rowCountHint() > 0
        item = self.itemFromIndex(parent)
        if not isinstance(item, LazyItemModel.LazyItem):
            return False
        return self.rowCountHint(item) > 0

    def canFetchMore(self, index):
        if not index.isValid():
            return self.rowCount() < self.rowCountHint()
        item = self.itemFromIndex(index)
        if not isinstance(item, LazyItemModel.LazyItem):
            return False
        return item.rowCount() < self.rowCountHint(item)

    def fetchMore(self, index):
        if not index.isValid():
            rows = self.createRows(parent=None, start=self.rowCount())
            for r in rows:
                self.appendRow(r)
        else:
            item = self.itemFromIndex(index)
            if not isinstance(item, LazyItemModel.LazyItem):
                return
            rows = self.createRows(parent=item, start=item.rowCount())
            for r in rows:
                assert isinstance(r[0], LazyItemModel.LazyItem)
                item.appendRow(r)

    def rowCountHint(self, parent=None):
        # type: (Optional[LazyItemModel.LazyItem]) -> int
        """
        Return the a true (or estimated) row count for a parent item

        Subclasses need to override this method
        """
        raise NotImplementedError

    def createRows(self, parent=None, start=0):
        # type: (Optional[LazyItemModel.LazyItem], int) -> List[QStandardItem]
        """
        Create and return new rows for parent from start row index forward

        Subclasses need to override this method
        """
        raise NotImplementedError


def graph_clone(g):
    # type: (Graph) -> Graph
    return {n: list(succ) for n, succ in g.items()}


class GraphModel(LazyItemModel):
    class NodeItem(LazyItemModel.LazyItem):
        def node(self):
            return self.data(Qt.UserRole)

    class Column(enum.IntEnum):
        Node = 0

    HeaderLabels = [
        (Column.Node, "Node")
    ]

    def __init__(self, parent=None, graph={}, **kwargs):
        super().__init__(parent, **kwargs)
        self.setHorizontalHeaderLabels([t for _, t in self.HeaderLabels])
        self.__graph = {}      # type: Graph
        self.__graph_inv = {}  # type: Graph
        self.__roots = []      # type: List[Any]

        if graph:
            self.setGraph(graph)

    def graph(self):
        # type: () -> Graph
        return graph_clone(self.__graph)

    def setGraph(self, graph):
        # type: (Dict[T, List[T]]) -> None
        self.setRowCount(0)
        self.__graph = graph_clone(graph)
        self.__graph_inv = graph_invert(self.__graph)
        self.__roots = [n for n, pred in self.__graph_inv.items()
                        if not pred]

    def rowCountHint(self, parent=None):
        if parent is None:
            return len(self.__roots)
        else:
            assert isinstance(parent, GraphModel.NodeItem)
            node = parent.node()
            return len(self.__graph.get(node, []))

    def createRows(self, parent=None, start=0):
        if parent is None:
            return [self.createRow(n) for n in self.__roots[start:]]
        else:
            assert isinstance(parent, GraphModel.NodeItem)
            node = parent.node()
            return [self.createRow(n)
                    for n in self.__graph.get(node, [])[start:]]

    def createRow(self, node):
        item = GraphModel.NodeItem()
        item.setData(node, Qt.DisplayRole)
        item.setData(node, Qt.UserRole)
        return [item]


class OntologyModel(LazyItemModel):
    """
    Lazily populated model adapter for a Tree view of an ontology
    """
    class OntologyTerm(LazyItemModel.LazyItem):
        def term(self):
            value = self.data(Qt.UserRole)
            if isinstance(value, ontology.OBOObject):
                return value
            else:
                return None

    class Column(enum.IntEnum):
        Name, TermID, Definition = range(3)

    HeaderLabels = [
        (Column.Name, "Name"),
        (Column.TermID, "Term"),
        (Column.Definition, "Definition"),
    ]

    def __init__(self, ontology, parent=None, **kwargs):
        super().__init__(parent, **kwargs)
        self.ontology = ontology  # type: ontology.OBOOntology
        self.setHorizontalHeaderLabels(
            [text for _, text in self.HeaderLabels])

        roots = [t for t in self.ontology.root_terms()
                 if not t.is_obsolete]
        self.__data = namespace(roots=roots)

    def makeRow(self, term):
        # type: (ontology.Term) -> List[QStandardItem]
        htmlsummary = Term_html(term)
        first = OntologyModel.OntologyTerm(term.name)
        first.setData(term, Qt.UserRole)
        first.setToolTip(htmlsummary)

        second = QStandardItem(term.id)
        second.setData(htmlsummary, Qt.ToolTipRole)
        second.setToolTip(htmlsummary)
        third = QStandardItem(getattr(term, "def_", ""))
        third.setToolTip(third.text())
        return [first, second, third]

    def rowCountHint(self, parent=None):
        if parent is None:
            return len(self.__data.roots)
        else:
            assert isinstance(parent, OntologyModel.OntologyTerm)
            return len(self.ontology.child_edges(parent.term()))

    def createRows(self, parent=None, start=0):
        if parent is None:
            # top level
            return [self.makeRow(term) for term in self.__data.roots]
        term = parent.term()
        assert isinstance(term, ontology.Term)
        edges = self.ontology.child_edges(term)
        return [self.makeRow(child) for rtype, child in edges[start:]]


class OntologyEnrichmentResultsModel(GraphModel):
    """
    Tree model like representation of ontology enrichment results
    """
    Column = EnrichmentTableModel.Column
    HeaderLabels = EnrichmentTableModel.HeaderLabels

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.results = {}  # type: Dict[ontology.Term, Result]
        self.setColumnCount(len(self.HeaderLabels))

    def setResults(self, results):
        # type: (List[Result]) -> None
        self.setRowCount(0)
        self.results = {r.node: r for r in results}

    def setGraph(self, graph):
        # type: (Dict[ontology.Term, List[ontology.Term]]) -> None
        super().setGraph(graph)

    def createRow(self, node):
        assert isinstance(node, ontology.Term)
        res = self.results.get(node, None)
        htmlsummary = Term_html(node)
        termitem = GraphModel.NodeItem(node.name)
        termitem.setData(node, Qt.UserRole)
        termitem.setToolTip(htmlsummary)
        termitem.setWhatsThis(htmlsummary)
        termid = QStandardItem(node.id)
        termid.setToolTip(htmlsummary)
        termid.setWhatsThis(htmlsummary)

        if res is None:
            return [termitem, termid]

        count = QStandardItem()
        count.setData(len(res.query_mapped), Qt.DisplayRole)
        ref = QStandardItem()
        ref.setData(len(res.reference_mapped), Qt.DisplayRole)

        pval = QStandardItem()
        pval.setData(res.p_value, Qt.DisplayRole)
        pval.setToolTip(str(res.p_value))
        fdrval = QStandardItem()
        fdrval.setData(res.fdr_value, Qt.DisplayRole)
        fdrval.setToolTip(str(res.fdr_value))

        querylist = list(res.query_mapped)

        def tooltip():
            return "<span>" + ", ".join(map(escape, querylist)) + "</span>"

        genelist = _Item(
            data={
                Qt.DisplayRole: lambda:
                    join_elided(", ", 80, querylist, " (and {} more)"),
                Qt.ToolTipRole: tooltip,
                Qt.UserRole: lambda: querylist
            }
        )
        enrich = QStandardItem()
        enrich.setData(res.enrichment_score, Qt.DisplayRole)

        logenrich = QStandardItem()
        logenrich.setData(float(numpy.log2(res.enrichment_score)),
                          Qt.DisplayRole)
        row = [termitem, termid, count, ref, fdrval, pval, genelist,
               enrich, logenrich]
        return row


class FilterProxyModel(QSortFilterProxyModel):
    """
    A simple filter proxy model with settable filter predicates

    Example
    -------
    >>> proxy = FilterProxyModel()
    >>> proxy.setFilters([
    ...     FilterProxyModel.Filter(0, Qt.DisplayRole, lambda value: value < 1)
    ... ])
    """
    Filter = namedtuple("Filter", ["column", "role", "predicate"])

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.__filters = []

    def setFilters(self, filters):
        # type: (Sequence[FilterProxyModel.Filter]) -> None
        filters = [FilterProxyModel.Filter(f.column, f.role, f.predicate)
                   for f in filters]
        self.__filters = filters
        self.invalidateFilter()

    def filterAcceptsRow(self, row, parent):
        source = self.sourceModel()

        def apply(f):
            index = source.index(row, f.column, parent)
            data = source.data(index, f.role)
            try:
                return f.predicate(data)
            except (TypeError, ValueError):
                return False

        return all(apply(f) for f in self.__filters)


def relation_list_to_multimap(rellist):
    """
    Convert a 'relations' list into a multimap

    Parameters
    ----------
    rellist : Iterable[Tuple[Hashable[K], Any]]

    Returns
    -------
    multimap : Dict[K, List[Any]]

    Example
    -------
    >>> relation_list_to_multimap([(1, "a"), (2, "c"), (1, "b")])
    {1: ['a', 'b'], 2: ['c']}
    """
    mmap = defaultdict(list)
    for key, val in rellist:
        mmap[key].append(val)
    return dict(mmap)


def multimap_to_relation_list(multimap):
    """
    Parameters
    ----------
    multimap : Dict[Any, List[Any]]

    Returns
    -------
    rellist : List[Tuple[Any, Any]]

    Example
    -------
    >>> multimap_to_relation_list({1: ['a', 'b'], 2: ['c']})
    [(1, 'a'), (1, 'b'), (2, 'c')]
    """
    rellist = []
    for key, values in multimap.items():
        rellist.extend([(key, val) for val in values])
    return rellist


def multimap_invert(multimap):
    """
    Invert a multimap

    Parameters
    ----------
    multimap : Dict[Hashable, List[Hashable]]

    Returns
    -------
    inverted : Dict[Hashable, List[Hashable]]

    Example
    -------
    >>> m = multimap_invert({1: ["a", "b"], 2: ["c"]})
    >>> assert m == {'a': [1], 'b': [1], 'c': [2]}
    """
    rellist = multimap_to_relation_list(multimap)
    return relation_list_to_multimap((v, k) for k, v in rellist)


def simple_graph(ontology, reltypes=["is_a"]):
    """
    Convert an OBOOntology instance into a simple Dict[str, List[str]]
    encoded successor graph.

    Parameters
    ----------
    ontology : ontology.OBOOntology
    reltypes : List[str]
        Specify the subset of relationship types to consider.

    Returns
    -------
    graph : Dict[str, List[str]]
        Ontology encoded as a graph. Keys are ontology term ids (e.g.
        'GO:0000001')
    """
    reltypes = set(reltypes)
    termids = [t.id for t in ontology.terms()]
    return {termid: [t.id for rtype, t in ontology.child_edges(termid)
                     if rtype in reltypes]
            for termid in termids}


def _dag_enrichment(graph, query, reference, querycount, refcount):
    """
    Internal. `dag_enrichment` dispatches to this function

    Parameters
    ----------
    graph : Dict[T, List[T]]
    query : Dict[T, List[M]]
    reference : Dict[T, List[M]]
    querycount : int
    refcount: int

    Returns
    -------
    results : List[Result[T, M]]
    """
    # map query/reference items onto the dag
    refset = dag_annotate_closure(graph, reference)
    queryset = dag_annotate_closure(graph, query)
    querynodes = [n for n, _ in queryset]
    queryset = dict(queryset)
    refset = dict(refset)

    probmodel = stats.Hypergeometric()

    results = []

    nan = float("nan")
    for node in querynodes:
        mapped_query = queryset.get(node, [])
        if not mapped_query:
            continue
        mapped_ref = refset.get(node, [])
        p_query = len(mapped_query) / querycount if querycount else nan
        p_ref = len(mapped_ref) / refcount if refcount else nan
        enrichment = p_query / p_ref if p_ref else nan
        pval = probmodel.p_value(len(mapped_query), refcount,
                                 len(mapped_ref), querycount)
        results.append(
            Result(node, mapped_query, mapped_ref, pval, nan,
                   enrichment)
        )
    fdr_vals = stats.FDR([r.p_value for r in results])
    results = [r._replace(fdr_value=min(f, 1.0), p_value=min(r.p_value, 1.0))
               for r, f in zip(results, fdr_vals)]
    return results


def dag_enrichment(graph, query, reference):
    # type: (Dict[T, List[T]], Dict[M, T], Dict[M, T]) -> Result[T, M]
    """
    Parameters
    ----------
    graph : Dict[T, List[T]]
        A DAG encoded as a successors mapping.
    query : Dict[M, List[T]]
        Query annotations.
    reference : Dict[M, List[T]]
        Reference annotations.

    Returns
    -------
    results : List[Result[T, M]]
    """
    querycount = len(query)
    refcount = len(reference)
    querymap = multimap_invert(query)
    refmap = multimap_invert(reference)
    return _dag_enrichment(graph, querymap, refmap, querycount, refcount)


_Result = namedtuple(
    "Result", [
        "node",              # type: T
        "query_mapped",      # type: list
        "reference_mapped",  # type: list
        "p_value",           # type: float
        "fdr_value",         # type: float
        "enrichment_score",  # type: float
    ]
)

T = TypeVar("T")
M = TypeVar("M")

Graph = Dict[T, List[T]]


class Result(_Result, Generic[T, M]):
    """
    Attributes
    ----------
    node : T
    query_mapped : List[M]
    reference_mapped : List[M]
    p_value : float
    fdr_value : float
    enrichement_score : float
    """
    node = _Result.node  # type: T
    query_mapped = _Result.query_mapped  # type: List[M]
    reference_mapped = _Result.reference_mapped  # type: List[M]
    p_value = _Result.p_value  # type: float
    fdr_value = _Result.fdr_value  # type: float
    enrichment_score = _Result.enrichment_score  # type: float


def unique(iterable):
    """
    Return an iterator over unique items of `iterable`.

    The order of the (first encountered) unique items is preserved.

    Parameters
    ----------
    iterable : Iterable[Any]

    Returns
    -------
    iter : Iterable
    """
    seen = set()
    for el in iterable:
        if el not in seen:
            seen.add(el)
            yield el


def graph_invert(graph):
    # type: (Graph) -> Graph
    """
    Return a graph with all edges reversed

    Parameters
    ----------
    graph : Dict[T, List[T]]

    Returns
    -------
    graph_inv: Dict[T, List[T]]
    """
    inverted = defaultdict(list)
    for node, successors in graph.items():
        for succ in successors:
            inverted[succ].append(node)

    # complete the relations mapping for leaf/root nodes
    for node in set(graph) - set(inverted):
        inverted[node] = []

    return dict(inverted)


def closure(graph):
    # type: (Graph) -> Graph
    """
    Compute the full transitive closure of a DAG.

    Parameters
    ----------
    graph : Dict[T, List[T]]
        DAG encoded as a successor mapping.

    Returns
    -------
    closure : Dict[T, List[T]]
        Full transitive closure of `graph`.

    Example
    -------
    >>> closure({1: [2, 3], 4: [2], 2: [3], 5: []})
    {1: [2, 3], 2: [3], 3: [], 4: [2, 3], 5: []}
    """
    rorder = reversed(topologicalsort(graph))
    result = {}
    for node in rorder:
        result[node] = list(unique(chain.from_iterable(
            ([succ] + result[succ]) for succ in graph.get(node, [])))
        )
    return result


def dag_annotate_closure(graph, annotations):
    # type: (Graph, Dict[T, Hashable]) -> List[Tuple[T, List[Hashable]]]
    """
    Propagate a set of annotations in a DAG onto its closure.

    Parameters
    ----------
    graph : Dict[T, List[T]]
        A DAG encoded as a node -> successors dictionary.
    annotations : Dict[T, List[Hashable]]
        A dictionary of annotations that map directly to a node.

    Returns
    -------
    result : List[Tuple[T, List[Hashable]]]
        For every node (`T`) in a graph a list of all annotations which are
        in the node's transitive closure.
    """
    rorder = reversed(topologicalsort(graph))
    result = []
    cache = {}
    for node in rorder:
        node_annot = [annotations.get(node, [])]
        for succ in graph.get(node, []):
            node_annot.append(cache[succ])
        node_annot = list(unique(chain.from_iterable(node_annot)))
        cache[node] = node_annot
        result.append((node, node_annot))
    return result


def topologicalsort(graph):
    """
    Return nodes in a DAG sorted topologically.

    Parameters
    ----------
    graph : Dict[T, List[T]]
        DAG expressed as a node -> successors mapping, that is for every key
        in graph graph[key] returns immediate successors to key

    Returns
    -------
    nodes : List[T]
    """
    successors = {node: list(successors)
                  for node, successors in graph.items()}
    incomming = graph_invert(successors)

    # complete the relations mapping for leaf/root nodes if missing
    for node in set(incomming) - set(successors):
        successors[node] = []

    topqueue = [node for node in incomming if len(incomming[node]) == 0]
    topqueue = deque(topqueue)

    result = []
    while topqueue:
        current = topqueue.popleft()
        for succ in successors[current]:
            incomming[succ].remove(current)
            if len(incomming[succ]) == 0:
                topqueue.append(succ)
        assert len(incomming[current]) == 0
        del incomming[current]
        result.append(current)
    if incomming:
        raise ValueError("Cycles on graph")
    return result


# TODO: wrong name, its not a subgraph (path contraction?)
def dag_subgraph(graph, nodeset):
    """
    Iteratively remove nodes from graph while preserving connectivity.

    That is if b is removed and a a -> b -> c is a path then an explicit
    a -> c edge is inserted. (i.e. perform path contraction on all paths
    through b)

    Parameters
    ----------
    graph : Graph
        Input graph
    nodeset : set
        A set of nodes to *preserve*

    Returns
    -------
    graph : Graph
    """
    order = topologicalsort(graph)
    successors = {n: list(e) for n, e in graph.items()}
    predeccessor = graph_invert(graph)
    subgraph = {}

    for n in reversed(order):
        if n in nodeset:
            subgraph[n] = successors[n]
        else:
            for pred in predeccessor[n]:
                # replace the pred -> n edge with pred -> successors[n]
                successors[pred].remove(n)
                for k in successors[n]:
                    if k in nodeset and k not in successors[pred]:
                        successors[pred].append(k)
    return subgraph


def graph_select_closure(graph, nodeset):
    """
    Return a subgraph from graph such that it contains all nodes in `nodeset`
    as well as nodes in the transitive closure of `nodeset`

    Parameters
    ----------
    graph : Graph
    nodeset : Iterable

    Returns
    -------
    graph : Graph
    """
    # TODO: same as subgraph(graph, dag_closure(graph, nodest))
    # where subgraph extracts only nodes and edged between them
    newgraph = {}
    queue = deque(unique(nodeset))
    visited = set()
    while queue:
        n = queue.popleft()
        if n in visited:
            continue
        newgraph[n] = list(graph.get(n, []))
        queue.extend(newgraph[n])
        visited.add(n)
    return newgraph


def mapman_mapping(stream):
    """
    Read a MapMan formatted excel-tab mapping file.

    Parameters
    ----------
    stream : Iterable[str]
        An iterable yielding text lines of the mapping file.

    Returns
    -------
    r : Tuple[List[mapman_mapping.Node], List[mapman_mapping.Record]]
    """
    reader = csv.reader(stream, dialect="excel-tab")
    return mapman_mapping_from_rows(reader)


def mapman_mapping_from_rows(rowsiter):
    try:
        header = next(rowsiter)
        if header[:4] != ["BINCODE", "NAME", "IDENTIFIER", "DESCRIPTION"]:
            warnings.warn("First row is not an expected header for a "
                          "MapMan mapping ({})"
                          .format(header), UserWarning)
    except StopIteration:
        return []

    mapping = []  # type: List[mapman_mapping.Record]
    # child -> parent relationship (is a tree)
    nodes = []  # type: List[mapman_mapping.Node]
    bincode_to_node = {}  # type: Dict[str, mapman_mapping.Node]
    graph = OrderedDict()
    for row in rowsiter:
        bincode, name = row[:2]
        parent, _, _ = bincode.rpartition('.')
        if bincode in nodes and nodes[bincode].name != name:
            raise ValueError("Inconsistent names ('{}' != '{}')"
                             .format(nodes[bincode].name, name))
        if bincode not in bincode_to_node:
            # first time seen
            node = mapman_mapping.Node(bincode, name, [])
            bincode_to_node[bincode] = node
            nodes.append(node)

            if parent:
                parent_node = bincode_to_node[parent]
                parent_node.children.append(node)

        if len(row) >= 3:
            mapping.append(
                mapman_mapping.Record(*((row + ["", ""])[:5]))
            )
    return nodes, mapping

mapman_mapping.Node = namedtuple(
    "Node", [
        "bincode",
        'name',
        'children',
    ]
)

mapman_mapping.Record = namedtuple(
    "Record", [
        "bincode",
        "name",
        "identifier",
        "description",
        "type",
    ]
)


def mapman_to_obo(nodes, idprefix):
    # type: (List[mapman_mapping.Node], str) -> ontology.OBOOntology
    """
    Convert/import an mapman nodes list to an OBOOntology instance

    Parameters
    ----------
    nodes : List[mapman_mapping.Node]
    idprefix : str
    """
    ont = ontology.OBOOntology()
    ont.add_header_tag("format-version", "1.2")
    ont.add_header_tag("auto-generated-by", __name__ + ".mapman_to_obo")
    terms = []
    for node in nodes:
        parent, _, _ = node.bincode.rpartition(".")
        term = ontology.Term(id=idprefix + ":" + node.bincode, name=node.name)
        terms.append(term)
        if parent:
            term.add_tag("is_a", idprefix + ":" + parent)
        ont.add_object(term)
    return ont


import unittest


class TestUtils(unittest.TestCase):
    def test_toposort(self):
        E1 = {}
        r = topologicalsort(E1)
        self.assertEqual(r, [])

        E2 = {1: []}
        r = topologicalsort(E2)
        self.assertEqual(r, [1])

        T1 = {1: [2, 3], 2: [4, 5], 3: [6, 7], 4: [], 5: [], 6: [], 7: []}
        r = topologicalsort(T1)
        self.assertEqual(r, [1, 2, 3, 4, 5, 6, 7])

        C1 = {1: [2], 2: [3], 3: [1]}
        with self.assertRaises(ValueError):
            r = topologicalsort(C1)

        C1 = {1: [1]}
        with self.assertRaises(ValueError):
            topologicalsort(C1)

        DAG1 = {1: [2, 3], 2: [3, 5], 3: [5], 5: []}
        r = topologicalsort(DAG1)
        self.assertEqual(r, [1, 2, 3, 5])

        DAG2 = {1: [2, 3], 2: [3], 3: [4], 4: []}
        r = topologicalsort(DAG2)
        self.assertEqual(r, [1, 2, 3, 4])

        G2 = {1: [3], 2: [3, 4], 3: [4], 4: [5], 5: [6, 7, 8],
              6: [7], 7: [8], 8: [9], 9: []}
        r = topologicalsort(G2)
        self.assertEqual(r, [1, 2, 3, 4, 5, 6, 7, 8, 9])

    def test_graph_invert(self):
        self.assertEqual(graph_invert({}), {})
        self.assertEqual(graph_invert({1: [], 2: []}), {1: [], 2: []})
        self.assertEqual(
            graph_invert({1: [2, 3], 4: [2]}),
            {1: [], 2: [1, 4], 3: [1], 4: []})
        self.assertEqual(
            graph_invert({1: [2, 3], 4: [2], 2: [3], 5: []}),
            {1: [], 2: [1, 4], 3: [1, 2], 4: [], 5: []})

    def test_closure(self):
        self.assertEqual(closure({}), {})
        self.assertEqual(closure({1: []}), {1: []})
        self.assertEqual(closure({1: [], 2: []}), {1: [], 2: []})

        T = {1: [2, 4], 2: [3], 4:[5]}
        self.assertEqual(closure(T), {1: [2, 3, 4, 5], 2: [3], 4: [5],
                                      3: [], 5: []})
        G = {1: [2, 3], 4: [2], 2: [3], 5: []}
        CG = {1: [2, 3], 2: [3], 3: [], 4: [2, 3], 5: []}
        self.assertEqual(closure(G), CG)

        G = {1: [2, 3], 2: [4], 3: [5], 4: [6], 5: [6]}
        GC = {1: [2, 4, 6, 3, 5], 2: [4, 6], 3: [5, 6], 4:[6], 5: [6], 6: []}
        self.assertEqual(closure(G), GC)

    def test_dag_annotate_closure(self):
        T1 = {1: [2], 2: [3, 4]}
        r = dag_annotate_closure(T1, {3: ["3"], 2: ["2"]})
        expected = {3: ["3"], 4: [], 2: ["3", "2"], 1: ["3", "2"]}
        for nid, annot in r:
            self.assertSetEqual(set(annot), set(expected[nid]))

        T2 = {1: [2]}
        r = dag_annotate_closure(T2, {2: ["2a", "2b"]})
        expected = {1: ["2a", "2b"], 2: ["2a", "2b"]}
        for nid, annot in r:
            self.assertSetEqual(set(annot), set(expected[nid]))

        G2 = {1: [2, 3], 4: [2], 2: [3], 5: []}
        r = dag_annotate_closure(G2, {3: ["3"], 2: ["2"], 5: ["5"]})
        expected = {3: ['3'], 5: ['5'], 2: ['2', '3'], 1: ['2', '3'],
                    4: ['2', '3']}
        for nid, annot in r:
            self.assertSetEqual(set(annot), set(expected[nid]))

        G3 = {1: [2, 6], 2: [], 3: [4, 5], 4: [2], 5: [6], 6: []}
        annot = {1: ["1"], 4: ["4|5"], 5: ["4|5", "5"], 6: []}
        expected = {1: ["1"], 2: [], 3: ["4|5", "5"], 4: ["4|5"],
                    5: ["4|5", "5"], 6: []}
        r = dag_annotate_closure(G3, annot)
        for nid, annot in r:
            self.assertSetEqual(set(annot), set(expected[nid]))

    def test_dag_enrichhment(self):
        E1 = {}
        res = _dag_enrichment(E1, {1: ["1a", "1b"]}, {}, 2, 0)
        self.assertEqual(res, [])

        S1 = {1: []}
        query = {1: ["a"]}
        ref = {1: ["a"]}
        res = _dag_enrichment(S1, query, ref, 1, 1)
        self.assertEqual(len(res), 1)
        res = res[0]
        assert isinstance(res, Result)
        self.assertEqual(res.node, 1)
        self.assertEqual(res.query_mapped, ["a"])
        self.assertEqual(res.reference_mapped, ["a"])
        self.assertEqual(res.p_value, 1.0)
        self.assertEqual(res.enrichment_score, 1.0)

        P1 = {1: [2]}
        query = {2: ["2a"]}
        ref = {1: ["1a"], 2: ["2a"]}
        res = _dag_enrichment(P1, query, ref, 1, 2)
        self.assertEqual(len(res), 2)
        [r2, r1] = res
        # TODO: Should this order be required?
        self.assertEqual((r1.node, r2.node), (1, 2), "not in rev. topological order")
        self.assertEqual(r1.query_mapped, ["2a"])
        self.assertSetEqual(set(r1.reference_mapped), {"1a", "2a"})
        self.assertEqual(r1.p_value, 1.0)
        self.assertEqual(r1.enrichment_score, 1.0)

        self.assertEqual(r2.query_mapped, ["2a"])
        self.assertEqual(r2.reference_mapped, ["2a"])
        self.assertEqual(r2.p_value, 0.5)
        self.assertEqual(r2.enrichment_score, 2.0)

        res = _dag_enrichment(
            P1, {1: ["1a", "1b"]}, {2: ["2a", "2b", "2c", "2d"]}, 2, 4)
        G1 = {1: [2, 4], 2: [3], 4: [3, 5]}
        res = _dag_enrichment(G1, {1: ["1a"], 4: ["4a", "4b"]},
                             {5: ["5a", "5b"], 3: ["3a", "3b"]}, 3, 4)

        G1 = {1: [2, 3], 4: [2], 2: [3], 5: []}


class TestMapMan(unittest.TestCase):
    def test_mapping_parse(self):
        mapman_mapping_f = io.StringIO(
            "BINCODE\tNAME\tIDENTIFIER\tDESCRIPTION\n"
            "1\tA\n"
            "1\tA\tG1\tsomething\n"
            "1\tA\tG2\tsomething else\n"
            "1.1\tA.B\n"
            "1.1\tA.B\tG1\tsomething\n"
        )

        n11 = mapman_mapping.Node("1.1", "A.B", [])
        n1 = mapman_mapping.Node("1", "A", [n11])
        expected_nodes = [n1, n11]
        expected_rcs = [
            mapman_mapping.Record("1", "A", "G1", "something", ""),
            mapman_mapping.Record("1", "A", "G2", "something else", ""),
            mapman_mapping.Record("1.1", "A.B", "G1", "something", ""),
        ]
        nodes, rcs = mapman_mapping(mapman_mapping_f)
        self.assertSequenceEqual(nodes, expected_nodes)
        self.assertSequenceEqual(rcs, expected_rcs)

    def test_mapman_to_obo(self):
        mmtree = [
            mapman_mapping.Node(
                "1", "a", [mapman_mapping.Node("1.1", "a.a", [])]),
            mapman_mapping.Node(
                "2", "b", [mapman_mapping.Node("2.1", "b.b", [])]),
        ]
        from Orange.clustering.hierarchical import preorder
        mmnodes = chain.from_iterable(
            preorder(mn, branches=lambda n: n.children) for mn in mmtree)
        mmnodes = list(mmnodes)
        ont = mapman_to_obo(mmnodes, "MM")
        for mn in mmtree:
            oterm = ont.term("MM:" + mn.bincode)
            self.assertEqual(oterm.name, mn.name)

            for ch in mn.children:
                chterm = ont.term("MM:" + ch.bincode)
                self.assertIn(("is_a", oterm), ont.parent_edges(chterm))


def example_ontology_view():
    from AnyQt.QtWidgets import QWidget, QVBoxLayout, QTextBrowser
    app = QApplication(list(sys.argv))
    argv = app.arguments()
    if len(argv) > 1:
        obofile = argv[1]
    else:
        res = GeneOntologyResource()
        res.fetch(progress=print)
        obofile = res.localpath()

    with open(obofile, "rb") as f:
        if obofile.endswith(".gz"):
            f = gzip.GzipFile(fileobj=f)
        ftext = io.TextIOWrapper(f, encoding="utf-8")
        ont = ontology.OBOOntology(ftext)
        ftext.detach()

    view = QTreeView()
    view.show()
    tb = QTextBrowser()
    w = QWidget()
    w.setLayout(QVBoxLayout())
    w.layout().setContentsMargins(0, 0, 0, 0)
    splitter = QSplitter(orientation=Qt.Vertical)
    splitter.addWidget(view)
    splitter.addWidget(tb)
    splitter.setSizes([100, 75])

    w.layout().addWidget(splitter)
    model = OntologyModel(ont)
    view.setModel(model)

    def update_displayed_term(selected, deselected):
        rows = view.selectionModel().selectedRows(0)
        if rows:
            term = model.data(rows[0], Qt.UserRole)
        else:
            term = None
        if not isinstance(term, ontology.Term):
            tb.setHtml("")
        else:
            tb.setHtml(Term_html(term))

    view.selectionModel().selectionChanged.connect(update_displayed_term)
    w.resize(700, 600)
    view.header().resizeSection(0, 300)
    w.show()
    w.raise_()
    app.exec_()


def main(argv=sys.argv):
    app = QApplication(list(argv))
    argv = app.arguments()
    w = OWMapManEnrichment()

    w.show()
    w.raise_()

    if len(argv) > 1:
        dataset = Orange.data.Table(argv[1])
    else:
        genelist = numpy.array(
            [[
                "ACS6",
                "AT1G72520",
                "AT3G02840",
                "BCB",
                "CML37",
                "ERF6",
                "HSFA2",
                "HSP17.6A",
                "RHL41",
                "STZ",
                "WRKY53",
            ]],
            dtype=object).T
        dataset = Orange.data.Table.from_numpy(
            Orange.data.Domain(
                [], [], [Orange.data.StringVariable("gene")]
            ),
            numpy.zeros((genelist.shape[0], 0)), metas=genelist
        )

    w.set_query_data(dataset)
    w.handleNewSignals()
    app.exec_()
    w.set_query_data(None)
    w.handleNewSignals()
    w.saveSettings()
    w.onDeleteWidget()
    return 0

if __name__ == "__main__":
    # example_ontology_view()
    main()
