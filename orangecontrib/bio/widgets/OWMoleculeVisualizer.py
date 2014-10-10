
from __future__ import absolute_import

import math

from cStringIO import StringIO
from collections import namedtuple

from PyQt4.QtGui import (
    QWidget, QFrame, QLabel, QListWidget, QScrollArea, QSplitter,
    QListView, QVBoxLayout, QGridLayout, QSizePolicy,
    QPainter, QBrush, QPen, QPalette
)

from PyQt4.QtCore import Qt, QSize, QRectF, QByteArray, QTimer
from PyQt4.QtCore import pyqtSlot as Slot, pyqtSignal as Signal

from PyQt4.QtSvg import QSvgWidget

import Orange.data
from Orange.utils import environ

from Orange.OrangeWidgets import OWGUI
from Orange.OrangeWidgets.OWWidget import \
    OWWidget, DomainContextHandler, ContextField, Default

from Orange.OrangeWidgets.OWItemModels import VariableListModel

NAME = "Molecule Visualizer"
DESCRIPTION = "Rendering of 2D structure of molecules based on " \
              "their SMILES description."
ICON = "icons/MoleculeVisualizer.svg"
PRIORITY = 2050

INPUTS = [("Molecules", Orange.data.Table, "setMoleculeTable", Default),
          ("Molecule subset", Orange.data.Table, "setMoleculeSubset"),
          ("Fragments", Orange.data.Table, "setFragmentTable")]
OUTPUTS = [("Selected Molecules", Orange.data.Table)]

REPLACES = ["_bioinformatics.widgets.OWMoleculeVisualizer.OWMoleculeVisualizer"]

try:
    import oasa
except ImportError:
    oasa = None

try:
    import pybel
except ImportError:
    pybel = None
except Exception, ex:
    pybel = None

svg_error_string = """<?xml version="1.0" ?>
<svg height="185" version="1.0" width="250" xmlns="http://www.w3.org/2000/svg">
    <g stroke="#000" stroke-width="1.0">
        <text font-family="Arial" font-size="16" stroke="rgb(0, 0, 0)" x="0.0" y="16.0">
            %s
        </text>
    </g>
</svg>
"""


def molecule_from_smiles(smiles, ):
#     converter = oasa.smiles.converter()
#     converter.configuration["R_BOND_LENGTH"] = 30
#     mols = converter.read_text(smiles)
#     if len(mols) > 1:
#         pass
#     return mols[0]

    return oasa.smiles.text_to_mol(smiles, calc_coords=30)


class svg_out(oasa.svg_out.svg_out):
    def __init__(self, ):
        self._highlight = set()

    def mol_with_substructure_to_svg(self, mol, embedding):
        self._highlight = set(embedding)
        res = self.mol_to_svg(mol)
        self._highlight = set()
        return res

    def _draw_edge(self, e):
        v1, v2 = e.vertices
        matched = v1 in self._highlight and v2 in self._highlight
        if matched:
            self._curr_color = "rgb(255, 0, 0)"
        else:
            self._curr_color = "rgb(0, 0, 0)"

        oasa.svg_out.svg_out._draw_edge(self, e)

    def _draw_vertex(self, v):
        if v in self._highlight:
            self._curr_color = "rgb(255, 0, 0)"
        else:
            self._curr_color = "rgb(0, 0, 0)"
        oasa.svg_out.svg_out._draw_vertex(self, v)

    def _draw_line(self, parent, start, end, line_width=1, capstyle=""):
        x1, y1 = start
        x2, y2 = end

        attrs = [
            ("x1", str(x1)), ("y1", str(y1)),
            ("x2", str(x2)), ("y2", str(y2)),
            ("stroke", self._curr_color)
        ]
        line = self._element("line", attrs)

        parent.appendChild(line)

    def _draw_text(self, parent, xy, text, font_name="Arial", font_size=16):
        x, y = xy
        attrs = [
            ("x", str(x)),
            ("y", str(y)),
            ("font-family", font_name),
            ("font-size", str(font_size)),
            ("stroke", self._curr_color)
        ]
        text_el = self._element("text", attrs)
        textnode = self.document.createTextNode(text)
        text_el.appendChild(textnode)
        parent.appendChild(text_el)

    def _element(self, tag, attrs=[]):
        el = self.document.createElement(tag)
        for key, val in attrs:
            el.setAttribute(key, val)
        return el


def molecule_to_svg(mol):
    """
    :param oasa.molecule.molecule mol:
        Input molecule.
    """
    if any(atom.x is None or atom.y is None for atom in mol.atoms):
        generator = oasa.coords_generator.coords_generator()
        generator.calculate_coords(mol, bond_length=30)

    writer = oasa.svg_out.svg_out()
    xmldoc = writer.mol_to_svg(mol)
    return xmldoc.toprettyxml()


def molecule_to_svg_with_substructure(mol, embedding):
    """
    :param oasa.molecule.molecule mol:
        Input molecule.
    :param list embedding:
        A list of atom to highlight.
    """
    if any(atom.x is None or atom.y is None for atom in mol.atoms):
        generator = oasa.choord_generator.coord_generator()
        generator.calculate_coords(mol, bond_length=30, force=True)

    writer = svg_out()
    xmldoc = writer.mol_with_substructure_to_svg(mol, embedding)
    return xmldoc.toprettyxml()


def substructure_embedding(mol, pattern):
    obmol, atom_map = (oasa.pybel_bridge.PybelConverter
                           .oasa_to_pybel_molecule_with_atom_map(mol))

    if not isinstance(pattern, pybel.Smarts):
        pattern = pybel.Smarts(pattern)

    matches = pattern.findall(obmol)

    rev_map = dict(map(reversed, atom_map.iteritems()))

    return [map(rev_map.get, match) for match in matches]


def is_valid_smiles(string):
    try:
        oasa.smiles.text_to_mol(string, calc_coords=False,
                                localize_aromatic_bonds=False)
    except Exception:
        return False
    else:
        return True


class ThumbnailWidget(QFrame):

    def __init__(self, parent, svg_data="", text=""):
        QFrame.__init__(self, parent)

        self.setFrameStyle(QFrame.StyledPanel)

        self.label = QLabel(text, self)
        self.label.setWordWrap(True)

        self.label.setSizePolicy(QSizePolicy.Ignored,
                                 QSizePolicy.Expanding)
        self.label.setAlignment(Qt.AlignCenter)
        self.label.setAttribute(Qt.WA_TransparentForMouseEvents)

        self.svg = SVGThumbnailWidget(svg_data, parent=self)
        self.svg.setAttribute(Qt.WA_TransparentForMouseEvents)

        layout = QVBoxLayout()
        layout.addWidget(self.label, stretch=1, alignment=Qt.AlignCenter)
        layout.addWidget(self.svg, stretch=4, alignment=Qt.AlignCenter)

        self.setLayout(layout)
        self.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)

        self._selected = False

    def setSelected(self, val):
        if self._selected != val:
            self._selected = val
            self.svg.setSelected(val)

    def getSelected(self):
        return self._selected

    selected = property(getSelected, setSelected)

    def setData(self, data):
        self.svg.setData(data)

    def setImageSize(self, width, height):
        self.svg.setFixedSize(width, height)

    def paintEvent(self, event):
        QFrame.paintEvent(self, event)

    def mousePressEvent(self, event):
        if event.button() & Qt.LeftButton:
            parent = self.parent()
            if event.modifiers() & Qt.ControlModifier:
                self.setSelected(not self._selected)
            else:
                for child in parent.findChildren(ThumbnailWidget):
                    child.setSelected(False)

                self.setSelected(True)

            parent.selectionChanged.emit()

        super(ThumbnailWidget, self).mousePressEvent(event)

    def mouseDoubleClickEvent(self, event):
        self._bigimage = BigSvgWidget()
        self._bigimage.load(QByteArray(self.svg._data))
        self._bigimage.show()

    def customEvent(self, event):
        self.layout().invalidate()
        self.update()


class SVGThumbnailWidget(QSvgWidget):

    def setSelected(self, val):
        if val != self._selected:
            self._selected = val
            self.update()

    def getSelected(self):
        return self._selected

    selected = property(getSelected, setSelected)

    @property
    def highlighted(self):
        return self._highlighted

    def __init__(self, svg_data="", parent=None):
        super(SVGThumbnailWidget, self).__init__(parent)

        self._data = svg_data
        self._cached = None
        self._selected = False
        self._highlighted = False

        self.resize(128, 128)

        if svg_data:
            self.load(QByteArray(svg_data))

        self.setBackgroundRole(QPalette.Base)

    def setData(self, svgdata):
        self._data = svgdata
        self.load(QByteArray(svgdata))

    def paintEvent(self, event):
        crect = self.rect()
        painter = QPainter(self)
        painter.setRenderHint(QPainter.Antialiasing)

        painter.setBrush(QBrush(Qt.white))
        painter.setPen(QPen(Qt.lightGray, 1.2))
        painter.drawRoundedRect(QRectF(crect).adjusted(2, 2, -2, -2), 2, 2,
                                Qt.AbsoluteSize)

        if self._selected:
            painter.setPen(QPen(QBrush(Qt.red), 2))
        if self._highlighted:
            painter.setBrush(QBrush(Qt.gray,  Qt.FDiagPattern))
        else:
            painter.setBrush(Qt.NoBrush)

        painter.drawRoundedRect(QRectF(crect).adjusted(2, 2, -2, -2), 2, 2,
                                Qt.AbsoluteSize)

        defsize = self.renderer().defaultSize()
        margin = 5

        bound = QSize(defsize)
        bound.scale(crect.width() - margin, crect.height() - margin,
                    Qt.KeepAspectRatio)

        svgrect = QRectF(0, 0, bound.width(), bound.height())
        svgrect.moveCenter(crect.center())
        self.renderer().render(painter, svgrect)

    def sizeHint(self):
        sh = self.renderer().defaultSize()
        size = max(sh.width(), sh.height())
        return QSize(size, size)

    def mouseDoubleClickEvent(self, event):
        self._bigimage = BigSvgWidget()
        self._bigimage.load(QByteArray(self._data))
        self._bigimage.show()


class BigSvgWidget(QSvgWidget):

    def paintEvent(self, event):
        painter = QPainter(self)
        painter.setBrush(QBrush(Qt.white))
        painter.setPen(Qt.NoPen)
        painter.drawRect(0, 0, self.width() - 1, self.height() - 1)
        painter.end()
        QSvgWidget.paintEvent(self, event)


class ScrollArea(QScrollArea):

    def __init__(self, *args, **kwargs):
        QScrollArea.__init__(self, *args, **kwargs)
        self.viewport().setMouseTracking(True)
        self.setMouseTracking(True)
        self.setWidget(QWidget())
        self.widget().setLayout(QGridLayout())

    def resizeEvent(self, event):
        QScrollArea.resizeEvent(self, event)
        self.reflow(self.viewport().width())

    def reflow(self, width):
        widget = self.widget()
        grid = widget.layout()

        items = [grid.itemAt(i) for i in range(grid.count())]
        sizes = [item.sizeHint() for item in items]
        wmax = reduce(max, (sh.width() for sh in sizes), 0)
        hspacing = grid.horizontalSpacing()
        left, _, right, _ = grid.getContentsMargins()

        # width >= wmax * ncols + (ncols - 1) * hspacing + left + right
        ncols = (width - left - right + hspacing) / (wmax + 1)
        ncols = max(1, math.floor(ncols))

        if ncols != grid.columnCount():
            for i in range(len(items) - 1, -1, -1):
                grid.takeAt(i)

            for i, item in enumerate(items):
                grid.addItem(item, i // ncols, i % ncols)


class GridWidget(QWidget):
    selectionChanged = Signal()

    def __init__(self, *args):
        super(GridWidget, self).__init__(*args)
        self.setLayout(QGridLayout())
        self.layout().setSizeConstraint(QGridLayout.SetMinAndMaxSize)
        self.__autoreflow = False

    def resizeEvent(self, event):
        super(GridWidget, self).resizeEvent(event)
        if self.__autoreflow:
            self.reflow(self.width())

    def appendWidget(self, widget):
        count = self.layout().count()
        ncol = self.layout().columnCount()
        self.layout().addWidget(widget, count // ncol, count % ncol)

    def count(self):
        return self.layout().count()

    def clear(self):
        for i in reversed(range(self.count())):
            item = self.layout().takeAt(i)
            if item.widget() is not None and item.widget().parent() is self:
                widget = item.widget()
                widget.setParent(None)
                widget.deleteLater()

    def reflow(self, width):
        grid = self.layout()

        items = [grid.itemAt(i) for i in range(grid.count())]
        sizes = [item.sizeHint() for item in items]
        wmax = reduce(max, (sh.width() for sh in sizes), 0)
        hspacing = grid.horizontalSpacing()
        left, _, right, _ = grid.getContentsMargins()

        # width >= wmax * ncols + (ncols - 1) * hspacing + left + right
        ncols = (width - left - right + hspacing) / (wmax + 1)
        ncols = max(1, math.floor(ncols))

        if ncols != grid.columnCount():
            for i in range(len(items) - 1, -1, -1):
                grid.takeAt(i)

            for i, item in enumerate(items):
                grid.addItem(item, i // ncols, i % ncols)


ParseError, StripedSalts, OK = 1, 2, 3

Item = namedtuple(
    "Item",
    ["index",
     "smiles",
     "status",
     "molecule",   # molecule | None
     ]
)


class OWMoleculeVisualizer(OWWidget):
    settingsList = ["colorFragmets", "showFragments"]

    contextHandlers = {
        "": DomainContextHandler(
            "",
            [ContextField("selected_title_indices"),
             ContextField("moleculeTitleAttributeList",
                          (DomainContextHandler.List +
                           DomainContextHandler.SelectedRequired +
                           DomainContextHandler.IncludeMetaAttributes),
                          selected="selectedMoleculeTitleAttrs"),
             ContextField("smiles_var",
                          DomainContextHandler.Required +
                          DomainContextHandler.IncludeMetaAttributes)],
             maxAttributesToPickle=1000)
    }

    def __init__(self, parent=None, signalManager=None,
                 title="Molecule visualizer"):
        super(OWMoleculeVisualizer, self).__init__(parent, signalManager, title)

        self.colorFragments = 1
        self.showFragments = 0
        self.selectedFragment = ""
        self.moleculeSmiles = []
        self.fragmentSmiles = []
        self.defFragmentSmiles = []
        self.smiles_var = 0
        self.moleculeTitleAttr = 0
        self.moleculeTitleAttributeList = []
        self.selectedMoleculeTitleAttrs = []
        self.fragmentSmilesAttr = 0
        self.imageSize = 200
        self.numColumns = 4
        self.commitOnChange = 0

        ## GUI
        box = OWGUI.widgetBox(self.controlArea, "Info", addSpace=True)
        self.infoLabel = OWGUI.label(box, self, "Chemicals:")
        box = OWGUI.radioButtonsInBox(
            self.controlArea, self, "showFragments",
            ["Show molecules", "Show fragments"], "Show",
            callback=self.updateitems
        )

        self.showFragmentsRadioButton = box.buttons[-1]
        self.markFragmentsCheckBox = OWGUI.checkBox(
            box, self, "colorFragments", "Mark fragments",
            callback=self._update
        )
        box.setSizePolicy(
            QSizePolicy(QSizePolicy.Minimum, QSizePolicy.Maximum))
        OWGUI.separator(self.controlArea)

        self.moleculeSmilesCombo = OWGUI.comboBox(
            self.controlArea, self, "smiles_var",
            "Molecule SMILES Attribute",
            callback=self.updateitems
        )
        self.moleculeSmilesCombo.box.setSizePolicy(
            QSizePolicy(QSizePolicy.Minimum, QSizePolicy.Maximum)
        )
        self.smiles_var_model = VariableListModel(parent=self)
        self.moleculeSmilesCombo.setModel(self.smiles_var_model)

        OWGUI.separator(self.controlArea)
        box = OWGUI.widgetBox(self.controlArea, "Molecule Title Attributes",
                              addSpace=True)

        self.title_var_view = QListView(
            selectionMode=QListView.ExtendedSelection
        )
        self.title_var_model = VariableListModel(parent=self)
        self.title_var_view.setModel(self.title_var_model)
        self.title_var_view.selectionModel().selectionChanged.connect(
            self._title_selection_changed
        )
        box.layout().addWidget(self.title_var_view)

        OWGUI.separator(self.controlArea)
        self.fragmentSmilesCombo = OWGUI.comboBox(
            self.controlArea, self, "fragmentSmilesAttr",
            "Fragment SMILES Attribute",
            callback=self.updateFragmentsListBox
        )

        self.fragmentSmilesCombo.setModel(VariableListModel(parent=self))
        self.fragmentSmilesCombo.box.setSizePolicy(
            QSizePolicy(QSizePolicy.Minimum, QSizePolicy.Maximum)
        )
        OWGUI.separator(self.controlArea)
        box = OWGUI.spin(self.controlArea, self, "imageSize", 50, 500, 10,
                         box="Image Size", callback=self._image_size_changed)

        box.setSizePolicy(
            QSizePolicy(QSizePolicy.Minimum, QSizePolicy.Maximum))

        OWGUI.separator(self.controlArea)
        box = OWGUI.widgetBox(self.controlArea, "Selection", addSpace=True)
        OWGUI.checkBox(box, self, "commitOnChange", "Commit on change")

        self.selectMarkedMoleculesButton = OWGUI.button(
            box, self, "Select &matched molecules", self.select_marked
        )
        OWGUI.button(box, self, "&Commit", callback=self.commit, default=True)
        OWGUI.separator(self.controlArea)
        OWGUI.rubber(self.controlArea)

        spliter = QSplitter(Qt.Vertical)
        self.scrollArea = ScrollArea(spliter)

        self.grid = GridWidget()
        self.grid.selectionChanged.connect(self._on_selection_changed)

        self.scrollArea.setWidget(self.grid)
        self.scrollArea.setWidgetResizable(True)
        self.mainArea.layout().addWidget(spliter)

        if pybel:
            self.listBox = QListWidget(spliter)
        else:
            self.listBox = QListWidget(None)
            self.listBox.setHidden(True)

        self.listBox.itemClicked.connect(self.fragmentSelection)

        self.fragmentSmilesCombo.box.setDisabled(not pybel)

        self.data = None
        self.data_subset = []
        self.fragment_data = None
        self.resize(800, 600)
        self.listBox.setMaximumHeight(150)
        self.fragmentSmilesCombo.setDisabled(True)
        self.selectMarkedMoleculesButton.setDisabled(True)
        self.markFragmentsCheckBox.setDisabled(True)
        self.showFragmentsRadioButton.setDisabled(True)

        self.loadSettings()

        if not pybel:
            self.showFragments = 0
            self.warning(10,
                         "Pybel module not installed. To view molecule fragments\n"
                         "please install openbabel python extension.")

        self.__loop = None

    def setMoleculeTable(self, data):
        self.closeContext()
        self.clear()

        self.data = data
        if data is not None:
            all_vars = data.domain.variables + data.domain.get_metas().values()
            text_vars = filter(
                lambda v: isinstance(v, (Orange.feature.Discrete,
                                         Orange.feature.String)),
                all_vars)
            var_scored = score_smiles_variables(data, text_vars)
            self.smiles_var_model[:] = [var for var, _ in var_scored]
            self.smiles_var = max(range(len(var_scored)),
                                  key=lambda i: var_scored[i][1])
            self.title_var_model[:] = all_vars

            self.setFragmentSmilesCombo()
            self.updateFragmentsListBox()
            if self.data_subset:
                try:
                    self.data_subset = self.data_subset.select(self.data.domain)
                except Exception:
                    self.data_subset = []
            self.openContext("", data)
        else:
            self.defFragmentSmiles = []
            if not self.fragmentSmilesAttr:
                self.listBox.clear()

            self.openContext("", data)
            self.send("Selected Molecules", None)

    def setMoleculeSubset(self, data):
        self.data_subset = data
        try:
            self.data_subset = self.data_subset.select(self.data.domain)
        except Exception:
            self.data_subset = []

    def setFragmentTable(self, data):
        self.fragment_data = data
        if data is not None:
            self.setFragmentSmilesCombo()
            self.updateFragmentsListBox()
            self.selectedFragment = ""
        else:
            self.setFragmentSmilesCombo()
            self.updateFragmentsListBox()

        self.fragmentSmilesCombo.setEnabled(data is not None)

    def handleNewSignals(self):
        self.updateitems()

    def clear(self):
        self.smiles_var_model[:] = []
        self.title_var_model[:] = []

        self.fragmentSmilesCombo.clear()
        self.grid.clear()
        self._widgets = []
        self._items = []

        if self.__loop is not None:
            self.__loop.close()
            self.__loop = None

    def cleargrid(self):
        self.grid.clear()
        self._widgets = []

    def _update_titles(self):
        if self.data is None:
            return

        title_vars = [self.title_var_model[ind.row()]
                      for ind in self.title_var_view.selectedIndexes()]

        for item, widget in zip(self._items, self._widgets):
            inst = self.data[item.index]
            text = " / ".join(map(str, (inst[var] for var in title_vars)))
            widget.label.setText(text)

    def setFragmentSmilesCombo(self):
        if self.fragment_data:
            candidates = score_smiles_variables(self.fragment_data)
        else:
            candidates = []

        self.fragmentSmilesCombo.model()[:] = [v for v, _ in candidates]

        if self.fragmentSmilesAttr > len(candidates):
            self.fragmentSmilesAttr = 0

    def updateFragmentsListBox(self):
        if pybel is None:
            return

        fragvars = self.fragmentSmilesCombo.model()
        if 0 <= self.fragmentSmilesAttr < len(fragvars):
            fvar = fragvars[self.fragmentSmilesAttr]
        else:
            fvar = None

        if fvar:
            frags = [str(e[fvar]) for e in self.fragment_data
                     if not e[fvar].is_special()]
            self.fragmentSmiles = [""] + frags
        else:
            self.fragmentSmiles = [""] + self.defFragmentSmiles

        self.listBox.clear()
        self.listBox.addItems(self.fragmentSmiles)

        self.showFragmentsRadioButton.setDisabled(len(self.fragmentSmiles) == 1)
        self.markFragmentsCheckBox.setDisabled(len(self.fragmentSmiles) == 1)
        self.selectMarkedMoleculesButton.setDisabled(True)

    def fragmentSelection(self, item):
        if pybel is None:
            return

        index = self.listBox.indexFromItem(item).row()
        if index == -1:
            index = 0
        self.selectedFragment = self.fragmentSmiles[index]
        self.selectMarkedMoleculesButton.setEnabled(bool(self.selectedFragment))
        self.markFragmentsCheckBox.setEnabled(bool(self.selectedFragment))
        if not self.showFragments and self.colorFragments:
            self._update()

    def _title_text(self, index):
        title_vars = [self.title_var_model[ind.row()]
                      for ind in self.title_var_view.selectedIndexes()]
        inst = self.data[index]
        return " / ".join(map(str, (inst[var] for var in title_vars)))

    def _items_from_var(self, var):
        if self.data is None:
            return None

        values = [(i, str(inst[var])) for i, inst in enumerate(self.data)
                  if not inst[var].is_special()]
        return [Item(i, smiles, *self._parse_smiles(smiles))
                for i, smiles in values]

    def _parse_smiles(self, smiles):
        try:
            return (OK, molecule_from_smiles(smiles))
        except Exception:
            return (ParseError, None)

    def updateitems(self):
        if self.showFragments and self.fragmentSmiles:
            values = [(None, frag) for frag in self.fragmentSmiles[1:]]
            items = [Item(i, smiles, *self._parse_smiles(smiles))
                     for i, smiles in values]
        else:
            smilesvar = self.smiles_var_model[self.smiles_var]
            items = self._items_from_var(smilesvar)

        self._items = items
        self.setupgrid()

    def setupgrid(self):
        self.cleargrid()

        layout = self.grid
        widgets = []

        for item in self._items:
            thumb = ThumbnailWidget(self.grid)
            thumb.setImageSize(self.imageSize, self.imageSize)
            if item.index is not None:
                text = self._title_text(item.index)
            else:
                text = ""
            thumb.label.setText(text)

            widgets.append(thumb)
            layout.appendWidget(thumb)

        self._widgets = widgets
        self.infoLabel.setText("Chemicals %i" % len(self._items))
        self._update()

    def __update_items(self, items, widgets, pattern=None):
        for i, item, widget in zip(range(len(items)), items, widgets):
            if item.status != ParseError:
                if pattern is not None:
                    emb = substructure_embedding(item.molecule, pattern)
                    emb = reduce(list.__iadd__, emb, [])
                    svg = molecule_to_svg_with_substructure(item.molecule, emb)
                else:
                    svg = molecule_to_svg(item.molecule)
            else:
                svg = ""

            widget.setData(svg)
            widget.setEnabled(True)
            yield i * 100.0 / len(items)

    def _update(self):
        if self.showFragments and self.fragmentSmiles:
            loop = self.__update_items(self._items, self._widgets)
        elif self.colorFragments and self.selectedFragment:
            pattern = pybel.Smarts(self.selectedFragment)
            loop = self.__update_items(self._items, self._widgets, pattern)
        else:
            loop = self.__update_items(self._items, self._widgets)
        self.__schedule(loop)

    def __schedule(self, coroutine):
        if self.__loop is not None:
            self.progressBarFinished()
            self.__loop.close()
            self.__loop = None

        self.__loop = coroutine

        self.progressBarInit()
        QTimer.singleShot(0, self.__loop_update)

    @Slot()
    def __loop_update(self):
        if self.__loop is None:
            return

        try:
            progress = next(self.__loop)
        except StopIteration:
            self.__loop = None
            self.progressBarFinished()
        else:
            self.progressBarSet(progress)
            QTimer.singleShot(0, self.__loop_update)

    def _title_selection_changed(self):
        self._update_titles()

    def _image_size_changed(self):
        for widget in self._widgets:
            widget.setImageSize(self.imageSize, self.imageSize)

        self.grid.layout().invalidate()

    def select_marked(self):
        if not pybel:
            return

        if not self.showFragments:
            pattern = pybel.Smarts(self.selectedFragment)
            for item, widget in zip(self._items, self._widgets):
                if item.status != ParseError:
                    emb = substructure_embedding(item.molecule, pattern)
                    widget.setSelected(bool(emb))
                else:
                    widget.setSelected(False)

            if self.commitOnChange:
                self.commit()

    def _on_selection_changed(self):
        if self.commitOnChange:
            self.commit()

    def commit(self):
        if self.showFragments:
            svar = self.smiles_var_model[self.smiles_var]
            items = self._items_from_var(svar)
            frags = [item for item, w in zip(self._items, self._widgets)
                     if w.selected]
            patterns = [pybel.Smarts(item.smiles) for item in frags]

            def test(molecule, patterns):
                return any(bool(substructure_embedding(molecule, patt))
                           for patt in patterns)

            matched = filter(
                lambda item: item.status != ParseError and
                             test(item.molecule, patterns),
                items
            )
            instances = [self.data[item.index] for item in matched]

            if instances:
                table = Orange.data.Table(instances)
                self.send("Selected Molecules", table)
            else:
                self.send("Selected Molecules", None)
        else:
            items = [item for item, w in zip(self._items, self._widgets)
                     if w.selected]
            instances = [self.data[item.index] for item in items]

            if instances:
                table = Orange.data.Table(instances)
                self.send("Selected Molecules", table)
            else:
                self.send("Selected Molecules", None)

    def onDeleteWidget(self):
        OWWidget.onDeleteWidget(self)
        if self.__loop is not None:
            self.__loop.close()
            self.__loop = None


def score_smiles_variables(data, candidates=None, max_sample=20):
    if candidates is None:
        candidates = data.domain.variables + data.domain.getmetas().values()
        candidates = filter(lambda v: isinstance(v, (Orange.feature.Discrete,
                                                     Orange.feature.String)),
                            candidates)

    if len(data) > max_sample:
        indices = Orange.data.sample.SubsetIndices2(data, max_sample)
        data = data.select(indices, 0)

    scored = [(var, sum([is_valid_smiles(str(e[var])) for e in data]))
              for var in candidates]
    return scored


def main():
    import os
    from PyQt4.QtGui import QApplication
    app = QApplication([])
    w = OWMoleculeVisualizer()
    data = Orange.data.Table("Chemogenomics")
    w.setMoleculeTable(Orange.data.Table(data[:10]))
    frag = Orange.data.Table(os.path.expanduser("~/Documents/Datasets/fragments.tab"))
    w.setFragmentTable(frag)
    w.handleNewSignals()
    w.show()
    w.raise_()
    return app.exec_()


if __name__ == "__main__":
    main()
