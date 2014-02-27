"""
<name>MeSH Browser</name>
<description>Browse MeSH ontology.</description>
<icon>icons/MeSHBrowser.svg</icon>
<contact>Crt Gorup (crt.gorup@gmail.com)</contact> 
<priority>2040</priority>
"""

from __future__ import absolute_import

from PyQt4.QtCore import *
from PyQt4.QtGui import *

from Orange.OrangeWidgets import OWGUI
from Orange.OrangeWidgets.OWWidget import *

from ..obiMeSH import *

NAME = "MeSH Browser"
DESCRIPTION = "Browse MeSH ontology."
ICON = "icons/MeSHBrowser.svg"
PRIORITY = 2040

INPUTS = [("Reference data", Orange.data.Table, "getReferenceData"),
          ("Cluster data", Orange.data.Table, "getClusterData")]

OUTPUTS = [("Selected examples", Orange.data.Table)]

REPLACES = ["_bioinformatics.widgets.OWMeSHBrowser.OWMeSHBrowser"]


class MyQTableWidgetItem(QTableWidgetItem):

    """ Our implementation of QTable item allowing numerical sorting.  """

    def __init__(self, table, text):
        QTableWidgetItem.__init__(self, table, QTableWidgetItem.Never, text)
        self.data = text

    def key(self):	  # additional setting to correctly handle text and numerical sorting
        try:	# sorting numerical column
            e = float(self.data)
            tdata = pval = "%.4g" % e
            l = len(self.data)
            pre = ""

            offset = 0	# additional parameter to hadle exponent '2e-14' float format
            if(tdata.count('e') > 0 or tdata.count('E') > 0):
                pre = "*"
                offset = 1

            for i in range(0, 40 - l - offset):
                pre = pre + "0"
            return pre + tdata
        except ValueError:	  # sorting text column
            return self.data


class ListViewToolTip(QToolTip):

    """brief A class to allow tooltips in a listview."""

    def __init__(self, view, column, data):
        """brief ListViewToolTip constructor.
        \param view		  QListView instance
        \param column		Listview column
        \param truncatedOnly Only display the tooltip if the column data is truncated
        """
        QToolTip.__init__(self, view.viewport())
        self.__view = view
        self.__col = column
        self.__data = data	  # mapping from name -> description
        # self.setWakeUpDelay(400)

    def appendData(self, data):
        self.__data = data

    def maybeTip(self, pos):
        """brief Draw the tooltip.
        \param pos Tooltip position.
        """
        item = self.__view.itemAt(pos)
        if item is not None:
            if(self.__data.has_key(str(item.text(self.__col)))):
                tipString = self.__data[str(item.text(self.__col))]
                counter = 45
                newTipString = ""
                for i in tipString.split(" "):
                    if counter < 0:
                        newTipString = newTipString + i + "\n"
                        counter = 45
                    else:
                        newTipString = newTipString + i + " "
                        counter -= len(i)
                tipString = newTipString
                cr = self.__view.itemRect(item)
                headerPos = self.__view.header().sectionPos(self.__col)
                cr.setLeft(headerPos)
                cr.setRight(headerPos + self.__view.header().sectionSize(self.__col))
                self.tip(cr, tipString)
            else:
                print item.text(self.__col) + " not found in toDecs"


class OWMeSHBrowser(OWWidget):
    settingsList = ["multi", "maxPValue", "minExamplesInTerm"]

    def __init__(self, parent=None, signalManager=None):
        OWWidget.__init__(self, parent, signalManager, "MeshBrowser")
        self.inputs = [("Reference data", ExampleTable, self.getReferenceData), ("Cluster data", ExampleTable, self.getClusterData)]
        self.outputs = [("Selected examples", ExampleTable)]

        # widget variables
        self.loadedRef = 0
        self.loadedClu = 0
        self.maxPValue = 0.05
        self.minExamplesInTerm = 5
        self.multi = 1
        self.reference = None
        self.cluster = None
        self.loadSettings()
        self.mesh = obiMeSH() # main object is created
        self.dataLoaded = self.mesh.dataLoaded

        # left pane
        box = OWGUI.widgetBox(self.controlArea, "Info")
        #box = QGroupBox("Info", self.controlArea)
        self.infoa = OWGUI.label(box, self, "No reference data.")
        self.infob = OWGUI.label(box, self, "No cluster data.")
        self.ratio = OWGUI.label(box, self, "")
        self.ref_att = OWGUI.label(box, self, "")
        self.clu_att = OWGUI.label(box, self, "")
        self.resize(960, 600)
        OWGUI.separator(self.controlArea)

        self.optionsBox = OWGUI.widgetBox(self.controlArea, "Options")
        self.maxp = OWGUI.lineEdit(self.optionsBox, self, "maxPValue", label="threshold:", orientation="horizontal", labelWidth=120, valueType=float)
        self.minf = OWGUI.lineEdit(self.optionsBox, self, "minExamplesInTerm", label="min. frequency:", orientation="horizontal", labelWidth=120, valueType=int)
        #OWGUI.checkBox(self.optionsBox, self, 'multi', 'Multiple selection', callback= self.checkClicked)
        OWGUI.button(self.optionsBox, self, "Refresh", callback=self.refresh)

        # right pane
        self.col_size = [280, 84, 84, 100, 110]
        self.sort_col = 0
        self.sort_dir = True
        self.columns = ['MeSH term', '# reference', '# cluster', 'p value', 'fold enrichment'] # both datasets

        self.splitter = QSplitter(Qt.Vertical, self.mainArea)
        self.mainArea.layout().addWidget(self.splitter)

        # list view
        self.meshLV = QTreeWidget(self.splitter)
        #self.meshLV.setSelectionMode(QAbstractItemView.MultiSelection)
        self.meshLV.setAllColumnsShowFocus(1)
        self.meshLV.setColumnCount(len(self.columns))
        self.meshLV.setHeaderLabels(self.columns)

        self.meshLV.header().setClickable(True)
        #self.meshLV.header().setSortIndicatorShown(True)
        #self.meshLV.setSortingEnabled(True)
        self.meshLV.setRootIsDecorated(True)
        self.connect(self.meshLV, SIGNAL("itemSelectionChanged()"), self.viewSelectionChanged)
        #self.meshLV.setItemDelegateForColumn(3, EnrichmentColumnItemDelegate(self))
        #self.tooltips = ListViewToolTip(self.meshLV,0, self.mesh.toDesc)

        # table of significant mesh terms
        self.sigTermsTable = QTableWidget(self.splitter)
        self.sigTermsTable.setColumnCount(len(self.columns))
        self.sigTermsTable.setRowCount(4)
        ## hide the vertical header
        self.sigTermsTable.verticalHeader().hide()
        #self.sigTermsTable.setLeftMargin(0)
        #self.sigTermsTable.setSelectionMode(QAbstractItemView.MultiSelection)

        for i in range(0, len(self.columns)):
            self.sigTermsTable.horizontalHeader().resizeSection(i, self.col_size[i])
            self.meshLV.header().resizeSection(i, self.col_size[i])

        self.sigTermsTable.setHorizontalHeaderLabels(self.columns)

        self.connect(self.sigTermsTable, SIGNAL("itemSelectionChanged()"), self.tableSelectionChanged)
        self.connect(self.sigTermsTable, SIGNAL("clicked(int,int,int,const QPoint&)"), self.tableClicked)
        self.splitter.show()
        self.optionsBox.setDisabled(1)

    def tableSelectionChanged(self):
        return True

    def tableClicked(self, row, col, button, point):
        if self.sort_col == col:
            self.sort_dir = not self.sort_dir
        else:
            self.sort_col = col
            self.sort_dir = True

        self.sigTermsTable.sortItems(self.sort_col, Qt.DescendingOrder)
        #print "sortiram ", col, " ",row

    def checkClicked(self):
        if self.multi == 0:
            self.meshLV.clearSelection()
        self.meshLV.setSelectionMode(QAbstractItemView.MultiSelection)

    def viewSelectionChanged(self):
        """
        Function viewSelectionChanged is used to handle the widget output once the user clicks or selects MeSH term inside the tree view.
        """
        items = list()
        self.progressBarInit()

        itms = self.meshLV.selectedItems()
        for i in itms:
            items.append(i.term)

        #for i in self.lvItem2Mesh.iterkeys():
        #	if i.isSelected():
        #		items.append(self.lvItem2Mesh[i][1])
        #print "selecting ", items

        if self.reference:
            data = self.mesh.findSubset(self.reference, items, callback=self.progressBarSet, MeSHtype='term')
        else:
            data = self.mesh.findSubset(self.cluster, items, callback=self.progressBarSet, MeSHtype='term')

        #print items
        self.send("Selected examples", data)
        self.progressBarFinished()

    def __updateData__(self):
        """
        Function __updateData__ is used to display the results of the MeSH term enrichment analysis inside the widget components.
        """
        self.lvItem2Mesh = dict()

        if(self.reference and self.cluster):
            if(len(self.cluster) > len(self.reference)):
                self.optionsBox.setDisabled(1)
                return False
            # everything is ok, now we can calculate enrichment and update labels, tree view and table data
            #self.warning()
            self.optionsBox.setDisabled(0)
            self.progressBarInit()
            self.treeInfo, self.results = self.mesh.findEnrichedTerms(self.reference, self.cluster, self.maxPValue, treeData= True, callback= self.progressBarSet)
            self.progressBarFinished()
            self.ratio.setText("ratio = %.4g" % self.mesh.ratio)
            self.clu_att.setText("cluster MeSH att: " + self.mesh.clu_att)
            self.ref_att.setText("reference MeSH att: " + self.mesh.ref_att)

            # table data update
            self.sigTermsTable.setRowCount(len(self.results))
            index = 0
            for i in self.results.iterkeys(): ## sorted by the p value
               # mTerm = i[0]
                mID = self.mesh.toName[i] + " (" + i + ")"
                rF = self.results[i][0]
                cF = self.results[i][1]
                pval = self.results[i][2]
                fold = self.results[i][3]
                pval = "%.4g" % pval
                fold = "%.4g" % fold
                vals = [mID, rF, cF, pval, fold]
                for j in range(len(vals)):
                    self.sigTermsTable.setItem(index, j, QTableWidgetItem(str(vals[j])))
                index = index + 1

            # initial sorting - p value
            self.sort_col = 3
            self.sort_dir = True
            self.sigTermsTable.sortItems(self.sort_col, Qt.DescendingOrder)

            # tree view update - The most beautiful part of this widget!
            starters = self.treeInfo["tops"]	   # we get a list of possible top nodes
            self.meshLV.clear()

            for e in starters:	  # we manualy create top nodes
                f = QTreeWidgetItem(self.meshLV);
                #f.setOpen(1)
                rfr = str(self.results[e][0])
                cfr = str(self.results[e][1])
                pval = "%.4g" % self.results[e][2]
                fold = "%.4g" % self.results[e][3]
                self.lvItem2Mesh[f] = (self.mesh.toName[e], e)
                data = [self.mesh.toName[e], rfr, cfr, pval, fold]
                f.term = self.mesh.toName[e]
                for t in range(len(data)):
                    f.setText(t, data[t])
                self.__treeViewMaker__(f, e, False)
            self.meshLV.expandAll()

        elif self.reference or self.cluster:
            if self.reference:
                current_data = self.reference
            else:
                current_data = self.cluster

            self.optionsBox.setDisabled(0)
            self.progressBarInit()
            self.treeInfo, self.results = self.mesh.findFrequentTerms(current_data, self.minExamplesInTerm, treeData=True, callback=self.progressBarSet)
            self.progressBarFinished()
            if self.reference:
                self.ref_att.setText("reference MeSH att: " + self.mesh.solo_att)
            else:
                self.clu_att.setText("cluster MeSH att: " + self.mesh.solo_att)

            # table data update
            self.sigTermsTable.setRowCount(len(self.results))
            index = 0
            for i in self.results.iterkeys():
                mID = self.mesh.toName[i] + " (" + i + ")"
                rF = self.results[i]
                vals = [mID, rF]
                for j in range(len(vals)):
                    self.sigTermsTable.setItem(index, j, QTableWidgetItem(str(vals[j])))
                index = index + 1

            # initial sorting - frequency
            self.sort_col = 1
            self.sort_dir = True
            self.sigTermsTable.sortItems(self.sort_col, Qt.DescendingOrder)

            # tree view update - The most beautiful part of this widget!
            starters = self.treeInfo["tops"]	   # we get a list of possible top nodes
            self.meshLV.clear()

            for e in starters:	  # we manualy create top nodes
                f = QTreeWidgetItem(self.meshLV);
                #f.setOpen(1)
                self.lvItem2Mesh[f] = (self.mesh.toName[e], e)
                rfr = str(self.results[e])
                data = [self.mesh.toName[e], rfr]
                f.term = self.mesh.toName[e]
                for t in range(len(data)):
                    f.setText(t, data[t])
                self.__treeViewMaker__(f, e, True)
            self.meshLV.expandAll()

    def __treeViewMaker__(self, parentLVI, parentID, soloMode):
        """
        Function __treeViewMaker__ is used to build the tree in treeListView. When soloMode=True function only displays \
        first two columns inside the tree view (suitable when only one dataset is present).
        """
        for i in self.treeInfo[parentID]:   # for each succesor
            f = QTreeWidgetItem(parentLVI);
            f.term = self.mesh.toName[i]
            data = [self.mesh.toName[i]]
            if soloMode:
                rfr = str(self.results[i])
                data.append(rfr)
            else:		   # when we have referece and cluster dataset we have to print additional info
                rfr = str(self.results[i][0])
                cfr = str(self.results[i][1])
                pval = "%.4g" % self.results[i][2]
                fold = "%.4g" % self.results[i][3]
                data.extend([rfr, cfr, pval, fold])
            self.lvItem2Mesh[f] = (data[0], i)		# mapping   QTreeWidgetItem <-> mesh id
            for t in range(len(data)):
                f.setText(t, data[t])
            self.__treeViewMaker__(f, i, soloMode)
        return True

    def refresh(self):
        """
        Function refresh is executed as number of widget inputs changes.
        """
        if self.reference or self.cluster:
            self.__switchGUI__()
            self.__updateData__()

    def __clearGUI__(self):
        """
        Function __clearGUI__ sets tree view, table and labels to their default values.
        """
        self.meshLV.clear()
        self.sigTermsTable.setRowCount(0)

    def __switchGUI__(self):
        """
        Function __switchGUI__ is capable of changing GUI based on number of connected inputs.
        """
        if(len(self.cluster) > len(self.reference)):
            self.optionsBox.setDisabled(1)
            #self.warning("Cluster dataset is greater than reference dataset. Please check the widget inputs.")
            QMessageBox.warning(None, "Invalid input dataset length", "Cluster dataset is longer than the reference dataset. Please check the widget inputs.", QMessageBox.Ok)
            return False
        if not self.reference and not self.cluster:
            self.optionsBox.setDisabled(1)
            return
        self.optionsBox.setDisabled(0)
        solo = True
        if self.reference and self.cluster:
            solo = False
        if solo:
            self.maxp.setDisabled(1)
            self.minf.setDisabled(0)
            for i in range(2, len(self.columns)):
                self.meshLV.hideColumn(i)
                self.sigTermsTable.hideColumn(i)

            self.sigTermsTable.setHorizontalHeaderLabels(["MeSH term", "frequency"])
            self.meshLV.setHeaderLabels(["MeSH term", "frequency"])
            self.ratio.setText("")
        else:
            self.maxp.setDisabled(0)
            self.minf.setDisabled(1)
            for i in range(0, len(self.columns)):
                self.meshLV.showColumn(i)
                self.sigTermsTable.showColumn(i)
            self.sigTermsTable.setHorizontalHeaderLabels(self.columns)
            self.meshLV.setHeaderLabels(self.columns)
            for i in range(0, len(self.columns)):
                self.meshLV.header().resizeSection(i, self.col_size[i])
                self.sigTermsTable.horizontalHeader().resizeSection(i, self.col_size[i])
            self.ratio.setText("ratio = %.4g" % self.mesh.ratio)

    def getReferenceData(self, data):
        """
        Function getReferenceData is executed once reference signal is connected to the widget.
        """
        if data:
            self.reference = data
            self.infoa.setText('%d reference examples' % len(data))
        else:
            self.reference = None
            self.infoa.setText('No reference data.')
            self.ref_att.setText('')

        if self.reference or self.cluster:
            self.__switchGUI__()
            self.__updateData__()
        else:
            self.__clearGUI__()

    def getClusterData(self, data):
        """
        Function getClusterData is executed once cluster signal is connected to the widget.
        """
        if data:
            self.cluster = data
            self.infob.setText('%d cluster examples' % len(data))
        else:
            self.cluster = None
            self.infob.setText('No cluster data.')
            self.clu_att.setText('')

        if self.reference or self.cluster:
            self.__switchGUI__()
            self.__updateData__()
        else:
            self.__clearGUI__()
