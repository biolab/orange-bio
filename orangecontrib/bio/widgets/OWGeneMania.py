"""<name>GeneMANIA</name>
<icon>icons/GeneMANIA.svg</icon>
"""

from __future__ import absolute_import

import os, sys
import multiprocessing
import random

import Orange
from Orange.OrangeWidgets import OWGUI
from Orange.OrangeWidgets.OWWidget import *

from .. import obiGeneMania

NAME = "GeneMANIA"
DESCRIPTION = ""
ICON = "icons/GeneMANIA.svg"

INPUTS = [("Input Genes", Orange.data.Table, "setData")]
OUTPUTS = [("Network", Orange.network.Graph, Default),
           ("Items", Orange.data.Table)]

REPLACES = ["_bioinformatics.widgets.OWGeneMania.OWGeneMania"]


class BarItemDelegate(QStyledItemDelegate):
    BarRole = OWGUI.OrangeUserRole.next() 
    BarForegroundRole = OWGUI.OrangeUserRole.next()
    
    def __init__(self, parent, brush=QBrush(QColor(255, 170, 127)), scale=(0.0, 1.0)):
        QStyledItemDelegate.__init__(self, parent) 
        self.brush = brush
        self.scale = scale
        
        
    def paint(self, painter, option, index):
        painter.save()
        qApp.style().drawPrimitive(QStyle.PE_PanelItemViewRow, option, painter)
        qApp.style().drawPrimitive(QStyle.PE_PanelItemViewItem, option, painter)
        rect = option.rect
        try:
            val, ok = index.data(self.BarRole).toDouble()
            if ok:
                color = index.data(self.BarForegroundRole)
                if color.isValid() and color.type() == QVariant.Color:
                    brush = QBrush(color)
                else:
                    brush = self.brush
                    
                minval, maxval = self.scale
                val = (val - minval) / (maxval - minval)
                painter.save()
                if option.state & QStyle.State_Selected:
                    painter.setOpacity(0.75)
                painter.setBrush(brush)
                painter.drawRect(rect.adjusted(1, 1, max(-rect.width() * (1.0 - val) - 2, - rect.width() + 2), -2))
                painter.restore()
        except Exception, ex:
            print >> sys.stderr, ex
            
        text = index.data(Qt.DisplayRole).toString()
        if text:
            align, ok = index.data(Qt.TextAlignmentRole).toInt()
            if not ok:
                align = Qt.AlignVCenter | Qt.AlignLeft
                
            painter.drawText(option.rect, align, text)
        painter.restore()
        
        
    def sizeHint(self, option, index):
        size = QStyledItemDelegate.sizeHint(self, option, index)
        metrics = QFontMetrics(option.font)
        height = metrics.lineSpacing() + 2
        return QSize(size.width(), height)
    

class OWGeneMania(OWWidget):
    contextHandlers = {"": DomainContextHandler("", ["selectedOrganismIndex", "selectedGeneAttrIndex", "genesInColumns"])}
    settingsList = ["serverAddress", "selectedOrganismIndex", "selectedGeneAttrIndex", "genesInColumns", "selectedMethodIndex", "resultCount"]
    def __init__(self, parent=None, signalManager=None, name="GeneMANIA"):
        OWWidget.__init__(self, parent, signalManager, name, wantMainArea=True)
        
        self.inputs = [("Input Genes", ExampleTable, self.setData)]
        self.outputs = [("Network", Orange.network.Graph, Default), ("Items", ExampleTable)]
        
        self.serverAddress = obiGeneMania.DEFAULT_SERVER
        self.selectedOrganismIndex = 0
        self.selectedGeneAttrIndex = 0
        self.genesInColumns = False
        self.selectedMethodIndex = 1
        self.resultCount = 10
        
        self.data = None
        
        self.loadSettings()
        
        #####
        # GUI
        #####
        
        self.organisms = [("H. sapiens", "9606"),
                          ("A. thaliana", "3702"),
                          ("C. elegans", "6239"),
                          ("D. melanogaster", "7227"),
                          ("M. musculus", "10090"),
                          ("S. cerevisiae", "4932")]
        
        box = OWGUI.widgetBox(self.controlArea, "Info", addSpace=True)
        self.info = OWGUI.widgetLabel(box, "No genes on input.")
        self.infoState = OWGUI.widgetLabel(box, "")
        self.infoState.setWordWrap(True)
        
        box = OWGUI.widgetBox(self.controlArea, "GeneMANIA server address", addSpace=True)
        OWGUI.lineEdit(box, self, "serverAddress")
        
        self.organismCombo = OWGUI.comboBox(self.controlArea, self, "selectedOrganismIndex", "Organims",
                                            items=[t[0] for t in self.organisms],
                                            tooltip="Select the organism",
                                            callback=self.updateInfo)
        
        box = OWGUI.widgetBox(self.controlArea, "Genes", addSpace=True)
        self.geneAttrCombo = OWGUI.comboBox(box, self, "selectedGeneAttrIndex",
                                            tooltip="Select the attribute with gene names",
                                            callback=self.updateInfo)
        
        cb = OWGUI.checkBox(box, self, "genesInColumns", "Use attribute names",
                            tooltip="Use attribute names as gene names instead.",
                            callback=self.updateInfo,
                            disables=[(-1, self.geneAttrCombo)])
        cb.makeConsistent()
        self.methodItems = [("Automatic relevance", "automatic_relevance"),
                            ("Automatic", "automatic"),
                            ("Biological process", "bp"),
                            ("Molecular function", "mf"),
                            ("Cellular component", "cc"),
                            ("Average", "average"),
                            ("Average category", "average_category")]
        
        toolTips = ["Assigned based on query genes",
                    "Automatically selected weighting method",
                    "Biological process based",
                    "Molecular function based",
                    "Cellular component based",
                    "Equal by data type",
                    "Equal by network"]
        
        OWGUI.comboBox(self.controlArea, self, "selectedMethodIndex", 
                       box="Net combining method",
                       items=[t[0] for t in self.methodItems],
                       callback=self.updateInfo
                       )
        
        OWGUI.spin(self.controlArea, self, "resultCount", 1, 100, 
                   box="Number of gene results",
                   callback=self.updateInfo
                   )
        
        self.geneManiaLinkLabel = OWGUI.widgetLabel(self.controlArea, "")
        self.geneManiaLinkLabel.setOpenExternalLinks(True)
        OWGUI.button(self.controlArea, self, "Retrieve", callback=self.retrieve)
        OWGUI.rubber(self.controlArea)
        
        self.networksReport = QTreeView()
        self.networksReport.setEditTriggers(QTreeView.NoEditTriggers)
        box = OWGUI.widgetBox(self.mainArea, "Networks")
        box.layout().addWidget(self.networksReport)
        
        self.resize(100, 200)
        
        self.connect(self, SIGNAL("widgetStateChanged(QString, int, QString)"), self.updateInfo)


    def setData(self, data=None):
        self.error([0,1,2])
        self.warning([0])
        self.data = data
        self.closeContext("")
        self.geneAttrCombo.clear()
        self.candidateGeneAttrs = []
        self.selectedGeneAttrIndex = 0
        
        if data is not None:
            self.candidateGeneAttrs = data.domain.variables + data.domain.getmetas().values()
            self.candidateGeneAttrs = [attr for attr in self.candidateGeneAttrs if attr.varType != orange.VarTypes.Continuous]
            self.geneAttrCombo.addItems([attr.name for attr in self.candidateGeneAttrs])
            self.openContext("", data)
        else:
            self.send("Network", None)
            self.send("Items", None)
            
        self.updateInfo()
            
    def onGeneSourceSelection(self):
        genes = self.getSelectedGenes()
        self.info.setText("")
        
    def updateInfo(self, *args):
        if self.data is not None:
            genes = self.getSelectedGenes()
            htmlState = self.widgetStateToHtml()
            self.info.setText("%i genes on input." % len(genes))
            self.infoState.setText(htmlState)
        else:
            self.info.setText("No data on input.")
            self.infoState.setText("")
        
        if self.data:
            org = self.organisms[self.selectedOrganismIndex][1]
            genes = self.getSelectedGenes()
            method = self.methodItems[self.selectedMethodIndex][1]
            conn = obiGeneMania.Connection(self.serverAddress)
            self.geneManiaLinkLabel.setText('<a href="%s">View network in external browser</a>' % conn._queryPage(org, genes, method, self.resultCount))
        else:
            self.geneManiaLinkLabel.setText('')
        
    def getSelectedGenes(self):
        if self.data is not None:
            if self.genesInColumns:
                names = [attr.name for attr in self.data.domain.attributes]
            else:
                attr = self.candidateGeneAttrs[self.selectedGeneAttrIndex]
                names = set([str(ex[attr]) for ex in self.data if not ex[attr].isSpecial()])
            return names
        else:
            return []
        
        
    def retrieve(self):
        org = self.organisms[self.selectedOrganismIndex][1]
        genes = self.getSelectedGenes()
        method = self.methodItems[self.selectedMethodIndex][1]
        self.error([0, 1, 2])
        self.warning([0])
        
        conn = obiGeneMania.Connection(self.serverAddress)
        errorCode, invalid, genes = conn.validate(org, genes)
        if not genes:
            self.error(2, "No valid gene names!")
            self.net, self.netTab = None, None
            self.updateNetworksReport()
            return
        elif invalid:
            self.warning(0, "There are invalid gene names on input:\n%s" % (",".join(invalid[:5])) + (" ..." if len(invalid) > 5 else ""))
            
#        print conn._queryPage(org, genes, method, self.resultCount)
        call = self.asyncCall(conn.retrieveXML, (org, genes), {"m": method, "r": self.resultCount})
        call()
        
        self.progressBarInit()
        self.setEnabled(False)
        try:
#            net = conn.retrieve(org, genes, m=method, r=self.resultCount)
            xml = call.get_result(processEvents=True)
            # Parse the xml in the main thread (pyexpat frequently crashes in 
            # Qt's threading model)
            dom = obiGeneMania.minidom.parseString(xml)
            net = obiGeneMania.parse(dom)
            items = net.items()
        except Exception, ex:
            self.error(0, "Failed to retrieve network from server!\n" + str(ex))
            sys.excepthook(*sys.exc_info())
            net = None
            items = None
        
        if net:
            try:
                self.netTab = obiGeneMania.parsePage(conn._page)
            except Exception:
                self.error(1, "Failed to parse network tabs!\n" + str(ex))
                self.netTab = None
        else:
            self.netTab = None
        self.setEnabled(True)
        self.progressBarFinished()
        
        self.net = net
        
        self.send("Network", net)
        self.send("Items", items)
        
        self.updateNetworksReport()
        
    def updateNetworksReport(self):
        model = QStandardItemModel(self)
        model.setHorizontalHeaderLabels(["Networks", "Weight"])
        root = model.invisibleRootItem()
        def toFloat(s):
            if s.strip().endswith("%"):
                return float(s.strip()[:-1])
            else:
                return float(s)
            
        if self.netTab and self.net and self.net.links:
            for group in self.netTab:
                groupItem = QStandardItem(group.name)
                groupWeightItem = QStandardItem("%.2f %%" % toFloat(group.weight))
                groupWeightItem.setData(QVariant(toFloat(group.weight) / 100.0), BarItemDelegate.BarRole)
                root.appendRow([groupItem, groupWeightItem])
                for net in group.networks:
                    netItem = QStandardItem(net.name)
                    netItem.setData(QVariant(net.description), Qt.ToolTipRole)
                    
                    netWeightItem = QStandardItem("%.2f %%" % toFloat(net.weight))
                    netWeightItem.setData(QVariant(toFloat(net.weight) / 100.0), BarItemDelegate.BarRole)
                    netWeightItem.setData(QVariant(net.description), Qt.ToolTipRole)
                    groupItem.appendRow([netItem, netWeightItem])
                
        self.networksReport.setModel(model)
        self.networksReport.setItemDelegateForColumn(1, BarItemDelegate(self))
        self.networksReport.resizeColumnToContents(0)
            
        
if __name__ == "__main__":
    app = QApplication(sys.argv)
    w = OWGeneMania()
    w.show()
    data = orange.ExampleTable("../../../doc/datasets/brown-selected.tab")
    w.setData(orange.ExampleTable(data[:3]))
    app.exec_()
    w.saveSettings()
        
        
