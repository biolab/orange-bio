"""<name>Gene Mania</name>
"""

from OWWidget import *
import OWGUI
import orngNetwork
import obiGeneMania

import os, sys
import multiprocessing
import random
        

class OWGeneMania(OWWidget):
    contextHandlers = {"": DomainContextHandler("", ["selectedOrganismIndex", "selectedGeneAttrIndex", "genesInColumns"])}
    settingsList = ["selectedOrganismIndex", "selectedGeneAttrIndex", "genesInColumns", "selectedMethodIndex", "resultCount"]
    def __init__(self, parent=None, signalManager=None, name="Gene Mania"):
        OWWidget.__init__(self, parent, signalManager, name, wantMainArea=False)
        
        self.inputs = [("Input Genes", ExampleTable, self.setData)]
        self.outputs = [("Network", orngNetwork.Network, Default), ("Items", ExampleTable)]
        
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
        
        self.info = OWGUI.widgetLabel(OWGUI.widgetBox(self.controlArea, "Info", addSpace=True), "No genes on input.\n")
        self.organismCombo = OWGUI.comboBox(self.controlArea, self, "selectedOrganismIndex", "Organims",
                                            items=[t[0] for t in self.organisms],
                                            tooltip="Select the organism",)
#                                            callback=self.onOrganismSelecion)
        
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
#                       toolTips=toolTips,
                       )
        
        OWGUI.spin(self.controlArea, self, "resultCount", 1, 100, 
                   box="Number of gene results")
        
        
        OWGUI.button(self.controlArea, self, "Retrieve", callback=self.retrieve)
        OWGUI.rubber(self.controlArea)
        
        self.resize(100, 200)


    def setData(self, data=None):
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
        
    def updateInfo(self):
        if self.data is not None:
            genes = self.getSelectedGenes()
            self.info.setText("%i genes on input.\n" % len(genes))
        else:
            self.info.setText("No data on input.\n")
        
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
        conn = obiGeneMania.Connection(self.serverAddress)
        
        call = self.asyncCall(conn.retrieveXML, (org, genes), {"m": method, "r": self.resultCount})
        call()

        self.error(0)
        self.progressBarInit()
        self.setEnabled(False)
        try:
#            net = conn.retrieve(org, genes, m=method, r=self.resultCount)
            xml = call.get_result(processEvents=True)
            # Parse the xml in the main thread (pyexpat frequently crashes in 
            # Qt's threading model)
            dom = obiGeneMania.minidom.parseString(xml)
            net = obiGeneMania.parse(dom)
            items = net.items
        except Exception, ex:
            self.error(0, "Failed to retrieve network from server!\n  " + str(ex))
            sys.excepthook(*sys.exc_info())
            net = None
            items = None
        self.setEnabled(True)
        self.progressBarFinished()
        
        self.send("Network", net)
        self.send("Items", items)
        
        
if __name__ == "__main__":
    app = QApplication(sys.argv)
    w = OWGeneMania()
    w.show()
    data = orange.ExampleTable("../../../doc/datasets/brown-selected.tab")
    w.setData(data)
    app.exec_()
    w.saveSettings()
        
        
