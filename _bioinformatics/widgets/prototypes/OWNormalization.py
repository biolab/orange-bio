"""
<name>Normalization</name>
<description>Gene Expression data normalization</description>
<prototype>1</prototype>
"""

from __future__ import absolute_import

import os, sys

import Orange
from Orange.OrangeWidgets import OWGUI
from Orange.OrangeWidgets.OWWidget import *

from ... import obiDifscale

class OWNormalization(OWWidget):
    settingsList = []
    
    METHODS = ["Median", "Quantile"]
    def __init__(self, parent=None, signalManager=None, title="Differentiation Scale"):
        OWWidget.__init__(self, parent, signalManager, title, wantMainArea=False)
        
        self.inputs = [("Data 1", Orange.data.Table, self.set_data1), ("Data 2", Orange.data.Table, self.set_data2)]
        self.outputs = [("Normalized Data 1", Orange.data.Table), ("Normalized Data 2", Orange.data.Table)]
        
        self.selected_method = "Median"
        
        #####
        # GUI
        #####
        
        OWGUI.comboBox(self.controlArea, self, "selected_method",
                       box="Method",
                       tooltip="Normalization type",
                       items=self.METHODS,
                       sendSelectedValue=True,
                       callback=self.on_method_changed)
        
        self.data1 = None
        self.data2 = None
        
        self.resize(200, 50)
        
    def set_data1(self, data=None):
        self.data1 = data
    
    def set_data2(self, data=None):
        self.data2 = data
        
    def handleNewSignals(self):
        self.run_normalization()
    
    def on_method_changed(self):
        self.run_normalization()
        
    def run_normalization(self):
        if self.data1 is not None:
            norm1, norm2 = obiDifscale.normalize(self.data1, self.data2, type=self.selected_method.lower())
            self.send("Normalized Data 1", norm1)
            self.send("Normalized Data 2", norm2)
        else:
            self.send("Normalized Data 1", None)
            self.send("Normalized Data 2", None)
            
if __name__ == "__main__":
    app = QApplication([])
    w = OWNormalization()
    data = Orange.data.Table("GDS")
    w.show()
    app.exec_()
