"""<name>MA Plot</name>
<description>Normalize expression array data on a MA - plot</description>
<icon>icons/Normalize.png</icons>
"""

__widgetdoc__ = """\
.. name: MA Plot
.. description: Normalize expression array data on a MA - plot
.. icon: icons/Normalize.png


"""
from OWWidget import *
from OWGraph import *
import OWGUI
import numpy

import obiExpression
        
import OWConcurrent
        
class OWMAPlot(OWWidget):
    settingsList = ["appendZScore", "appendRIValues"]
    contextHandlers = {"": DomainContextHandler("", ["selectedGroup", "selectedCenterMethod",
                                                     "selectedMergeMethod", "zCutoff"])}
    
    CENTER_METHODS = [("Average", obiExpression.MA_center_average),
                      ("Lowess (fast - interpolated)", obiExpression.MA_center_lowess_fast),
                      ("Lowess", obiExpression.MA_center_lowess)]
    
    MERGE_METHODS = [("Average", numpy.ma.average),
                     ("Median", numpy.ma.median),
                     ("Geometric mean", obiExpression.geometric_mean)]
    
    def __init__(self, parent=None, signalManager=None, name="Normalize Expression Array"):
        OWWidget.__init__(self, parent, signalManager, name, wantGraph=True)
        
        self.inputs = [("Expression array", ExampleTable, self.setData)]
        self.outputs = [("Normalized expression array", ExampleTable, Default), ("Filtered expression array", ExampleTable)]
        
        self.selectedGroup = 0
        self.selectedCenterMethod = 0
        self.selectedMergeMethod = 0
        self.zCutoff = 1.96
        self.appendZScore = False
        self.appendRIValues = False
        self.autoCommit = False
        
        self.loadSettings()
        ## GUI
        self.infoBox = OWGUI.widgetLabel(OWGUI.widgetBox(self.controlArea, "Info", addSpace=True),
                                         "No data on input.")
        
        box = OWGUI.widgetBox(self.controlArea, "Split by", addSpace=True)
        self.groupCombo = OWGUI.comboBox(box, self, "selectedGroup", 
                                         callback=self.onGroupSelection
                                         )
        
        self.centerCombo = OWGUI.comboBox(self.controlArea, self, "selectedCenterMethod",
                                          box="Center Fold-change Using",
                                          items=[name for name, _ in self.CENTER_METHODS],
                                          callback=self.onCenterMethodChange,
                                          addSpace=True
                                          )
        
        self.mergeCombo = OWGUI.comboBox(self.controlArea, self, "selectedMergeMethod",
                                         box="Merge Replicates",
                                         items=[name for name, _ in self.MERGE_METHODS],
                                         tooltip="Select the method for replicate merging",
                                         callback=self.onMergeMethodChange,
                                         addSpace=True
                                         )
        
        box = OWGUI.doubleSpin(self.controlArea, self, "zCutoff", 0.0, 3.0, 0.01,
                               box="Z-Score Cutoff",
                               callback=[self.replotMA, self.commitIf])
        
        OWGUI.separator(self.controlArea)
        
        box = OWGUI.widgetBox(self.controlArea, "Ouput")
        OWGUI.checkBox(box, self, "appendZScore", "Append Z-Scores",
                       tooltip="Append calculated Z-Scores to output",
                       callback=self.commitIf
                       )
        
        OWGUI.checkBox(box, self, "appendRIValues", "Append Log Ratio and Intensity values",
                       tooltip="Append calculated Log Ratio and Intensity values to output data",
                       callback=self.commitIf
                       )
        
        cb = OWGUI.checkBox(box, self, "autoCommit", "Commit on change",
                       tooltip="Commit data on any change",
                       callback=self.commitIf
                       )
        
        b = OWGUI.button(box, self, "Commit", callback=self.commit)
        OWGUI.setStopper(self, b, cb, "changedFlag", callback=self.commit)
        
        self.connect(self.graphButton, SIGNAL("clicked()"), self.saveGraph)
        
        OWGUI.rubber(self.controlArea)
        self.graph = OWGraph(self.mainArea)
        self.graph.setAxisTitle(QwtPlot.xBottom, "Intensity: log<sub>10</sub>(R*G)")
        self.graph.setAxisTitle(QwtPlot.yLeft, "Log ratio: log<sub>2</sub>(R/G)")
        self.graph.showFilledSymbols = True
        self.mainArea.layout().addWidget(self.graph)
        self.groups = []
        self.split_data = None, None
        self.merged_splits = None, None
        self.centered = None, None
        self.changedFlag = False
        self.data = None
        
        self.resize(800, 600)
        
#        self.myThread = WorkerThread()
#        self.myThread.start()
    
        
    def onFinished(self, status):
        self.setEnabled(True)
    
    
    def onUnhandledException(self, ex_info):
        print >> sys.stderr, "Unhandled exception in non GUI thread"
        
        ex_type, ex_val, tb = ex_info
        if ex_type == numpy.linalg.LinAlgError and False:
            self.error(0, "Linear algebra error: %s" % repr(ex_val))
        else:
            sys.excepthook(*ex_info)
        self.setEnabled(True)
    
    
    def onGroupSelection(self):
        if self.data:
            self.updateInfoBox()
            self.splitData()
            self.runNormalization()
        
        
    def onCenterMethodChange(self):
        if self.data:
            self.runNormalization()
        
        
    def onMergeMethodChange(self):
        if self.data:
            self.splitData()
            self.runNormalization()
        
        
    def proposeGroups(self, data):
        col_labels = [attr.attributes.items() for attr in data.domain.attributes]
        col_labels = sorted(reduce(set.union, col_labels, set()))
        col_labels = [(key, value, 1) for key, value in col_labels]
        
        attrs = [attr for attr in data.domain.variables + data.domain.getmetas().values() \
                 if attr.varType == orange.VarTypes.Discrete]
        
        row_labels = [(attr.name, value, 0) for attr in attrs for value in attr.values]
        
        return col_labels + row_labels
    
    
    def setData(self, data):
        self.closeContext("")
        self.data = data
        self.error(0)
        if data is not None:
            self.infoBox.setText("%i genes on input" % len(data))
            self.groups = self.proposeGroups(data)
            self.groupCombo.clear()
            self.groupCombo.addItems(["%s: %s" % (key, value) for key, value, axis in self.groups])
            
            self.openContext("", data)
            self.selectedGroup = min(self.selectedGroup, len(self.groups) - 1)
            
            self.updateInfoBox()
            self.splitData()
            self.runNormalization()
        else:
            self.clear()
        
        
    def clear(self):
        self.groups = []
        self.data = None
        self.centered = None, None
        self.split_data = None, None
        self.merged_splits = None, None
        self.infoBox.setText("No data on input")
        self.send("Normalized expression array", None)
        self.send("Filtered expression array", None)
        
        
    def updateInfoBox(self):
        genes = self.getGeneNames()
        self.infoBox.setText("%i genes on input" % len(self.data))
        
        
    def getSelectedGroup(self):
        return self.groups[self.selectedGroup]
    
    
    def getSelectedGroupSplit(self):
        key, value, axis = self.getSelectedGroup()
        other_values = [v for k, v, a in self.groups if k == key and a == axis and v != value]
        return [(key, value), (key, other_values)], axis
    
    
    def getGeneNames(self):
        key, value, axis = self.getSelectedGroup()
        if axis == 0:
            genes = [str(ex[key]) for ex in self.data]
        else:
            genes = [attr.name for attr in self.data.domain.attributes]
            
        return genes
    
    
    def splitData(self): 
        groups, axis = self.getSelectedGroupSplit()
        self.split_ind = [obiExpression.select_indices(self.data, key, value, axis) for key, value in groups]
        self.split_data = obiExpression.split_data(self.data, groups, axis)
        
        
    def getMerged(self):
        split1, split2 = self.split_data
        (array1, _, _), (array2, _, _) = split1.toNumpyMA(), split2.toNumpyMA()
        
        _, _, axis = self.getSelectedGroup()
        merge_function = self.MERGE_METHODS[self.selectedMergeMethod][1]
        
        merged1 = obiExpression.merge_replicates(array1, axis, merge_function=merge_function)
        merged2 = obiExpression.merge_replicates(array2, axis, merge_function=merge_function)
        self.merged_splits = merged1, merged2
        
        return self.merged_splits
        
        
    def runNormalization(self):
        self.progressBarInit()
        self.progressBarSet(0.0)
        G, R = self.getMerged()
        self.progressBarSet(5.0)
        
        center_method = self.CENTER_METHODS[self.selectedCenterMethod][1]
        
        # TODO: progess bar , lowess can take a long time
        if self.selectedCenterMethod in [1, 2]: #Lowess
            Gc, Rc = center_method(G, R, f = 1./min(500., len(G)/100), iter=1)
        else:
            Gc, Rc = center_method(G, R)
        self.progressBarSet(70.0)
        self.centered = Gc, Rc
        self.z_scores = obiExpression.MA_zscore(Gc, Rc, 1./3.)
        self.progressBarSet(100.0)
        self.plotMA(Gc, Rc, self.z_scores, self.zCutoff)
        self.progressBarFinished()
        
        
    def runNormalizationAsync(self):
        """ Run MA centering and z_score estimation in a separate thread 
        """
        self.setEnabled(False)
        self.error(0)
        self.progressBarInit()
        self.progressBarSet(0.0)
        G, R = self.getMerged()
        self.progressBarSet(5.0)
        
        center_method = self.CENTER_METHODS[self.selectedCenterMethod][1] 
        
        def onCenterResult((Gc, Rc)):
            """ Handle results of MA_center* method
            """
            self.centered = Gc, Rc
            self.progressBarSet(70.0)
            def onZScores(z_scores):
                """ Handle results of MA_z_scores method
                """
                self.z_scores = z_scores
                self.progressBarFinished()
                self.setEnabled(True)
                QTimer.singleShot(50, lambda: self.plotMA(Gc, Rc, self.z_scores, self.zCutoff))
                
            self.z_scores_async = OWConcurrent.createTask(obiExpression.MA_zscore, (Gc, Rc, 1./3.),
                                                          onResult=onZScores,
                                                          onError=self.onUnhandledException)
            
        if self.selectedCenterMethod in [1, 2]: #Lowess
            async_center = OWConcurrent.createTask(center_method, (G, R), {"f": 2./3.,
                                                                           "iter": 1},
                                                   onResult=onCenterResult,
                                                   onError=self.onUnhandledException)
        else:
            async_center = OWConcurrent.createTask(center_method, (G, R),
                                                   onResult=onCenterResult,
                                                   onError=self.onUnhandledException)

        self.async_center = async_center
            
    ## comment out this line if threading creates any problems 
    runNormalization = runNormalizationAsync
    
    
    def plotMA(self, G, R, z_scores, z_cuttof):
        ratio, intensity = obiExpression.ratio_intensity(G, R)
        
        filter = numpy.isfinite(ratio) & numpy.isfinite(intensity) & numpy.isfinite(z_scores)
        for array in [ratio, intensity, z_scores]:
            if numpy.ma.is_masked(array):
                filter &= -array.mask
        
        filtered_ind = numpy.where(filter)
        ratio = numpy.take(ratio, filtered_ind)
        intensity = numpy.take(intensity, filtered_ind)
        z_scores = numpy.take(z_scores, filtered_ind)
        
        red_ind = numpy.where(numpy.ma.abs(z_scores) >= z_cuttof)
        blue_ind = numpy.where(numpy.ma.abs(z_scores) < z_cuttof)
        
        red_xdata, red_ydata = intensity[red_ind], ratio[red_ind]
        blue_xdata, blue_ydata = intensity[blue_ind], ratio[blue_ind]
        self.graph.removeDrawingCurves()
        
        c = self.graph.addCurve("", Qt.black, Qt.black, xData=[0.0, 1.0], yData=[0.0, 0.0], style=QwtPlotCurve.Lines, symbol=QwtSymbol.NoSymbol)
        c.setAxis(QwtPlot.xTop, QwtPlot.yLeft)
        
        self.graph.addCurve("Z >= %.2f" % z_cuttof, QColor(Qt.red), Qt.red, brushAlpha=100, size=6, enableLegend=True, xData=list(red_xdata), yData=list(red_ydata), autoScale=True)
        self.graph.addCurve("Z < %.2f" % z_cuttof, QColor(Qt.blue), Qt.blue, brushAlpha=100, size=6, enableLegend=True, xData=list(blue_xdata), yData=list(blue_ydata), autoScale=True)
        
        self.graph.setAxisScale(QwtPlot.xTop, 0.0, 1.0)
         
        self.graph.replot()
        
        
    def replotMA(self):
        if self.data and self.centered:
            Gc, Rc = self.centered
            self.plotMA(Gc, Rc, self.z_scores, self.zCutoff)
        
        
    def commitIf(self):
        if self.autoCommit and self.changedFlag:
            self.commit()
        else:
            self.changedFlag = True
            
            
    def commit(self):
        if not self.data:
            return
        
        G, R = self.merged_splits
        Gc, Rc = self.centered
        ind1, ind2 = self.split_ind
        
        gfactor = Gc / G
        
        domain = orange.Domain(self.data.domain.attributes, self.data.domain.classVar)
        domain.addmetas(self.data.domain.getmetas())
        
        _, _, axis = self.getSelectedGroup()
        if self.appendZScore and axis == 1:
            attr = orange.FloatVariable("Z-Score")
            if not hasattr(self, "z_score_mid"):
                self.z_score_mid = orange.newmetaid()
            mid = self.z_score_mid
            domain.addmeta(mid, attr)
            
        if self.appendRIValues and axis == 1:
            r_attr = orange.FloatVariable("Log Ratio")
            i_attr = orange.FloatVariable("Intensity")
            if not hasattr(self, "_r_mid"):
                self._r_mid = orange.newmetaid()
                self._i_mid = orange.newmetaid()
            r_mid, i_mid = self._r_mid, self._i_mid
            domain.addmeta(r_mid, r_attr)
            domain.addmeta(i_mid, i_attr)
            
        data = orange.ExampleTable(domain, self.data)
        
        def finite_nonmasked(g):
            return numpy.isfinite(g) and g is not numpy.ma.masked
        
        if axis == 0:
            for i, gf in zip(ind1, gfactor):
                for attr in domain.attributes:
                    if not data[i][attr].isSpecial() and finite_nonmasked(gf):
                        data[i][attr] = data[i][attr] * gf
        else:
            ratio, intensity = obiExpression.ratio_intensity(*self.centered) #*self.merged_splits)
            for ex, gf, r, inten, z in zip(data, gfactor, ratio, intensity, self.z_scores):
                for i in ind1:
                    if not ex[i].isSpecial() and finite_nonmasked(gf):
                        ex[i] = float(ex[i]) * gf
                if self.appendZScore:
                    ex[attr] = z if finite_nonmasked(z) else "?"
                    
                if self.appendRIValues:
                    ex[r_attr] = r if finite_nonmasked(r) else "?"
                    ex[i_attr] = inten if finite_nonmasked(inten) else "?"
        
        filtered_ind = list(numpy.ma.abs(self.z_scores.filled(0)) >= self.zCutoff)
        if axis == 0:
            attrs = [attr for attr, filtered in zip(domain.attributes, filtered_ind) if filtered]
            filtered_domain = orange.Domain(attrs, domain.classVar)
            filtered_domain.addmetas(domain.getmetas())
            filtered_data = orange.ExampleTable(filtered_domain, data)
        else:
            filtered_data = data.select([int(b) for b in filtered_ind])
            
        self.send("Normalized expression array", data)
        self.send("Filtered expression array", filtered_data)
        
        
    def saveGraph(self):
        from OWDlgs import OWChooseImageSizeDlg
        dlg = OWChooseImageSizeDlg(self.graph, parent=self)
        dlg.exec_()
        
        
if __name__ == "__main__":
    app = QApplication(sys.argv)
    w= OWMAPlot()
#    data = orange.ExampleTable(os.path.expanduser("~/GDS1210.tab"))
    data = orange.ExampleTable(os.path.expanduser("../../../doc/datasets/brown-selected.tab"))
    w.setData(data)
    w.show()
    app.exec_()
        
        