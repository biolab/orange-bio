"""
<name>Gene Selection</name>
<description>Gene scoring and selection.</description>
<priority>230</priority>
<icon>icons/FeatureSelection.png</icon>
"""
import orange

from obiExpression import *

from OWGraph import *
from OWGraphTools import PolygonCurve
from OWWidget import *
from OWHist import OWInteractiveHist
from OWToolbars import ZoomSelectToolbar
from obiGEO import transpose
from collections import defaultdict

import OWGUI

class ScoreHist(OWInteractiveHist):
    def __init__(self, master, parent=None, type="hiTail"):
        OWInteractiveHist.__init__(self, parent, type=type)
        self.master = master
        self.setAxisTitle(QwtPlot.xBottom, "Score")
        self.setAxisTitle(QwtPlot.yLeft, "Frequency")
        
    def setBoundary(self, low, hi):
        OWInteractiveHist.setBoundary(self, low, hi)
        self.master.UpdateSelectedInfoLabel(low, hi)
        if self.master.autoCommit:
            self.master.Commit()
        
class OWFeatureSelection(OWWidget):
    settingsList = ["methodIndex", "dataLabelIndex", "computeNullDistribution", "permutationsCount", "selectPValue", "autoCommit"]
##    contextHandlers={"":DomainContextHandler("",[])}
    def __init__(self, parent=None, signalManager=None, name="Gene selection"):
        OWWidget.__init__(self, parent, signalManager, name, wantGraph=True, showSaveGraph=True)
        self.inputs = [("Examples", ExampleTable, self.SetData)]
        self.outputs = [("Examples with selected attributes", ExampleTable), ("Examples with remaining attributes", ExampleTable), ("Selected attributes", ExampleTable)]

        self.methodIndex = 0
        self.dataLabelIndex = 0
        self.computeNullDistribution = False
        self.permutationsCount = 10
        self.autoCommit = False
        self.selectNBest = 20
        self.selectPValue = 0.01

        self.oneTailTestHi = oneTailTestHi = lambda attr, low, hi:self.scores.get(attr,0)>=hi
        self.oneTailTestLow = oneTailTestLow = lambda attr, low, hi:self.scores.get(attr,0)<=low
        self.twoTailTest = twoTailTest = lambda attr, low, hi:self.scores.get(attr,0)>=hi or self.scores.get(attr,0)<=low
        self.middleTest = middleTest = lambda attr, low, hi:self.scores.get(attr,0)<=hi and self.scores.get(attr,0)>=low
        self.histType = {oneTailTestHi:"hiTail", oneTailTestLow:"lowTail", twoTailTest:"twoTail", middleTest:"middle"}
##      
        self.scoreMethods = [("fold change", lambda attr, data: MA_fold_change()(attr, data), twoTailTest),
                             ("log2 fold change", lambda attr, data: math.log(max(min(MA_fold_change()(attr, data), 1e300), 1e-300), 2.0), twoTailTest),
                             ("t-test",lambda attr, data: MA_t_test()(attr, data), twoTailTest),
                             ("t-test p-value",lambda attr, data: MA_t_test(prob=True)(attr, data), oneTailTestLow),
                             ("anova", lambda attr, data: MA_anova()(attr, data), oneTailTestHi),
                             ("anova p-value", lambda attr, data: MA_anova(prob=True)(attr, data), oneTailTestLow),
                             ("signal to noise ratio", lambda attr, data: MA_signalToNoise()(attr, data), twoTailTest),
                             ("info gain", orange.MeasureAttribute_info, oneTailTestHi),
                             ("chi-square", orange.MeasureAttribute_chiSquare, oneTailTestHi)]

        boxHistogram = OWGUI.widgetBox(self.mainArea)
        self.histogram = ScoreHist(self, boxHistogram)
        boxHistogram.layout().addWidget(self.histogram)
        self.histogram.show()
        
        box = OWGUI.widgetBox(self.controlArea, "Info", addSpace=True)
        self.dataInfoLabel = OWGUI.widgetLabel(box, "\n\n")
        self.selectedInfoLabel = OWGUI.widgetLabel(box, "")

        #self.testRadioBox = OWGUI.radioButtonsInBox(self.controlArea, self, "methodIndex", [sm[0] for sm in self.scoreMethods], box="Scoring Method", callback=self.Update, addSpace=True)
        # OWGUI.comboBoxWithCaption(box, self, "ptype", items=nth(self.permutationTypes, 0), \
        #             tooltip="Permutation type.", label="Permutate")
        box1 = OWGUI.widgetBox(self.controlArea, "Scoring Method")
        self.testRadioBox = OWGUI.comboBox(box1, self, "methodIndex", items=[sm[0] for sm in self.scoreMethods], callback=self.Update)
        self.dataLabelComboBox = OWGUI.comboBox(box1, self, "dataLabelIndex", "Attribute labels", callback=self.Update, tooltip="Use attribute labels for score computation")
##        self.useClassCheck = OWGUI.checkBox(box1, self, "useClass", "Use class information", callback=self.Update, tooltip="Use class information for score computation", disables=[self.da)
    
        ZoomSelectToolbar(self, self.controlArea, self.histogram, buttons=[ZoomSelectToolbar.IconSelect, ZoomSelectToolbar.IconZoom, ZoomSelectToolbar.IconPan])
        OWGUI.separator(self.controlArea)
        
        box = OWGUI.widgetBox(self.controlArea, "Selection", addSpace=True)
        callback = lambda: self.histogram.setBoundary(self.histogram.lowerBoundary, self.histogram.upperBoundary) if self.data else None
        self.upperBoundarySpin = OWGUI.doubleSpin(box, self, "histogram.upperBoundary", min=-1e6, max=1e6, step= 1e-6, label="Upper threshold:", callback=callback, callbackOnReturn=True)
        self.lowerBoundarySpin = OWGUI.doubleSpin(box, self, "histogram.lowerBoundary", min=-1e6, max=1e6, step= 1e-6, label="Lower threshold:", callback=callback, callbackOnReturn=True)
        check = OWGUI.checkBox(box, self, "computeNullDistribution", "Compute null distribution", callback=self.Update)

        check.disables.append(OWGUI.spin(box, self, "permutationsCount", min=1, max=10, label="Permutations:", callback=self.Update, callbackOnReturn=True))

        box1 = OWGUI.widgetBox(box, orientation='horizontal')
        check.disables.append(OWGUI.doubleSpin(box1, self, "selectPValue" , min=2e-7, max=1.0, step=1e-7, label="P-value:"))
        check.disables.append(OWGUI.button(box1, self, "Select", callback=self.SelectPBest))
        check.makeConsistent()

        box1 = OWGUI.widgetBox(box, orientation='horizontal')
        OWGUI.spin(box1, self, "selectNBest", 0, 10000, step=1, label="Best Ranked:")
        OWGUI.button(box1, self, "Select", callback=self.SelectNBest)

        box = OWGUI.widgetBox(self.controlArea, "Commit")
        OWGUI.button(box, self, "&Commit", callback=self.Commit)
        OWGUI.checkBox(box, self, "autoCommit", "Commit on change")
        OWGUI.rubber(self.controlArea)

        self.connect(self.graphButton, SIGNAL("clicked()"), self.histogram.saveToFile)
        
        self.loadSettings()

        self.data = None
        self.discData = None
        self.scoreCache = {}
        self.nullDistCache = {}
        self.cuts = {}
        self.discretizer = orange.EquiNDiscretization(numberOfIntervals=5)
        self.transposedData = False
        self.nullDistribution = []
        self.scores = {}

        self.resize(700, 600)
        
    def SetData(self, data):
        self.error(0)
        self.warning(0)
        self.scoreCache = {}
        self.nullDistCache = {}
        self.discData = None
        self.data = data
        self.transposedData = None
        disabled = []
        if self.data:
            self.dataLabels = reduce(lambda dict, tags: [dict[key].add(value) for key, value in tags.items()] and False or dict,
                                       [attr.attributes for attr in self.data.domain.attributes],
                                       defaultdict(set))
            self.dataLabels = [key for key, value in self.dataLabels.items() if len(value) > 1]
        else:
            self.dataLabels = []
        self.dataLabelComboBox.clear()
        if self.data and data.domain.classVar:
            self.dataLabels = ["(None)"] + self.dataLabels
        self.dataLabelComboBox.addItems(self.dataLabels)
        
        self.dataLabelComboBox.setDisabled(len(self.dataLabels) == 1)
        self.dataLabelIndex = max(min(self.dataLabelIndex, len(self.dataLabels) - 1), 0)
        
        self.Update()
        self.Commit()

    def ComputeAttributeScore(self, data, scoreFunc, label="(None)", useCache=True, progressCallback=None):
##        if (scoreFunc, label) in self.scoreCache and useCache:
##            return self.scoreCache[scoreFunc, label]
        attributes = data.domain.attributes
##        if self.methodIndex in [2,3,4,5,6] and len(data.domain.classVar.values)>2:
        if self.methodIndex in [0, 1, 2, 3, 6] and len(data.domain.classVar.values) > 2:
            self.warning(0, self.scoreMethods[self.methodIndex][0]+" works only on two class data (using only first two class values for computation)")
        else:
            self.warning(0)
##        self.progressBarInit()
        if scoreFunc==orange.MeasureAttribute_info or scoreFunc==orange.MeasureAttribute_chiSquare:
            if self.discData and useCache:
                data = self.discData
                newAttrs = data.domain.attributes
            else:
                newAttrs = [self.discretizer(attr, data) for attr in attributes]
##                data = data.select(newAttrs + [data.domain.classVar])
                newDomain = orange.Domain(newAttrs + [data.domain.classVar])
                table = orange.ExampleTable(newDomain)
                for i, ex in enumerate(data):
                    table.append(orange.Example(newDomain, ex))
                    self.progressBarSet(100.0*i/len(data))
                if useCache:
                    self.discData = data = table
        else:
            newAttrs = attributes
        scores = {}
        milestones = set(range(0, len(attributes), max(len(attributes)/100, 1)))
        for i, (attr, newAttr) in enumerate(zip(attributes, newAttrs)):
            scores[attr] = scoreFunc(newAttr, data)
            if progressCallback and i in milestones:
                progressCallback(100.0*i/len(attributes))
##        self.progressBarFinished()
        if useCache:
            self.scoreCache[scoreFunc, label] = scores
        return scores

    def ComputeNullDistribution(self, data, scoreFunc, label="(None)", progressCallback=None):
        if (scoreFunc, self.permutationsCount, label) in self.nullDistCache:
            return self.nullDistCache[scoreFunc, self.permutationsCount, label]

        originalClasses = [ex.getclass() for ex in data]
        scores = []
        import random
        random.seed(0)
        for i in range(self.permutationsCount):
            permClasses = list(originalClasses)
            random.shuffle(permClasses)
            for ex, class_ in zip(data, permClasses):
                ex.setclass(class_)
            _progressCallback = lambda val: progressCallback(100.0*i/self.permutationsCount + float(val)/self.permutationsCount) if \
                               progressCallback else None
            scores.extend(self.ComputeAttributeScore(data, scoreFunc, useCache=False, progressCallback=_progressCallback).values())

        for ex, class_ in zip(data, originalClasses):
            ex.setclass(class_)

        self.nullDistCache[scoreFunc, self.permutationsCount, label] = scores
        return scores
        
    def Update(self):
        if not self.data:
            return
        self.error(0)
        if self.dataLabels[self.dataLabelIndex] != "(None)":
            try:
                data = transpose_labels_to_class(self.data, classlabel=self.dataLabels[self.dataLabelIndex])
                self.transposedData = data
##                print data.domain.classVar
            except Exception, ex:
                raise
        else:
            data = self.data
            self.transposedData = None
            
        if not (data and data.domain.classVar):
            self.error(0, "Class var missing!")
            self.histogram.clear()
            return
        
        if (data.domain.classVar.values) == 2:
            disabled = [4, 5]
        elif len(data.domain.classVar.values) > 2:
           disabled = [0, 1, 2, 3, 6]
##        for i, button in enumerate(self.testRadioBox.buttons):
##            button.setDisabled(i in disabled)
            
        if data and data.domain.classVar:
            label = self.dataLabels[self.dataLabelIndex]
            self.progressBarInit()
            self.scores = self.ComputeAttributeScore(data, self.scoreMethods[self.methodIndex][1], label, progressCallback=self.progressBarSet,)
            if self.computeNullDistribution:
                self.nullDistribution = self.ComputeNullDistribution(data, self.scoreMethods[self.methodIndex][1], label, progressCallback=self.progressBarSet)
            self.progressBarFinished()
            self.histogram.type = self.histType[self.scoreMethods[self.methodIndex][2]]
            self.histogram.setValues(self.scores.values())
            self.histogram.setBoundary(*self.cuts.get(self.methodIndex, (0, 0)))
            if self.computeNullDistribution:
##                nullY, nullX = numpy.histogram(self.nullDistribution, bins=100)
                nullY, nullX = numpy.histogram(self.nullDistribution, bins=self.histogram.xData)
                self.histogram.nullCurve = self.histogram.addCurve("nullCurve", Qt.black, Qt.black, 6, symbol = QwtSymbol.NoSymbol, style = QwtPlotCurve.Steps, xData = nullX, yData = nullY/self.permutationsCount)
                
                minx = min(min(nullX), self.histogram.minx)
                maxx = max(max(nullX), self.histogram.maxx)
                miny = min(min(nullY/self.permutationsCount), self.histogram.miny)
                maxy = max(max(nullY/self.permutationsCount), self.histogram.maxy)

                self.histogram.setAxisScale(QwtPlot.xBottom, minx - (0.05 * (maxx - minx)), maxx + (0.05 * (maxx - minx)))
                self.histogram.setAxisScale(QwtPlot.yLeft, miny - (0.05 * (maxy - miny)), maxy + (0.05 * (maxy - miny)))                                            
            state = dict(hiTail=(False, True), lowTail=(True, False), twoTail=(True, True))
            for spin, visible in zip((self.upperBoundarySpin, self.lowerBoundarySpin), state[self.histogram.type]):
                spin.setVisible(visible)
            
##            if self.methodIndex in [2, 3, 5, 6]:
            if self.methodIndex in [0, 2, 4, 6]:
                classValues = data.domain.classVar.values
                if self.methodIndex == 0: ## fold change is centered on 1.0
                    x1, y1 = (self.histogram.minx + 1) / 2 , self.histogram.maxy
                    x2, y2 = (self.histogram.maxx + 1) / 2 , self.histogram.maxy
                else:
                    x1, y1 = (self.histogram.minx) / 2 , self.histogram.maxy
                    x2, y2 = (self.histogram.maxx) / 2 , self.histogram.maxy
                self.histogram.addMarker(classValues[1], x1, y1)
                self.histogram.addMarker(classValues[0], x2, y2)
            self.histogram.replot()
            
        else:
            self.histogram.clear()
            
        self.UpdateDataInfoLabel()
            
    def UpdateDataInfoLabel(self):
        data = self.transposedData if self.transposedData else self.data
        if data:
            text = "%i samples, %i genes\n" % (len(data), len(data.domain.attributes))
            if data.domain.classVar:
                text = text+"Sample labels: %s" % data.domain.classVar.name
            else:
                text = text+"Info with sample lables"
        else:
            text = "No data on input\n"
        self.dataInfoLabel.setText(text)

    def UpdateSelectedInfoLabel(self, cutOffLower=0, cutOffUpper=0):
        self.cuts[self.methodIndex] = (cutOffLower, cutOffUpper)
        data = self.transposedData if self.transposedData else self.data
        if data:
            test = self.scoreMethods[self.methodIndex][2]
            self.selectedInfoLabel.setText("%i selected genes" %len([attr for attr in data.domain.attributes if test(attr, cutOffLower, cutOffUpper)])) #self.scores.get(attr,0)>cutOffUpper or self.scores.get(attr,0)<cutOffLower]))
        else:
            self.selectedInfoLabel.setText("0 selected genes")

    def SelectNBest(self):
        scores = self.scores.items()
        scores.sort(lambda a,b:cmp(a[1], b[1]))
        if not scores:
            return
        if self.scoreMethods[self.methodIndex][2]==self.oneTailTestHi:
            scores = scores[-max(self.selectNBest, 1):]
            self.histogram.setBoundary(scores[0][1], scores[0][1])
        elif self.scoreMethods[self.methodIndex][2]==self.oneTailTestLow:
            scores = scores[:max(self.selectNBest,1)]
            self.histogram.setBoundary(scores[-1][1], scores[-1][1])
        else:
            scoresHi = scores[-max(min(self.selectNBest, len(scores)/2), 1):]
            scoresLo = scores[:max(min(self.selectNBest, len(scores)/2), 1)]
            scores = [(abs(score), 1) for attr, score in scoresHi] + [(abs(score), -1) for attr, score in scoresLo]
            if self.scoreMethods[self.methodIndex][0]=="fold change": ## fold change is on a logaritmic scale
                scores =  [(abs(math.log(max(min(score, 1e300), 1e-300), 2.0)), sign) for score, sign in scores]
            scores.sort()
            scores = scores[-max(self.selectNBest, 1):]
            countHi = len([score for score, sign in scores if sign==1])
            countLo = len([score for score, sign in scores if sign==-1])
            cutHi = scoresHi[-countHi][1] if countHi else scoresHi[-1][1] + 1e-7
            cutLo = scoresLo[countLo-1][1] if countLo else scoresLo[0][1] - 1e-7
            self.histogram.setBoundary(cutLo, cutHi)

    def SelectPBest(self):
        if not self.nullDistribution:
            return
        nullDist = sorted(self.nullDistribution)
        test = self.scoreMethods[self.methodIndex][2]
        count = int(len(nullDist)*self.selectPValue)
        if test == self.oneTailTestHi:
            cut = nullDist[-count] if count else nullDist[-1] + 1e-7
            self.histogram.setBoundary(cut, cut)
##            print cut
        elif test == self.oneTailTestLow:
            cut = nullDist[count - 1] if count else nullDist[0] - 1e-7
            self.histogram.setBoundary(cut, cut)
##            print cut
        elif count:
            scoresHi = nullDist[-count:]
            scoresLo = nullDist[:count]
            scores = [(abs(score), 1) for score in scoresHi] + [(abs(score), -1) for score in scoresLo]
            if self.scoreMethods[self.methodIndex][0] == "fold change": ## fold change is on a logaritmic scale
                scores =  [(abs(math.log(max(min(score, 1e300), 1e-300), 2.0)), sign) for score, sign in scores]
            scores = sorted(scores)[-count:]
            countHi = len([score for score, sign in scores if sign==1])
            countLo = len([score for score, sign in scores if sign==-1])
            cutHi = scoresHi[-countHi] if countHi else scoresHi[-1] + 1e-7
            cutLo = scoresLo[countLo-1] if countLo else scoresLo[0] - 1e-7
            self.histogram.setBoundary(cutLo, cutHi)
        else:
            self.histogram.setBoundary(nullDist[0] - 1e-7, nullDist[-1] + 1e-7)
        
    def Commit(self):
        data =  self.transposedData if self.transposedData else self.data
        if data and data.domain.classVar:
            cutOffUpper = self.histogram.upperBoundary
            cutOffLower = self.histogram.lowerBoundary
            test = self.scoreMethods[self.methodIndex][2]
            
            selectedAttrs = [attr for attr in data.domain.attributes if  test(attr, cutOffLower, cutOffUpper)] #self.scores.get(attr,0)>cutOffUpper or self.scores.get(attr,0)<cutOffLower]
##            data = self.data.select(selectedAttrs+[self.data.domain.classVar])
            newdomain = orange.Domain(selectedAttrs, data.domain.classVar)
            newdomain.addmetas(data.domain.getmetas())
##            print data.domain["_$_$sample"], [d["_$_$sample"] for d in data]
            newdata = orange.ExampleTable(newdomain, data)
            if self.transposedData and selectedAttrs:
                newdata = transpose_class_to_labels(newdata, classlabel=self.dataLabels[self.dataLabelIndex])
            self.send("Examples with selected attributes", newdata if selectedAttrs else None)
            
            remainingAttrs = [attr for attr in data.domain.attributes if  not test(attr, cutOffLower, cutOffUpper)] #self.scores.get(attr,0)>cutOffUpper or self.scores.get(attr,0)<cutOffLower]
##            data = self.data.select(remainingAttrs+[self.data.domain.classVar])
            newdomain = orange.Domain(remainingAttrs, data.domain.classVar)
            newdomain.addmetas(data.domain.getmetas())
            newdata = orange.ExampleTable(newdomain, data)
            if self.transposedData and remainingAttrs:
                newdata = transpose_class_to_labels(newdata, classlabel=self.dataLabels[self.dataLabelIndex])
            self.send("Examples with remaining attributes", newdata if remainingAttrs else None)
            
            domain = orange.Domain([orange.StringVariable("label"), orange.FloatVariable(self.scoreMethods[self.methodIndex][0])], False)
            self.send("Selected attributes", orange.ExampleTable([orange.Example(domain, [attr.name, self.scores.get(attr, 0)]) for attr in selectedAttrs]) if selectedAttrs else None)
            
        else:
            self.send("Examples with selected attributes", None)
            self.send("Examples with remaining attributes", None)
            self.send("Selected attributes", None)


def _float_or_na(x):
    if x.isSpecial():
        return "?"
    return float(x)

def transpose_class_to_labels(data, attcol="_$_$sample", classlabel="group"):
    """Converts data with genes as attributes to data with genes in rows."""
    if attcol in [v.name for v in data.domain.getmetas().values()]:
        atts = [orange.FloatVariable(str(d[attcol])) for d in data]
    else:
        atts = [orange.FloatVariable("S%d" % i) for i in range(len(data))]
    for i, d in enumerate(data):
        if classlabel:
            atts[i].attributes[classlabel] = str(d.getclass())
        for meta in data.domain.getmetas().values():
            if meta.name != attcol:
                atts[i].attributes[meta.name] = str(d[meta])
    domain = orange.Domain(atts, False)
    
    newdata = []
    for a in data.domain.attributes:
        newdata.append([_float_or_na(d[a]) for d in data])

##    gene = orange.StringVariable("gene")
##    id = orange.newmetaid()
    new = orange.ExampleTable(domain, newdata)
##    new.domain.addmeta(id, gene)
    for label in reduce(set.union, [attr.attributes.keys() for attr in data.domain.attributes], set()):
        new.domain.addmeta(orange.newmetaid(), orange.StringVariable(label))
    for i, d in enumerate(new):
##        d[gene] = data.domain.attributes[i].name
        for label, value in data.domain.attributes[i].attributes.items():
            d[label] = value

    return new

def transpose_labels_to_class(data, classlabel="group"):
    """Converts data with genes in rows to data with genes as attributes."""
    if "gene" in [v.name for v in data.domain.getmetas().values()]:
        atts = [orange.FloatVariable(str(d["gene"])) for d in data]
    else:
        atts = [orange.FloatVariable("A%d" % i) for i in range(len(data))]

    metas = data.domain.getmetas().values()
    for a, d in zip(atts, data):
        for meta in metas:
            if not d[meta].isSpecial():
                a.attributes[meta.name] = str(d[meta])
    classvalues = list(set([a.attributes.get(classlabel) for a in data.domain.attributes]) - set([None]))
    classvar = orange.EnumVariable(classlabel, values=classvalues)
    domain = orange.Domain(atts + [classvar])

    labels = reduce(set.union, [a.attributes for a in data.domain.attributes], set()) - set([classlabel])
    domain.addmetas(dict([(orange.newmetaid(), orange.StringVariable(label)) for label in labels if label]))
    
    newdata = []
    attrs = []
    
    for a in data.domain.attributes:
        if classlabel in a.attributes:
            newdata.append([_float_or_na(d[a]) for d in data] + [a.attributes[classlabel]])
            attrs.append(a)
            
    sample = orange.StringVariable("_$_$sample")
    id = orange.newmetaid()
    new = orange.ExampleTable(domain, newdata)
    new.domain.addmeta(id, sample)
    for i, d in enumerate(new):
        d[sample] = attrs[i].name
##        print d[sample], attrs[i].name
        for label in labels:
            if label in attrs[i].attributes:
                d[label] = attrs[i].attributes[label]

    return new



if __name__=="__main__":
    import sys
    app = QApplication(sys.argv)
    data = orange.ExampleTable("E:\\out1.tab")
    w = OWFeatureSelection()
    w.show()
    w.SetData(data)
    app.exec_()
    w.saveSettings()
