"""
<name>Differential expression</name>
<description>Gene differential expression scoring and selection.</description>
<priority>1010</priority>
<icon>icons/GeneSelection.svg</icon>
"""

from __future__ import absolute_import, with_statement

from collections import defaultdict
from functools import wraps
from operator import add

import numpy as np
import numpy.ma as ma

import orange
from Orange.OrangeWidgets import OWGUI
from Orange.OrangeWidgets.OWGraph import *
from Orange.OrangeWidgets.OWGraphTools import PolygonCurve
from Orange.OrangeWidgets.OWHist import OWInteractiveHist
from Orange.OrangeWidgets.OWToolbars import ZoomSelectToolbar
from Orange.OrangeWidgets.OWWidget import *

from ..obiExpression import *


NAME = "Differential Expression"
DESCRIPTION = "Gene selection by differential expression analysis."
ICON = "icons/GeneSelection.svg"
PRIORITY = 1010

INPUTS = [("Examples", Orange.data.Table, "set_data")]
OUTPUTS = [("Example table with selected genes", Orange.data.Table, Default),
           ("Example table with remaining genes", Orange.data.Table),
           ("Selected genes", Orange.data.Table)]

REPLACES = ["_bioinformatics.widgets.OWFeatureSelection.OWFeatureSelection"]


class ExpressionSignificance_TTest_PValue(ExpressionSignificance_TTest):
    def __call__(self, *args, **kwargs):
        return [(key, pval) for key, (t, pval) in \
                ExpressionSignificance_TTest.__call__(self, *args, **kwargs)]


class ExpressionSignificance_TTest_T(ExpressionSignificance_TTest):
    def __call__(self, *args, **kwargs):
        return [(key, t) for key, (t, pval) in \
                ExpressionSignificance_TTest.__call__(self, *args, **kwargs)]


class ExpressionSignificance_ANOVA_PValue(ExpressionSignificance_ANOVA):
    def __call__(self, *args, **kwargs):
        return [(key, pval) for key, (t, pval) in \
                ExpressionSignificance_ANOVA.__call__(self, *args, **kwargs)]


class ExpressionSignificance_ANOVA_F(ExpressionSignificance_ANOVA):
    def __call__(self, *args, **kwargs):
        return [(key, f) for key, (f, pval) in \
                ExpressionSignificance_ANOVA.__call__(self, *args, **kwargs)]


class ExpressionSignificance_Log2FoldChange(ExpressionSignificance_FoldChange):
    def __call__(self, *args, **kwargs):
        return [(key, math.log(fold, 2.0) if fold > 1e-300 and fold < 1e300 else 0.0) \
                for key, fold in ExpressionSignificance_FoldChange.__call__(self, *args, **kwargs)]


class ExpressionSignigicance_MannWhitneyu_U(ExpressionSignificance_MannWhitneyu):
    def __call__(self, *args, **kwargs):
        return [(key, u) for key, (u, p_val) in \
                ExpressionSignificance_MannWhitneyu.__call__(self, *args, **kwargs)]


class ScoreHist(OWInteractiveHist):
    def __init__(self, master, parent=None, type="hiTail"):
        OWInteractiveHist.__init__(self, parent, type=type)
        self.master = master
        self.setAxisTitle(QwtPlot.xBottom, "Score")
        self.setAxisTitle(QwtPlot.yLeft, "Frequency")
        self.activateSelection()

    def setBoundary(self, low, hi):
        OWInteractiveHist.setBoundary(self, low, hi)
        self.master.on_boundary_change(low, hi)


def disable_controls(method):
    """ Disable the widget's control area during the duration of this call.
    """
    @wraps(method)
    def f(self, *args, **kwargs):
        self.controlArea.setDisabled(True)
        qApp.processEvents()
        try:
            return method(self, *args, **kwargs)
        finally:
            self.controlArea.setDisabled(False)
    return f

from Orange.orng.orngDataCaching import data_hints

from .OWGenotypeDistances import SetContextHandler
from .OWVulcanoPlot import LabelSelectionWidget


class OWFeatureSelection(OWWidget):
    settingsList = [
        "method_index", "dataLabelIndex", "compute_null",
        "permutations_count", "selectPValue", "auto_commit",
        "thresholds"]

    contextHandlers = {
        "Data": DomainContextHandler("Data", ["genes_in_columns"]),
        "TargetSelection":
            SetContextHandler("TargetSelection", findImperfect=False)
    }

    def __init__(self, parent=None, signalManager=None, name="Gene selection"):
        OWWidget.__init__(self, parent, signalManager, name, wantGraph=True, showSaveGraph=True)
        self.inputs = [("Examples", ExampleTable, self.set_data)]
        self.outputs = [("Example table with selected genes", ExampleTable), ("Example table with remaining genes", ExampleTable), ("Selected genes", ExampleTable)]

        self.method_index = 0
        self.genes_in_columns = False
        self.compute_null = False
        self.permutations_count = 10
        self.auto_commit = False
        self.selectNBest = 20
        self.selectPValue = 0.01
        self.data_changed_flag = False
        self.add_scores_to_output = True
        self.thresholds = {
            "fold change": (0.5, 2.),
            "log2 fold change": (-1, 1),
            "t-test": (-2, 2),
            "t-test p-value": (0.01, 0.01),
        }

        self.oneTailTestHi = oneTailTestHi = lambda array, low, hi: array >= hi
        self.oneTailTestLow = oneTailTestLow = lambda array, low, hi: array <= low
        self.twoTailTest = twoTailTest = lambda array, low, hi: (array >= hi) | (array <= low)
        self.middleTest = middleTest = lambda array, low, hi: (array <= hi) | (array >= low)
        
        self.histType = {oneTailTestHi:"hiTail", oneTailTestLow:"lowTail", twoTailTest:"twoTail", middleTest:"middle"}

        # [(name, func, tail test, two sample test), ...]
        self.score_methods = [("fold change", ExpressionSignificance_FoldChange, twoTailTest, True),
                             ("log2 fold change", ExpressionSignificance_Log2FoldChange, twoTailTest, True),
                             ("t-test", ExpressionSignificance_TTest_T, twoTailTest, True),
                             ("t-test p-value", ExpressionSignificance_TTest_PValue, oneTailTestLow, True),
                             ("anova", ExpressionSignificance_ANOVA_F, oneTailTestHi, False),
                             ("anova p-value", ExpressionSignificance_ANOVA_PValue, oneTailTestLow, False),
                             ("signal to noise ratio", ExpressionSignificance_SignalToNoise, twoTailTest, True),
                             ("info gain", ExpressionSignificance_Info, oneTailTestHi, True),
                             ("chi-square", ExpressionSignificance_ChiSquare, oneTailTestHi, True),
                             ("mann-whitney", ExpressionSignigicance_MannWhitneyu_U, oneTailTestLow, True)]

        boxHistogram = OWGUI.widgetBox(self.mainArea)
        self.histogram = ScoreHist(self, boxHistogram)
        boxHistogram.layout().addWidget(self.histogram)
        self.histogram.show()
        
        box = OWGUI.widgetBox(self.controlArea, "Info")
        self.dataInfoLabel = OWGUI.widgetLabel(box, "\n\n")
        self.dataInfoLabel.setWordWrap(True)
        self.selectedInfoLabel = OWGUI.widgetLabel(box, "")

        box1 = OWGUI.widgetBox(self.controlArea, "Scoring Method")
        self.testRadioBox = OWGUI.comboBox(box1, self, "method_index",
                                    items=[sm[0] for sm in self.score_methods],
                                    callback=[self.on_scoring_method_changed, self.update_scores])
        
        box = OWGUI.widgetBox(self.controlArea, "Target Labels")
        self.label_selection_widget = LabelSelectionWidget(self)
        self.label_selection_widget.setMaximumHeight(150)
        box.layout().addWidget(self.label_selection_widget)
        self.connect(self.label_selection_widget,
                     SIGNAL("selection_changed()"),
                     self.on_target_changed)
        
        self.connect(self.label_selection_widget,
                     SIGNAL("label_activated(int)"),
                     self.on_label_activated)
        
        self.genes_in_columns_check = OWGUI.checkBox(box, self, "genes_in_columns",
                                                  "Genes in columns",
                                                  callback=self.on_genes_in_columns_change)

        box = OWGUI.widgetBox(self.controlArea, "Selection")
        box.layout().setSpacing(0)

        self.upperBoundarySpin = OWGUI.doubleSpin(box, self, "histogram.upperBoundary",
                                                  min=-1e6, max=1e6, step= 1e-6,
                                                  label="Upper threshold:", 
                                                  callback=self.update_boundary, 
                                                  callbackOnReturn=True)
        
        self.lowerBoundarySpin = OWGUI.doubleSpin(box, self, "histogram.lowerBoundary", 
                                                  min=-1e6, max=1e6, step= 1e-6, 
                                                  label="Lower threshold:", 
                                                  callback=self.update_boundary, 
                                                  callbackOnReturn=True)
        
        check = OWGUI.checkBox(box, self, "compute_null", "Compute null distribution",
                               callback=self.update_scores)

        check.disables.append(OWGUI.spin(box, self, "permutations_count", min=1, max=10, 
                                         label="Permutations:", callback=self.update_scores, 
                                         callbackOnReturn=True))

        box1 = OWGUI.widgetBox(box, orientation='horizontal')
        check.disables.append(OWGUI.doubleSpin(box1, self, "selectPValue", 
                                               min=2e-7, max=1.0, step=1e-7, 
                                               label="P-value:"))
        check.disables.append(OWGUI.button(box1, self, "Select", callback=self.select_p_best))
        check.makeConsistent()

        box1 = OWGUI.widgetBox(box, orientation='horizontal')
        OWGUI.spin(box1, self, "selectNBest", 0, 10000, step=1, label="Best Ranked:")
        OWGUI.button(box1, self, "Select", callback=self.select_n_best)

        box = OWGUI.widgetBox(self.controlArea, "Output")
        b = OWGUI.button(box, self, "&Commit", callback=self.commit)
        cb = OWGUI.checkBox(box, self, "auto_commit", "Commit on change")
        OWGUI.setStopper(self, b, cb, "data_changed_flag", self.commit)
        OWGUI.checkBox(box, self, "add_scores_to_output", "Add gene scores to output",
                       callback=self.commit_if) 
        
        OWGUI.rubber(self.controlArea)

        self.connect(self.graphButton, SIGNAL("clicked()"), self.histogram.saveToFile)
        
        self.loadSettings()

        self.data = None
        self.discData = None
        self.scoreCache = {}
        self.nullDistCache = {}
        self.cuts = {}
        self.null_dist = []
        self.targets = []
        self.scores = {}
        self.genes_in_columns = True
        self.target_selections = None
        
        self.on_scoring_method_changed()
        
        self.resize(800, 600)
        
    def clear(self):
        """ Clear widget.
        """
        self.scoreCache = {}
        self.nullDistCache = {}
        self.discData = None
        self.data = None
        self.attribute_targets = []
        self.class_targets = []
        self.label_selection_widget.clear()
        self.clear_plot()
        
    def clear_plot(self):
        """ Clear the histogram plot.
        """
        self.histogram.removeDrawingCurves()
        self.histogram.clear()
        self.histogram.replot()
        
    def init_from_data(self, data):
        """ Init widget state from the data.
        """
        if data:
            items = [attr.attributes.items() for attr in data.domain.attributes]
            items = reduce(add, items, [])
            
            targets = defaultdict(set)
            for label, value in items:
                targets[label].add(value)
                
            targets = [(key, sorted(vals)) for key, vals in targets.items() \
                       if len(vals) >= 2]
            self.attribute_targets = targets
            
            class_var = data.domain.class_var
            if class_var and isinstance(class_var, orange.EnumVariable):
                targets = [(class_var.name,
                            list(class_var.values))]
                self.class_targets = targets
            else:
                self.class_targets = []
            
    def update_targets_widget(self):
        """ Update the contents of the targets widget.
        """
        if self.data:
            if self.genes_in_columns:
                targets = self.attribute_targets
            elif self.data.domain.classVar:
                targets = self.class_targets
            else:
                targets = []
        else:
            targets = []
            
        self.label_selection_widget.clear()
        self.label_selection_widget.set_labels(targets)
        self.data_labels = targets

    def set_data(self, data):
        self.closeContext("Data")
        self.closeContext("TargetSelection")
        self.error([0, 1])
        self.warning(0)
        self.clear()
        self.genes_in_columns_check.setEnabled(True)
        self.data = data
        self.init_from_data(data)

        if self.data:
            self.genes_in_columns = not data_hints.get_hint(data, "genesinrows", False)
            self.openContext("Data", data)

            # If only attr. labels or only class values then disable
            # the 'Genes in columns' control
            if not self.attribute_targets or not self.class_targets:
                self.genes_in_columns_check.setEnabled(False)
                self.genes_in_columns = bool(self.attribute_targets)

        self.update_targets_widget()

        if self.data is not None  and \
                not (self.attribute_targets or self.class_targets):
            # If both attr. labels and classes are missing, show an error
            self.error(1, "Cannot compute gene scores! Differential expression widget "
                          "requires a data-set with a discrete class variable "
                          "or attribute labels!")
            self.data = None

        if self.data:
            # Load context selection
            items = [(label, v) for label, values in self.data_labels for v in values]

            self.target_selections = [values[:1] for _, values in self.data_labels]
            label, values = self.data_labels[0]
            self.current_target_selection = label, values[:1] # Default selections

            self.openContext("TargetSelection", set(items)) # Load selections from context
            self.label_selection_widget.set_selection(*self.current_target_selection)

        if not self.data:
            self.send("Example table with selected genes", None)
            self.send("Example table with remaining genes", None)
            self.send("Selected genes", None)

    def set_targets(self, targets):
        """ Set the target groups for score computation.
        """
        self.targets = targets
        self.update_scores()
    
    def compute_scores(self, data, score_func, use_attribute_labels,
                       target=None, advance=lambda: None):
        score_func = score_func(data, use_attribute_labels)
        advance()
        score = score_func(target=target)
        score = [(key, val) for key, val in score if val is not ma.masked]
        return score
    
    def compute_null_distribution(self, data, score_func, use_attributes,
                                  target=None, perm_count=10, advance=lambda: None):
        score_func = score_func(data, use_attributes)
        dist = score_func.null_distribution(perm_count, target, advance=advance)
        return [score for run in dist for k, score in run if score is not ma.masked]
            
    @disable_controls
    def update_scores(self):
        """ Compute the scores and update the histogram.
        """
        self.clear_plot()
        self.error(0)
        label, values = self.current_target_selection
        if not self.data or label is None:
            return
        _, score_func, _, two_sample_test = self.score_methods[self.method_index]
        if two_sample_test:
            target = self.targets
            score_target = set(target)
            ind1, ind2 = score_func(self.data, self.genes_in_columns).test_indices(score_target)
            if not len(ind1) or not len(ind2):
                self.error(0, "Target labels most exclude/include at least one value.")
                return
            
        else: # ANOVA should use all labels. 
            target = dict(self.data_labels)[label]
            if self.genes_in_columns:
                target = [(label, t) for t in target]
            score_target = target
#            indices = score_func(self.data, self.genes_in_columns).test_indices(score_target)
            # TODO: Check that each label has more than one measurement, raise warning otherwise.
         
        pb = OWGUI.ProgressBar(self, 4 + self.permutations_count if self.compute_null else 3)
        self.scores = dict(self.compute_scores(self.data, score_func,
                    self.genes_in_columns, score_target, advance=pb.advance))
        pb.advance()
        if self.compute_null:
            self.null_dist = self.compute_null_distribution(self.data,
                        score_func, self.genes_in_columns, score_target,
                        self.permutations_count, advance=pb.advance)
        else:
            self.null_dist = []
        pb.advance()
        htype = self.histType[self.score_methods[self.method_index][2]]
        score_type = self.score_methods[self.method_index][0]
        self.histogram.type = htype
        if self.scores:
            self.histogram.setValues(self.scores.values())
            low, high = self.thresholds.get(score_type, (float("-inf"), float("inf")))
            minx, maxx = self.histogram.minx, self.histogram.maxx
            low, high = max(low, minx), min(high, maxx)

            if htype == "hiTail":
                low = high
            if htype == "lowTail":
                high = low

            self.histogram.setBoundary(low, high)

            if self.compute_null and self.null_dist:
                nullY, nullX = numpy.histogram(self.null_dist, bins=self.histogram.xData)
                nullY = nullY / self.permutations_count
                self.histogram.nullCurve = self.histogram.addCurve(
                    "nullCurve", Qt.black, Qt.black, 6,
                    symbol=QwtSymbol.NoSymbol, style=QwtPlotCurve.Steps,
                    xData=nullX, yData=nullY)

                minx = min(min(nullX), minx)
                maxx = max(max(nullX), maxx)
                miny = min(min(nullY), self.histogram.miny)
                maxy = max(max(nullY), self.histogram.maxy)
                spanx, spany = maxx - minx, maxy - miny
                self.histogram.setAxisScale(
                    QwtPlot.xBottom, minx - 0.05 * spanx, maxx + 0.05 * spanx)
                self.histogram.setAxisScale(
                    QwtPlot.yLeft, miny - 0.05 * spany, maxy + 0.05 * spany)

            state = dict(hiTail=(False, True), lowTail=(True, False), twoTail=(True, True))
            for spin, visible in zip((self.upperBoundarySpin, self.lowerBoundarySpin), state[self.histogram.type]):
                spin.setVisible(visible)

            # If this is a two sample test add markers to the left and right
            # plot indicating which target group is over-expressed in that
            # part
            if self.method_index in [0, 2, 6]:
                if self.method_index == 0: ## fold change is centered on 1.0
                    x1, y1 = (self.histogram.minx + 1) / 2 , self.histogram.maxy
                    x2, y2 = (self.histogram.maxx + 1) / 2 , self.histogram.maxy
                else:
                    x1, y1 = (self.histogram.minx) / 2 , self.histogram.maxy
                    x2, y2 = (self.histogram.maxx) / 2 , self.histogram.maxy
                if self.genes_in_columns:
                    label = target[0][0]
                    target_values = [t[1] for t in target]
                    values = dict(self.data_labels)[label]
                else:
                    target_values = target
                    values = self.data_labels[0][1]
                    
                left = ", ".join(v for v in values if v not in target_values)
                right = ", ".join(v for v in values if v in target_values)
                    
                self.histogram.addMarker(left, x1, y1)
                self.histogram.addMarker(right, x2, y2)
            self.warning(0)
        else:
            self.warning(0, "No scores obtained.")
        self.histogram.replot()
        pb.advance()
        pb.finish()
        self.update_data_info_label()

    def on_boundary_change(self, low, high):
        stype = self.score_methods[self.method_index][0]
        self.thresholds[stype] = (low, high)
        self.update_selected_info_label(low, high)
        self.commit_if()

    def update_data_info_label(self):
        if self.data:
            samples, genes = len(self.data), len(self.data.domain.attributes)
            if self.genes_in_columns:
                samples, genes = genes, samples
                target_labels = [t[1] for t in self.targets]
            else:
                target_labels = self.targets
            text = "%i samples, %i genes\n" % (samples, genes)
            text += "Sample target: '%s'" % (",".join(target_labels))
        else:
            text = "No data on input\n"
        self.dataInfoLabel.setText(text)

    def update_selected_info_label(self, cutOffLower=0, cutOffUpper=0):
        self.cuts[self.method_index] = (cutOffLower, cutOffUpper)
        if self.data:
            scores = np.array(self.scores.values())
            test = self.score_methods[self.method_index][2]
            self.selectedInfoLabel.setText("%i selected genes" % len(np.nonzero(test(scores, cutOffLower, cutOffUpper))[0]))
        else:
            self.selectedInfoLabel.setText("0 selected genes")

    def select_n_best(self):
        scores = self.scores.items()
        scores.sort(lambda a,b:cmp(a[1], b[1]))
        if not scores:
            return
        if self.score_methods[self.method_index][2]==self.oneTailTestHi:
            scores = scores[-max(self.selectNBest, 1):]
            self.histogram.setBoundary(scores[0][1], scores[0][1])
        elif self.score_methods[self.method_index][2]==self.oneTailTestLow:
            scores = scores[:max(self.selectNBest,1)]
            self.histogram.setBoundary(scores[-1][1], scores[-1][1])
        else:
            scoresHi = scores[-max(min(self.selectNBest, len(scores)/2), 1):]
            scoresLo = scores[:max(min(self.selectNBest, len(scores)/2), 1)]
            scores = [(abs(score), 1) for attr, score in scoresHi] + [(abs(score), -1) for attr, score in scoresLo]
            if self.score_methods[self.method_index][0]=="fold change": ## comparing fold change on a logaritmic scale
                scores =  [(abs(math.log(max(min(score, 1e300), 1e-300), 2.0)), sign) for score, sign in scores]
            scores.sort()
            scores = scores[-max(self.selectNBest, 1):]
            countHi = len([score for score, sign in scores if sign==1])
            countLo = len([score for score, sign in scores if sign==-1])
            cutHi = scoresHi[-countHi][1] if countHi else scoresHi[-1][1] + 1e-7
            cutLo = scoresLo[countLo-1][1] if countLo else scoresLo[0][1] - 1e-7
            self.histogram.setBoundary(cutLo, cutHi)

    def update_boundary(self):
        if self.data != None and self.scores.items():
            if self.score_methods[self.method_index][2]==self.oneTailTestHi:
                self.histogram.setBoundary(self.histogram.lowerBoundary, self.histogram.lowerBoundary)
            elif self.score_methods[self.method_index][2]==self.oneTailTestLow:
                self.histogram.setBoundary(self.histogram.upperBoundary, self.histogram.upperBoundary)
            else:
                self.histogram.setBoundary(self.histogram.lowerBoundary, self.histogram.upperBoundary)

    def select_p_best(self):
        if not self.null_dist:
            return
        nullDist = sorted(self.null_dist)
        test = self.score_methods[self.method_index][2]
        count = min(int(len(nullDist)*self.selectPValue), len(nullDist))
        if test == self.oneTailTestHi:
            cut = nullDist[-count] if count else nullDist[-1] # + 1e-7
            self.histogram.setBoundary(cut, cut)
        elif test == self.oneTailTestLow:
            cut = nullDist[count - 1] if count else nullDist[0] # - 1e-7
            self.histogram.setBoundary(cut, cut)
        elif count:
            scoresHi = nullDist[-count:]
            scoresLo = nullDist[:count]
            scores = [(abs(score), 1) for score in scoresHi] + [(abs(score), -1) for score in scoresLo]
            if self.score_methods[self.method_index][0] == "fold change": ## fold change is on a logaritmic scale
                scores =  [(abs(math.log(max(min(score, 1e300), 1e-300), 2.0)), sign) for score, sign in scores]
            scores = sorted(scores)[-count:]
            countHi = len([score for score, sign in scores if sign==1])
            countLo = len([score for score, sign in scores if sign==-1])
            cutHi = scoresHi[-countHi] if countHi else scoresHi[-1] + 1e-7
            cutLo = scoresLo[countLo-1] if countLo else scoresLo[0] - 1e-7
            self.histogram.setBoundary(cutLo, cutHi)
        else:
            self.histogram.setBoundary(nullDist[0] - 1e-7, nullDist[-1] + 1e-7)
            

    def commit_if(self):
        if self.auto_commit:
            self.commit()
        else:
            self.data_changed_flag = True
        
    def commit(self):
        if not self.data or not self.scores:
            return
        test = self.score_methods[self.method_index][2]
        
        cutOffUpper = self.histogram.upperBoundary
        cutOffLower = self.histogram.lowerBoundary
        
        scores = np.array(self.scores.items())
        scores[:, 1] = test(np.array(scores[:, 1], dtype=float), cutOffLower, cutOffUpper)
        selected = set([key for key, test in scores if test])
        remaining = set([key for key, test in scores if not test])
        if self.data and self.genes_in_columns:
            selected = sorted(selected)
            if selected:
                newdata = orange.ExampleTable(
                    orange.Domain(self.data.domain),
                    [self.data[int(i)] for i in selected],
                    name=self.data.name
                )
            else:
                newdata = None
            if self.add_scores_to_output:
                score_attr = orange.FloatVariable(self.score_methods[self.method_index][0])
                mid = orange.newmetaid()
                
            if self.add_scores_to_output and newdata is not None:
                newdata.domain.addmeta(mid, score_attr)
                for ex, key in zip(newdata, selected):
                    ex[mid] = self.scores[key]
                    
            self.send("Example table with selected genes", newdata)
            
            remaining = sorted(remaining)
            if remaining:
                newdata = orange.ExampleTable(
                    orange.Domain(self.data.domain),
                    [self.data[int(i)] for i in remaining],
                    name=self.data.name
                )
            else:
                newdata = None
            
            if self.add_scores_to_output and newdata is not None:
                newdata.domain.addmeta(mid, score_attr)
                for ex, key in zip(newdata, remaining):
                    ex[mid] = self.scores[key]
                    
            self.send("Example table with remaining genes", newdata)
            
        elif self.data and not self.genes_in_columns:
            method_name = self.score_methods[self.method_index][0]
            selected_attrs = [attr for attr in self.data.domain.attributes
                              if attr in selected or
                              attr.varType == orange.VarTypes.String]  # ?? why strings
            if self.add_scores_to_output:
                scores = [self.scores[attr] for attr in selected_attrs]
                attrs = [copy_descriptor(attr) for attr in selected_attrs]
                for attr, score in zip(attrs, scores):
                    attr.attributes[method_name] = str(score)
                selected_attrs = attrs

            newdomain = orange.Domain(selected_attrs, self.data.domain.classVar)
            newdomain.addmetas(self.data.domain.getmetas())
            newdata = orange.ExampleTable(
                newdomain, self.data, name=self.data.name
            )
            self.send("Example table with selected genes",
                      newdata if selected_attrs else None)

            remaining_attrs = [attr for attr in self.data.domain.attributes
                               if attr in remaining]
            if self.add_scores_to_output:
                scores = [self.scores[attr] for attr in remaining_attrs]
                attrs = [copy_descriptor(attr) for attr in remaining_attrs]
                for attr, score in zip(attrs, scores):
                    attr.attributes[method_name] = str(scores)
                remaining_attrs = attrs

            newdomain = orange.Domain(remaining_attrs, self.data.domain.classVar)
            newdomain.addmetas(self.data.domain.getmetas())
            newdata = orange.ExampleTable(
                newdomain, self.data, name=self.data.name
            )
            self.send("Example table with remaining genes",
                      newdata if remaining_attrs else None)

            domain = orange.Domain([orange.StringVariable("label"),
                                    orange.FloatVariable(self.score_methods[self.method_index][0])],
                                    False)
            if selected_attrs:
                selected_genes = orange.ExampleTable(domain,
                            [[attr.name, self.scores.get(attr, 0)] for attr in selected_attrs])
            else:
                selected_genes = None
            self.send("Selected genes",  selected_genes)
            
        else:
            self.send("Example table with selected genes", None)
            self.send("Example table with remaining genes", None)
            self.send("Selected genes", None)
        self.data_changed_flag = False

    def on_target_changed(self):
        label, values = self.label_selection_widget.current_selection()

        if values is None:
            values = []

        if self.genes_in_columns:
            targets = [(label, t) for t in values]
        else:
            targets = values

        self.targets = targets
        self.current_target_selection = label, values
        # Save target label selection
        labels = [l for l, _ in self.data_labels]
        if label in labels:
            label_index = labels.index(label)
            self.target_selections[label_index] = values
        self.set_targets(targets)

    def on_label_activated(self, index):
        selection = self.target_selections[index]
        if not selection:
            selection = self.data_labels[index][1][:1]
        self.label_selection_widget.set_selection(index, selection)

    def on_genes_in_columns_change(self):
        self.closeContext("TargetSelection")
        self.update_targets_widget()
        items = [(label, v) for label, values in self.data_labels \
                                for v in values]
        if self.data_labels:
            label, values = self.data_labels[0] 
            self.current_target_selection = label, values[:1]
            self.target_selections = [values[1:] for _, values in self.data_labels]
            self.openContext("TargetSelection", set(items))
            self.label_selection_widget.set_selection(*self.current_target_selection)
        
    def on_scoring_method_changed(self):
        two_sample = self.score_methods[self.method_index][3]
        self.label_selection_widget.values_view.setDisabled(not two_sample)
        
    def settingsFromWidgetCallbackTargetSelection(self, handler, context):
        context.target_selections = self.target_selections
        context.current_target_selection = self.current_target_selection
        
    def settingsToWidgetCallbackTargetSelection(self, handler, context):
        self.target_selections = getattr(context, "target_selections", self.target_selections)
        self.current_target_selection = getattr(context, "current_target_selection", self.current_target_selection)


def copy_descriptor(desc):
    d = desc.clone()
    d.attributes = dict(desc.attributes)
    c = Orange.core.ClassifierFromVar()
    c.whichVar = desc
    d.get_value_from = c
    return d


if __name__=="__main__":
    import sys
    app = QApplication(sys.argv)
    data = orange.ExampleTable(os.path.expanduser("~/Documents/GDS636"))
    w = OWFeatureSelection()
    w.show()
    w.set_data(data)
    app.exec_()
    w.saveSettings()
