import sys
import os

from functools import partial, reduce

import numpy
from AnyQt.QtGui import QPainter
from AnyQt.QtCore import QObject, QSize, QThread, QThreadPool, Slot

import pyqtgraph as pg

import Orange.data

from Orange.widgets import widget, gui, settings
from Orange.widgets.utils import concurrent

from ..utils import expression


def is_discrete(var):
    return isinstance(var, Orange.data.DiscreteVariable)


def group_mask_columns(data, key, values):
    """
    Return a boolean array mask of data columns (variables).

    The mask will be True wherever the `var.attributes[key]` contains
    one of `values`.

    Parameters
    ----------
    data : Orange.data.Table
        Source data table.
    key : str
        The variable label key (where `values` are looked for)
    values : sequence of str
        The values (for the corresponding `key`) selecting the columns.
    """
    target = set([(key, value) for value in values])
    mask = [not target.isdisjoint(var.attributes.items())
            for var in data.domain.attributes]
    return numpy.array(mask, dtype=bool)


def group_mask_rows(data, var, values):
    """
    Return a boolean array mask for data rows (instances).

    The mask will be True wherever the row's entry for `var` contains
    one of the `values`.

    Parameters
    ----------
    data : Orange.data.Table
        Source data table.
    var : Orange.data.DiscreteVariable
        The variable/column on which to match `values`.
    values : sequence of str
        The values to select (must be a subset of `var.values`)
    """
    var = data.domain[var]
    col_view, _ = data.get_column_view(var)
    target_ind = [var.values.index(t) for t in values]

    mask = numpy.zeros_like(col_view, dtype=bool)
    for i in target_ind:
        mask |= col_view == i

    return mask


def group_mask(data, key, values, axis=1):
    if axis == 1:
        return group_mask_columns(data, key, values)
    elif axis == 0:
        return group_mask_rows(data, key, values)
    else:
        raise ValueError("0 <= axis < 2")


def table_take(data, indices, axis=0):
    if axis == 0:
        return data[indices]
    elif axis == 1:
        return data[:, indices]


class ScatterPlotItem(pg.ScatterPlotItem):
    def paint(self, painter, option, widget):
        if self.opts["antialias"]:
            painter.setRenderHint(QPainter.Antialiasing, True)
        if self.opts["pxMode"]:
            painter.setRenderHint(QPainter.SmoothPixmapTransform, True)
        super().paint(painter, option, widget)


class ProgressBarDiscard(QObject):

    def __init__(self, parent, redirect):
        QObject.__init__(self, parent)
        self.redirect = redirect
        self._delay = False

    @Slot(float)
    def progressBarSet(self, value):
        # Discard OWBaseWidget.progressBarSet call, because it calls qApp.processEvents
        # which can result in 'event queue climbing' and max. recursion error if GUI thread
        # gets another advance signal before it finishes with this one
        if not self._delay:
            try:
                self._delay = True
                self.redirect.progressBarSet(value)
            finally:
                self._delay = False


def withexcepthook(func):
    def wrapped(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except BaseException as ex:
            sys.excepthook(*sys.exc_info())
            raise

    return wrapped


class OWMAPlot(widget.OWWidget):
    name = "MA Plot"
    description = "Normalize expression array data on a MA - plot"
    icon = "../widgets/icons/MAPlot.svg"
    priority = 5000

    inputs = [("Expression array", Orange.data.Table, "setData")]
    outputs = [("Normalized expression array", Orange.data.Table,
                widget.Default),
               ("Filtered expression array", Orange.data.Table)]

    CENTER_METHODS = [("Average", expression.MA_center_average),
                      ("Lowess (fast - interpolated)", expression.MA_center_lowess_fast),
                      ("Lowess", expression.MA_center_lowess)]

    MERGE_METHODS = [("Average", numpy.ma.average),
                     ("Median", numpy.ma.median),
                     ("Geometric mean", expression.geometric_mean)]

    settingsHandler = settings.DomainContextHandler()

    appendZScore = settings.Setting(False)
    appendRIValues = settings.Setting(False)

    selectedGroup = settings.ContextSetting(0)
    selectedCenterMethod = settings.Setting(0)
    selectedMergeMethod = settings.Setting(0)
    zCutoff = settings.Setting(1.96)

    autoCommit = settings.Setting(False)

    def __init__(self, parent=None):
        super().__init__(parent)

        # GUI
        box = gui.widgetBox(self.controlArea, "Info", addSpace=True)
        self.infoBox = gui.widgetLabel(box, "No data on input.")

        box = gui.widgetBox(self.controlArea, "Split by", addSpace=True)
        self.groupCombo = gui.comboBox(
            box, self, "selectedGroup", callback=self.onGroupSelection)

        gui.comboBox(self.controlArea, self, "selectedCenterMethod",
                     box="Center Fold-change Using",
                     items=[name for name, _ in self.CENTER_METHODS],
                     callback=self.onCenterMethodChange,
                     addSpace=True)

        gui.comboBox(self.controlArea, self, "selectedMergeMethod",
                     box="Merge Replicates",
                     items=[name for name, _ in self.MERGE_METHODS],
                     tooltip="Select the method for replicate merging",
                     callback=self.onMergeMethodChange,
                     addSpace=True)

        box = gui.doubleSpin(self.controlArea, self, "zCutoff", 0.0, 3.0, 0.01,
                             box="Z-Score Cutoff",
                             callback=[self.replotMA, self.commitIf])

        gui.separator(self.controlArea)

        box = gui.widgetBox(self.controlArea, "Ouput")
        gui.checkBox(box, self, "appendZScore", "Append Z-Scores",
                     tooltip="Append calculated Z-Scores to output",
                     callback=self.commitIf)

        gui.checkBox(box, self, "appendRIValues",
                     "Append Log Ratio and Intensity values",
                     tooltip="Append calculated Log Ratio and Intensity "
                             "values to output data",
                     callback=self.commitIf)

        gui.rubber(self.controlArea)

        gui.auto_commit(self.controlArea, self, "autoCommit", "Commit")

        self.graph = pg.PlotWidget(background="w")
        self.graph.getAxis("bottom").setLabel("Intensity: log<sub>10</sub>(R*G)")
        self.graph.getAxis("left").setLabel("Log ratio: log<sub>2</sub>(R/G)")

        self.mainArea.layout().addWidget(self.graph)
        self.groups = []
        self.split_data = None, None
        self.merged_splits = None, None
        self.centered = None, None
        self.changedFlag = False
        self.data = None

        self._executor = concurrent.ThreadExecutor(
            threadPool=QThreadPool(maxThreadCount=1))

    def sizeHint(self):
        return QSize(800, 600)

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
        col_labels = [attr.attributes.items()
                      for attr in data.domain.attributes]
        col_labels = sorted(reduce(set.union, col_labels, set()))
        col_labels = [(key, value, 1) for key, value in col_labels]

        attrs = [attr for attr in data.domain.variables + data.domain.metas
                 if is_discrete(attr)]

        row_labels = [(attr.name, value, 0)
                      for attr in attrs for value in attr.values]

        def filterSingleValues(labels):
            ret = []
            for name, value, axis in labels:
                match = [(n, v, a) for n, v, a in labels if n == name]
                if len(match) > 1:
                    ret.append((name, value, axis))
            return ret

        col_labels = filterSingleValues(col_labels)
        row_labels = filterSingleValues(row_labels)

        return col_labels + row_labels

    def setData(self, data):
        self.closeContext()
        self.clear()
        self.unconditional_commit()
        self.data = data
        self.error([0, 1])
        if data is not None:
            self.infoBox.setText("%i genes on input" % len(data))
            self.groups = self.proposeGroups(data)
            self.groupCombo.clear()
            self.groupCombo.addItems(["%s: %s" % (key, value)
                                      for key, value, axis in self.groups])

            if not self.groups:
                self.error(1, "Input data has no class attribute or attribute labels!")
                self.clear()
                return

            self.openContext(data)
            self.selectedGroup = min(self.selectedGroup, len(self.groups) - 1)

            self.updateInfoBox()
            self.splitData()
            self.runNormalization()
        else:
            self.clear()

    def clear(self):
        self.error([0, 1])
        self.data = None
        self.groups = []
        self.centered = None, None
        self.split_data = None, None
        self.merged_splits = None, None
        self.groupCombo.clear()

        self.graph.clear()
        self.infoBox.setText("No data on input")

    def updateInfoBox(self):
        # TODO: report on number of columns or rows depending
        # on the genes 'axis' 
        self.infoBox.setText("%i genes on input" % len(self.data))

    def getSelectedGroup(self):
        return self.groups[self.selectedGroup]

    def getSelectedGroupSplit(self):
        key, value, axis = self.getSelectedGroup()
        other_values = [v for k, v, a in self.groups
                        if k == key and a == axis and v != value]
        return [(key, [value]), (key, other_values)], axis

    def getGeneNames(self):
        key, value, axis = self.getSelectedGroup()
        if axis == 0:
            genes = [str(ex[key]) for ex in self.data]
        else:
            genes = [attr.name for attr in self.data.domain.attributes]

        return genes

    def splitData(self):
        groups, axis = self.getSelectedGroupSplit()
        self.split_ind = [expression.select_indices(self.data, key, value, axis)
                          for key, value in groups]
        self.split_ind = [numpy.flatnonzero(
                              group_mask(self.data, key, value, axis))
                          for key, value in groups]

        self.split_data = [table_take(self.data, indices, axis)
                           for indices in self.split_ind]

    def getMerged(self):
        split1, split2 = self.split_data
        array1 = numpy.ma.array(split1.X, copy=False, mask=numpy.isnan(split1.X))
        array2 = numpy.ma.array(split2.X, copy=False, mask=numpy.isnan(split2.X))

        _, _, axis = self.getSelectedGroup()
        merge_function = self.MERGE_METHODS[self.selectedMergeMethod][1]

        merged1 = expression.merge_replicates(
            array1, axis, merge_function=merge_function)
        merged2 = expression.merge_replicates(
            array2, axis, merge_function=merge_function)
        self.merged_splits = merged1, merged2

        return self.merged_splits

    def runNormalization(self):
        self.progressBarInit()
        self.progressBarSet(0.0)
        G, R = self.getMerged()
        self.progressBarSet(5.0)

        center_method = self.CENTER_METHODS[self.selectedCenterMethod][1]

        # TODO: progess bar , lowess can take a long time
        if self.selectedCenterMethod in [1, 2]:  # Lowess
            Gc, Rc = center_method(G, R, f=1. / min(500., len(G) / 100), iter=1)
        else:
            Gc, Rc = center_method(G, R)
        self.progressBarSet(70.0)
        self.centered = Gc, Rc
        self.z_scores = expression.MA_zscore(Gc, Rc, 1. / 3.)
        self.progressBarSet(100.0)
        self.plotMA(Gc, Rc, self.z_scores, self.zCutoff)
        self.progressBarFinished()

    def runNormalizationAsync(self):
        """ Run MA centering and z_score estimation in a separate thread 
        """
        self.error(0)
        self.progressBarInit(processEvents=None)
        self.progressBarSet(0.0, processEvents=None)
        G, R = self.getMerged()

        center_method = self.CENTER_METHODS[self.selectedCenterMethod][1]
        use_lowess = self.selectedCenterMethod in [1, 2]

        @withexcepthook
        def run(progressCallback=lambda value: None):
            if use_lowess:
                Gc, Rc = center_method(
                    G, R, f=2. / 3., iter=1,
                    progressCallback=lambda val: progressCallback(val / 2))
            else:
                Gc, Rc = center_method(G, R)
            progressCallback(50)
            z_scores = expression.MA_zscore(
                Gc, Rc, 1. / 3.,
                progressCallback=lambda val: progressCallback(50 + val / 2))

            return Gc, Rc, z_scores

        self.progressDiscard = ProgressBarDiscard(self, self)

        progress = concurrent.methodinvoke(
            self.progressDiscard, "progressBarSet", (float,))

        self._task = concurrent.Task(function=partial(run, progress))
        self._task.resultReady.connect(self.onResultsReady)
        self._task.exceptionReady.connect(self.onException)

        self.setEnabled(False)
        self.setBlocking(True)

        self._executor.submit(self._task)

    # comment out this line if threading creates any problems
    runNormalization = runNormalizationAsync

    def onResultsReady(self, results):
        """Handle the results of centering and z-scoring
        """
        assert(QThread.currentThread() is self.thread())
        Gc, Rc, z_scores = results

        self.setEnabled(True)
        self.setBlocking(False)
        self.progressBarFinished()
        self.centered = Gc, Rc
        self.z_scores = z_scores
        self.plotMA(Gc, Rc, z_scores, self.zCutoff)
        self.commitIf()

    def onException(self, error):
        print("Error: ", error, file=sys.stderr)
        self.error(0, "Error: {!s}".format(error))
        self.setEnabled(True)
        self.setBlocking(False)

    def plotMA(self, G, R, z_scores, z_cuttof):
        ratio, intensity = expression.ratio_intensity(G, R)

        validmask = (numpy.isfinite(ratio) &
                     numpy.isfinite(intensity) &
                     numpy.isfinite(z_scores))
        for array in [ratio, intensity, z_scores]:
            if numpy.ma.is_masked(array):
                validmask &= ~array.mask

        filtered_ind = numpy.ma.where(validmask)
        ratio = numpy.take(ratio, filtered_ind)
        intensity = numpy.take(intensity, filtered_ind)
        z_scores = numpy.take(z_scores, filtered_ind)

        red_ind = numpy.ma.where(numpy.ma.abs(z_scores) >= z_cuttof)
        blue_ind = numpy.ma.where(numpy.ma.abs(z_scores) < z_cuttof)

        red_xdata, red_ydata = intensity[red_ind], ratio[red_ind]
        blue_xdata, blue_ydata = intensity[blue_ind], ratio[blue_ind]

        self.graph.clear()
        line = pg.InfiniteLine(pos=(0, 0), angle=0, pen=pg.mkPen((0, 0, 0)))
        self.graph.addItem(line)

        red_points = ScatterPlotItem(
            x=red_xdata, y=red_ydata,
            brush=pg.mkBrush((255, 0, 0, 100)), size=6, antialias=True)
        blue_points = ScatterPlotItem(
            x=blue_xdata, y=blue_ydata,
            brush=pg.mkBrush((0, 0, 255, 100)), size=6, antialias=True)

        self.graph.addItem(red_points)
        self.graph.addItem(blue_points)

    def replotMA(self):
        if self.data and self.centered:
            Gc, Rc = self.centered
            self.plotMA(Gc, Rc, self.z_scores, self.zCutoff)

    def commitIf(self):
        self.commit()

    def commit(self):
        if not self.data:
            self.send("Normalized expression array", None)
            self.send("Filtered expression array", None)
            return

        G, R = self.merged_splits
        Gc, Rc = self.centered
        ind1, ind2 = self.split_ind

        gfactor = Gc / G
        domain = self.data.domain
        newmetas = []
        M = []

        _, _, axis = self.getSelectedGroup()
        if self.appendZScore and axis == 1:
            attr = Orange.data.ContinuousVariable("Z-Score")
            newmetas.append(attr)
            M.append(self.z_scores.filled(numpy.nan))

        if self.appendRIValues and axis == 1:
            r_attr = Orange.data.ContinuousVariable("Log Ratio")
            i_attr = Orange.data.ContinuousVariable("Intensity")
            ratio, intensity = expression.ratio_intensity(Gc, Rc)
            newmetas.extend([r_attr, i_attr])
            M.extend([ratio.filled(numpy.nan),
                      intensity.filled(numpy.nan)])

        if newmetas:
            domain = Orange.data.Domain(
                self.data.domain.attributes, self.data.domain.class_vars,
                self.data.domain.metas + tuple(newmetas))

        data = Orange.data.Table.from_table(domain, self.data)
        data.ensure_copy()

        if axis == 0:
            data.X[ind1, :] *= gfactor.reshape((1, -1))
        else:
            data.X[:, ind1] *= gfactor.reshape((-1, 1))

        for i, mcol in enumerate(reversed(M)):
            data.metas[:, -i - 1] = mcol

        selected_indices = numpy.flatnonzero(
            numpy.abs(self.z_scores.filled(0)) >= self.zCutoff)

        if axis == 0:
            attrs = [data.domain[i] for i in selected_indices]
            domain = Orange.data.Domain(attrs, data.domain.class_vars,
                                        data.domain.metas)
            filtered_data = data.from_table(domain, data)
        else:
            filtered_data = data[selected_indices]

        self.send("Normalized expression array", data)
        self.send("Filtered expression array", filtered_data)


def test_main(argv=sys.argv):
    from AnyQt.QtWidgets import QApplication
    if len(argv) > 1:
        filename = argv[1]
    else:
        filename = os.path.expanduser("~/GDS1210.tab")
    data = Orange.data.Table(filename)

    app = QApplication(argv)
    w = OWMAPlot()
    w.setData(data)
    w.show()
    w.raise_()
    r = app.exec_()
    w.saveSettings()
    return r

if __name__ == "__main__":
    sys.exit(test_main())
