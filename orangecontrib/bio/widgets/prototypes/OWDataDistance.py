## Automatically adapted for numpy.oldnumeric Oct 04, 2007 by 

"""
<name>Data Distance</name>
<description>Computes a distance matrix between data files.</description>
<icon>icons/ChipDistance.png</icon>
<priority>1160</priority>
<contact>Peter Juvan (peter.juvan@fri.uni-lj.si)</contact>
<prototype>1</prototype>
"""

from __future__ import absolute_import

from qt import *
from qtcanvas import *

import numpy.oldnumeric as Numeric, numpy.oldnumeric.ma as MA

import orange, statc
from Orange.OrangeWidgets import OWGUI
from Orange.OrangeWidgets.OWWidget import *

from .OWDataFiles import DataFiles

import warnings
warnings.filterwarnings("ignore", "'strain'", orange.AttributeWarning)
warnings.filterwarnings("ignore", "'dirname'", orange.AttributeWarning)

##############################################################################
# main class

class OWDataDistance(OWWidget):
    settingsList = ["Metrics"]

    def __init__(self, parent=None, signalManager = None):
        OWWidget.__init__(self, parent, signalManager, 'Data Distance') 
        
        self.inputs = [("Structured Data", DataFiles, self.chipdata)]
        self.outputs = [("Distance Matrix", orange.SymMatrix)]

        self.Metrics = 0
        self.loadSettings()
        self.data = []
##        self.metrics = [("Euclidean", orange.ExamplesDistanceConstructor_Euclidean),
##            ("Manhattan", orange.ExamplesDistanceConstructor_Manhattan),
##            ("Hamming", orange.ExamplesDistanceConstructor_Hamming)]
        self.metrics = [("Manhattan", distManhattan), ("Euclidean", distEuclidean), ("1 - (Pearson correlation coefficient)", distPearson), ("1 - (Spearman rank correlation coefficient)", distSpearman)]

        # GUI
        self.mainArea.setFixedWidth(0)
        # Info box
        box = QVGroupBox("Info", self.controlArea)
        self.infoa = QLabel('No data on input.', box)
        self.infob = QLabel('', box)
        OWGUI.separator(self.controlArea)

        # Distance metrics selection
        items = [x[0] for x in self.metrics]
        OWGUI.comboBox(self.controlArea, self, "Metrics", box="Distance Metrics", items=items,
            tooltip="Metrics to measure distance between data sets.",
            callback=self.onMetricsChange)

        self.resize(384, 138)
        

    ##########################################################################
    # handling of input/output signals

##    def computeDistance(self, d1, d2, dist):
##        """employs orange to compute distances (slower)
##        """
##        d = 0
##        for i in range(len(d1)):
##            d += dist(d1[i], d2[i])
##        d = d / len(d1)
##        return d

    def computeDistance(self, d1, d2):
        """employs MA to cumpute distances (faster)
        """
        return dist(d1.toNumpyMA("a")[0], d2.toNumpyMA("a")[0])


    def computeMatrix(self):
        if not self.data:
            self.send("Distance Matrix", None)
            return
##        if self.Metrics == 0: # bug in orange, correct (remove normalize) once it is fixed
##            dist = self.metrics[self.Metrics][1](self.data[0], normalize=0)
##        else:
##            dist = self.metrics[self.Metrics][1](self.data[0])            
        matrix = orange.SymMatrix(len(self.data))
        matrix.setattr('items', self.data)
        self.progressBarInit()
        pbStep = 100./(len(self.data)**2/2. - len(self.data)/2.)
        for i in range(len(self.data)-1):
            for j in range(i+1, len(self.data)):
##                matrix[i, j] = self.computeDistance(self.data[i], self.data[j], dist)
                matrix[i, j] = self.metrics[self.Metrics][1](MA.ravel(self.data[i].toNumpyMA("a")[0]), MA.ravel(self.data[j].toNumpyMA("a")[0]))
                self.progressBarAdvance(pbStep)
        self.progressBarFinished()
        self.send("Distance Matrix", matrix)


    def chipdata(self, data):
        self.data = []
        if data:
            self.infob.setText("")
            numFiles = reduce(lambda a,b: a+len(b[1]), data, 0)
            lenSD = len(data)
            self.infoa.setText("%d set%s, total of %d data file%s." % (lenSD, ["","s"][lenSD!=1], numFiles, ["","s"][numFiles!=1]))
            numExamplesList = []
            # construct a list of ExampleTable lengths and a list of attribute names
            for (name, etList) in data:
                for et in etList:
                    setattr(et,"dirname",name)
                    setattr(et,"strain",name)
                    self.data.append(et)
                    numExamplesList.append(len(et))
            if len(self.data)>1:
                # test that files contain the same attributes and equal number of examples
                attrSorted = self.data[0].domain.attributes
                attrSorted.sort()
                numEx = len(self.data[0])
                for et in self.data[1:]:
                    attrSorted2 = et.domain.attributes
                    attrSorted2.sort()
                    if map(lambda x: x.name, attrSorted) != map(lambda x: x.name, attrSorted2):
                        self.data = []
                        self.infob.setText("Error: data files contain different attributes, aborting distance computation.")
                        return
                    if len(et) != numEx:
                        self.data = []
                        self.infob.setText("Error: data files contain unequal number of examples, aborting distance computation.")
                        return
                # compute distances
                pb = OWGUI.ProgressBar(self, iterations=len(self.data))
                self.computeMatrix()
                pb.finish()

            else:
                self.data = []
                self.infob.setText('Error: not enough data, aborting distance computation.')
        else:
            self.infoa.setText('No data on input.')


    def onMetricsChange(self):
        if self.data and len(self.data)>1:
            self.computeMatrix()



###########################################################################
# Distance Metrics
###########################################################################

def distManhattan(x,y):
    """normalized Manhattan distance
    """
    x = MA.asarray(x)
    y = MA.asarray(y)
    assert MA.rank(x) == MA.rank(y) == 1
    sumWeights = MA.add.reduce(MA.logical_not(MA.logical_or(MA.getmaskarray(x), MA.getmaskarray(y))).astype(Numeric.Float))
    return MA.add.reduce(MA.absolute(x-y)) / sumWeights


def distManhattanW(x,y,w):
    """normalized weighted Manhattan distance
    """
    x = MA.asarray(x)
    y = MA.asarray(y)
    w = MA.asarray(w)
    assert MA.rank(x) == MA.rank(y) == MA.rank(w) == 1
    sumWeights = MA.add.reduce(w * MA.logical_not(MA.logical_or(MA.getmaskarray(x), MA.getmaskarray(y))).astype(Numeric.Float))
    return MA.add.reduce(w * MA.absolute(x-y)) / sumWeights


def distEuclidean(x,y):
    """normalized euclidean distance
    """
    x = MA.asarray(x)
    y = MA.asarray(y)
    assert MA.rank(x) == MA.rank(y) == 1
    sumWeights = MA.add.reduce(MA.logical_not(MA.logical_or(MA.getmaskarray(x), MA.getmaskarray(y))).astype(Numeric.Float))
    return MA.sqrt(MA.add.reduce((x-y)**2) / sumWeights)


def distEuclideanW(x,y,w):
    """normalized weighted euclidean distance
    """
    x = MA.asarray(x)
    y = MA.asarray(y)
    w = MA.asarray(w)
    assert MA.rank(x) == MA.rank(y) == MA.rank(w) == 1
    sumWeights = MA.add.reduce(w * MA.logical_not(MA.logical_or(MA.getmaskarray(x), MA.getmaskarray(y))).astype(Numeric.Float))
    return MA.sqrt(MA.add.reduce(w * (x-y)**2) / sumWeights)


def distPearson(x,y):
    """distance corresponding to 1 - pearson's correlation coefficient for arrays x,y
    returns distance: 1 - pearson_r
    """
    x = MA.asarray(x)
    y = MA.asarray(y)
    assert MA.rank(x) == MA.rank(y) == 1
    cond = MA.logical_not(MA.logical_or(MA.getmaskarray(x), MA.getmaskarray(y)))
    return 1 - statc.pearsonr(MA.compress(cond,x).tolist(), MA.compress(cond,y).tolist())[0]


def distPearsonW(x,y,w):
    """weighted distance corresponding to 1 - pearson's correlation coefficient for arrays x,y and weights w
    returns distance: 1 - pearson_r
    """
    #TINY = 1.0e-20
    # ones for non-masked places at x,y and w
    x = MA.asarray(x)
    y = MA.asarray(y)
    w = MA.asarray(w)
    assert MA.rank(x) == MA.rank(y) == MA.rank(w) == 1
    mask = MA.logical_or(MA.logical_or(MA.getmaskarray(x), MA.getmaskarray(y)), MA.getmaskarray(w))
    # set mask to w that is equal to the mask from x, y and w
    w = MA.masked_array(w, mask=mask)
    n_w_mean = MA.add.reduce(w)    # n * mean(w)
    x_w = x*w       # x * w
    y_w = y*w       # y * w
    x_wmean = MA.divide(MA.add.reduce(x_w), n_w_mean)     # weighted_mean(x)
    y_wmean = MA.divide(MA.add.reduce(y_w), n_w_mean)     # weighted_mean(x)    
    r_num = MA.add.reduce(x*y*w) - n_w_mean*x_wmean*y_wmean
    r_den = MA.sqrt((MA.add.reduce(x_w*x) - n_w_mean*x_wmean**2) * (MA.add.reduce(y_w*y) - n_w_mean*y_wmean**2))
    return 1 - MA.divide(r_num, r_den)


def distSpearman(x,y):
    """distance corresponding to 1 - spearman's correlation coefficient for arrays x,y
    returns distance: 1 - spearman_r
    """
    x = MA.asarray(x)
    y = MA.asarray(y)
    assert MA.rank(x) == MA.rank(y) == 1
    cond = MA.logical_not(MA.logical_or(MA.getmaskarray(x), MA.getmaskarray(y)))
    return 1 - statc.spearmanr(MA.compress(cond,x).tolist(), MA.compress(cond,y).tolist())[0]

def distSpearmanW(x,y,w):
    """weighted distance corresponding to 1 - spearman's correlation coefficient for arrays x,y and weights w
    returns distance: 1 - spearman_r
    """
    distSpearFunc = _distSpearmanW_NU
    for var in (x,y,w):
        if type(var) == MA.array and MA.count(var) != Numeric.multiply.reduce(var.shape):
            distSpearFunc = _distSpearmanW_MA
            break
    return distSpearFunc(x,y,w)


def _distSpearmanW_NU(x,y,w):
    """x,y,w must be Numeric
    """
    x = Numeric.asarray(x)
    y = Numeric.asarray(y)
    w = Numeric.asarray(w)
    assert Numeric.rank(x) == Numeric.rank(y) == Numeric.rank(w) == 1
    rankx = Numeric.array(statc.rankdata(x.tolist()))
    ranky = Numeric.array(statc.rankdata(y.tolist()))
    return distPearsonW(rankx,ranky,w)


def _distSpearmanW_MA(x,y,w):
    """if any of x,y,w is a MA array containing masked values
    """
    x = MA.asarray(x)
    y = MA.asarray(y)
    w = MA.asarray(w)
    assert MA.rank(x) == MA.rank(y) == MA.rank(w) == 1
    cond = MA.logical_not(MA.logical_or(MA.logical_or(MA.getmaskarray(x), MA.getmaskarray(y)), MA.getmaskarray(w))) 
    # with MA use compress before tolist() !
    rankx = Numeric.array(statc.rankdata(MA.compress(cond, x).tolist()))
    ranky = Numeric.array(statc.rankdata(MA.compress(cond, y).tolist()))
    return distPearsonW(rankx,ranky,MA.compress(cond,w))

###########################################################################
# testing
###########################################################################

if __name__=="__main__":
    from . import OWDataFiles
    from Orange.orng import orngSignalManager
    signalManager = orngSignalManager.SignalManager(0)
    a=QApplication(sys.argv)
    ow=OWDataDistance(signalManager = signalManager)
    signalManager.addWidget(ow)
    a.setMainWidget(ow)
    ow.show()
    ds = OWDataFiles.OWDataFiles(signalManager = signalManager)
    signalManager.addWidget(ds)
    ds.loadData("potato.sub100")
    signalManager.setFreeze(1)
    signalManager.addLink(ds, ow, 'Structured Data', 'Structured Data', 1)
    signalManager.setFreeze(0)
    a.exec_loop()
    ow.saveSettings()
