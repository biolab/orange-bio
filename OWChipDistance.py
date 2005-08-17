"""
<name>Chip Distance</name>
<description>Computes a distance matrix from a set of chips.</description>
<icon>icons/ChipDistance.png</icon>
<priority>1160</priority>
"""

import orange, OWGUI
from qt import *
from qtcanvas import *
from OWWidget import *
from OWChipDataFiles import ChipData

##############################################################################
# main class

class OWChipDistance(OWWidget):
    settingsList = ["Metrics"]

    def __init__(self, parent=None, signalManager = None):
        OWWidget.__init__(self, parent, signalManager, 'Chip Distance') 
        
        self.inputs = [("Structured Chip Data", ChipData, self.chipdata, 1)]
        self.outputs = [("Distance Matrix", orange.SymMatrix)]

        self.Metrics = 0
        self.loadSettings()
        self.data = []
        self.metrics = [("Euclidean", orange.ExamplesDistanceConstructor_Euclidean),
            ("Manhattan", orange.ExamplesDistanceConstructor_Manhattan),
            ("Hamiltonian", orange.ExamplesDistanceConstructor_Hamiltonian)]

        # Info box
        box = QVGroupBox("Info", self.controlArea)
        self.infoa = QLabel('No data on input.', box)
        self.infob = QLabel('', box)
        OWGUI.separator(self.controlArea)

        # Distance metrics selection
        items = [x[0] for x in self.metrics]
        OWGUI.comboBox(self.controlArea, self, "Metrics", box="Distance Metrics", items=items,
            tooltip="Metrics to measure distance between data sets.",
            callback=self.computeMatrix)

        self.resize(100, 50)

    ##########################################################################
    # handling of input/output signals

    def computeDistance(self, d1, d2, dist):
        d = 0
        for i in range(len(d1)):
            d += dist(d1[i], d2[i])
        d = d / len(d1)
        return d

    def computeMatrix(self, pb):
        if not self.data:
            self.send("Distance Matrix", None)
            return
        data = self.data
        if self.Metrics == 0: # bug in orange, correct (remove normalize) once it is fixed
            dist = self.metrics[self.Metrics][1](data[0], normalize=0)
        else:
            dist = self.metrics[self.Metrics][1](data[0])            
        
        matrix = orange.SymMatrix(len(data))
        matrix.setattr('items', self.data)
        for i in range(len(data)-1):
            for j in range(i+1, len(data)):
                matrix[i, j] = self.computeDistance(data[i], data[j], dist)
                pb.advance()
        self.send("Distance Matrix", matrix)

    def chipdata(self, data):
        self.data = []
        if not data:
            return
        for (strainname, ds) in data:
            for d in ds:
                d.setattr("strain", strainname)
                self.data.append(d)
        pb = OWGUI.ProgressBar(self, iterations=len(self.data))
        self.infoa.setText('%d data file sets' % len(data))
        self.infob.setText('%d data files' % len(self.data))
        self.computeMatrix(pb)
        pb.finish()

##########################################################################
# test script

def loadData(root):
    chipdata = [] # structured [(dirname0, [d00, d01, ...]), ...]
    datasets = []  # flat list containing all the data sets
    dirs = os.listdir(root)
    for d in dirs:
        dirname = root+'\\'+d
        if os.path.isdir(dirname):
            dirdata = []   
            files  = os.listdir(dirname)
            for f in files:
                name, ext = os.path.splitext(f)
                if ext in ['.tab', '.txt', '.data']:
                    try:
                        data = None
                        data = orange.ExampleTable(dirname+'\\'+f)
                        data.name = name
                        dirdata.append(data)
                    except orange.KernelException:
                        print 'eee Exception, file not in appropriate format'
            if len(dirdata):
                chipdata.append((os.path.split(dirname)[1], dirdata))
                datasets = datasets + dirdata
    return chipdata

if __name__=="__main__":
    import orange
    a = QApplication(sys.argv)
    ow = OWChipDistance()
    a.setMainWidget(ow)
    ow.show()
#    data = loadData(r"c:\Python23\Lib\site-packages\orange\OrangeWidgets\Genomics\smallchipdata")
    data = loadData(r"c:\Python23\Lib\site-packages\orange\OrangeWidgets\Genomics\testchip")
    ow.chipdata(data)

    a.exec_loop()
    ow.saveSettings()
