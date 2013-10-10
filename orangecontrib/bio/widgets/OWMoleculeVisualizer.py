"""<name>Molecule Visualizer</name>
<description>Rendering of 2D structure of molecules based on their SMILES description.</description>
<icon>icons/MoleculeVisualizer.svg</icon>
<contact>Ales Erjavec (ales.erjavec(@at@)fri.uni-lj.si)</contact> 
<priority>2050</priority>
"""

from __future__ import absolute_import, with_statement

import sys, os, urllib2, urllib
import warnings
import shelve
from cStringIO import StringIO
import pickle

from PyQt4.QtSvg import *

import orange
from Orange.orng import orngEnviron
from Orange.OrangeWidgets import OWGUI
from Orange.OrangeWidgets.OWWidget import *

from .. import obiChem

NAME = "Molecule Visualizer"
DESCRIPTION = "Rendering of 2D structure of molecules based on their SMILES description."
ICON = "icons/MoleculeVisualizer.svg"
PRIORITY = 2050

INPUTS = [("Molecules", Orange.data.Table, "setMoleculeTable", Default),
		  ("Molecule subset", Orange.data.Table, "setMoleculeSubset"),
		  ("Fragments", Orange.data.Table, "setFragmentTable")]
OUTPUTS = [("Selected Molecules", Orange.data.Table)]

REPLACES = ["_bioinformatics.widgets.OWMoleculeVisualizer.OWMoleculeVisualizer"]


class dummy_module(object):
	def __init__(self, name):
		self.name = name
	def __getattr__(self, name):
		raise ImportError(self.name)
	def __nonzero__(self):
		return False

try:
	import pybel
except ImportError:
	pybel = dummy_module("pybel")
except Exception, ex:
	pybel = dummy_module("pybel")

svg_error_string = """<?xml version="1.0" ?>
<svg height="185" version="1.0" width="250" xmlns="http://www.w3.org/2000/svg">
	<g stroke="#000" stroke-width="1.0">
		<text font-family="Arial" font-size="16" stroke="rgb(0, 0, 0)" x="0.0" y="16.0">
			%s
		</text>
	</g>
</svg>
"""

oasaLocal = True

#  try import from bkchem 
try:
    from bkchem import oasa
except ImportError:
    #  try stanalone import 
	try:
	    import oasa
	except ImportError:
	    oasaLocal = False

if oasaLocal and pybel:
    def mol_to_svg(molSmiles, fragSmiles):
        s = StringIO()
        obiChem.mol_to_svg(molSmiles, fragSmiles, s)
        s.seek(0)
        return s.read()
else:
    def mol_to_svg(molSmiles, fragSmiles):
        params = urllib.urlencode({'molSmiles': molSmiles, 'fragSmiles': fragSmiles})
        f = urllib2.urlopen("http://asterix.fri.uni-lj.si/misc/bkchem/drawMol_oasa.py", params)
        return f.read()
  
class DrawContext(object):
    def __init__(self, molecule="", fragment="", size=200, imageprefix="", imagename="", title="", grayedBackground=False, useCached=False):
        self.molecule=molecule
        self.fragment=fragment
        self.size=size
        self.imageprefix=imageprefix
        self.imagename=imagename
        self.title=title
        self.grayedBackground=grayedBackground
        self.useCached=useCached

    def __hash__(self):
        return (self.molecule+"%%"+self.fragment+"%%"+str(self.size)+"%%"+self.title+str(self.grayedBackground)).__hash__()

    def __eq__(self, other):
        return self.molecule==other.molecule and self.fragment==other.fragment and \
               self.size==other.size and self.title==other.title and self.grayedBackground==other.grayedBackground

    def __ne__(self, other):
        return not self==other
            
from threading import RLock

def synchronized(lock):
    def syncfunc(func):
        def f(*args, **kw):
            with lock:
                ret = func(*args, **kw)
            return ret
        return f
    return syncfunc

class ImageCache(object):
    __shared_state = {"shelve": None}
    lock = RLock()
    @synchronized(lock)
    def __init__(self):
        self.__dict__ = self.__shared_state
        if self.shelve is None:
            try:
                os.mkdir(os.path.join(orngEnviron.bufferDir, "molimages"))
            except OSError:
                pass
            try:
                self.shelve = shelve.open(os.path.join(orngEnviron.bufferDir, "molimages", "cache.shelve"))
            except Exception, ex:
                warnings.warn("Cannot open molecule images cache! " + str(ex))
                self.shelve = {}
    
    
    @synchronized(lock)
    def __getitem__(self, key):
        val_str = str(key)
        if val_str in self.shelve:
            return self.shelve[val_str]
        else:
            mol_smiles = key[0]
            frag_smiles = key[1] if len(key) > 1 else None
            self.shelve[val_str] = mol_to_svg(mol_smiles, frag_smiles)
            return self.shelve[val_str]

    @synchronized(lock)
    def __setitem__(self, key, value):
        if len(self.shelve) > 1000:
            self.sync()
        self.shelve[str(key)] = value

            
    @synchronized(lock)
    def __contains__(self, key):
        return str(key) in self.shelve
    
    @synchronized(lock)
    def sync(self):
        if len(self.shelve.keys()) > 1000:
            for key in self.shelve.keys()[:-900]:
                del self.shelve[key]
        if hasattr(self.shelve, "sync"):
            self.shelve.sync()

    def __del__(self):
    	if hasattr(self.shelve, " sync"):
            self.sync()

class MolWidget(QFrame):
    def setSelected(self, val):
        self.image.selected = val
    def getSelected(self):
        return self.image.selected
    selected = property(getSelected, setSelected)

    def __init__(self, master, parent, context):
        QFrame.__init__(self, parent)
        self.master=master
        self.context=context
        self.label=QLabel()
        self.label.setSizePolicy(QSizePolicy.Minimum, QSizePolicy.Fixed)
        try:
        	self.from_cache = (context.molecule, context.fragment) in ImageCache()
        	s = ImageCache()[context.molecule, context.fragment]
        	self.state = 1
        except Exception, ex:
            from traceback import print_exc
            print_exc()
            s = svg_error_string % "Error loading: "+str(ex)
            self.from_cache = False
            self.state = 0
        self.image=SVGImageThumbnail(s, self)
        self.image.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
        self.label.setText(context.title)
        self.label.setAlignment(Qt.AlignHCenter | Qt.AlignVCenter)
        self.label.setMaximumWidth(context.size)
        layout = QVBoxLayout()
        layout.addWidget(self.label)
        layout.addWidget(self.image)
        self.setLayout(layout)
        self.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
        self.show()

    def repaint(self):
        QFrame.repaint(self)
        self.label.repaint()
        self.image.repaint()
        
class SVGImageThumbnail(QFrame):
    def setSelected(self, val):
        self._selected = val
        self.update()
    def getSelected(self):
        return self._selected
    selected = property(getSelected, setSelected)
    
    @property
    def highlight(self):
        return self.parent().context.grayedBackground
    
    def __init__(self, file, parent):
        QWidget.__init__(self, parent)
        self.doc = file if type(file) == str else file.read()
        self.renderer = QSvgRenderer(QByteArray(self.doc))
        self.buffer = None
        self.resize(200, 200)
        self.selected = False
        
    def paintEvent(self, event):
        if not self.buffer or self.buffer.size() != self.size():
            defSize = self.renderer.defaultSize()
            scaleFactor = float(self.width()) / max(defSize.width(), defSize.height())
            self.buffer = QImage(self.size(), QImage.Format_ARGB32_Premultiplied)
            self.buffer.fill(0)
            painter = QPainter(self.buffer)
            painter.setViewport(self.width()/2 - scaleFactor*defSize.width()/2, self.height()/2 - scaleFactor*defSize.height()/2,
                                scaleFactor*defSize.width(), scaleFactor*defSize.height())
            self.renderer.render(painter)
        painter = QPainter(self)            
        painter.setBrush(QBrush(Qt.white))
        painter.drawRect(0, 0, self.width(), self.height())
        if self.selected:
            painter.setPen(QPen(QBrush(Qt.red), 2))
        if self.highlight:
            painter.setBrush(QBrush(Qt.gray,  Qt.FDiagPattern))
        else:
            painter.setBrush(Qt.NoBrush)
        painter.drawRect(0, 0, self.width(), self.height())
        painter.drawImage(0, 0, self.buffer)

    def sizeHint(self):
        return QSize(self.parent().master.imageSize, self.parent().master.imageSize)

    def mouseDoubleClickEvent(self, event):
        self._bigimage = BigSvgWidget()
        self._bigimage.load(QByteArray(self.doc))
        self._bigimage.show()

    def mousePressEvent(self, event):
        self.parent().master.mouseAction(self.parent(), event)

class BigSvgWidget(QSvgWidget):
    def paintEvent(self, event):
        painter = QPainter(self)
        painter.setBrush(QBrush(Qt.white))
        painter.drawRect(0, 0, self.width(), self.height())
        painter.end()
        QSvgWidget.paintEvent(self, event)
        
class ScrollArea(QScrollArea):
    def __init__(self, master, *args):
        QScrollArea.__init__(self, *args)
        self.master=master
        self.viewport().setMouseTracking(True)
        self.setMouseTracking(True)
        
    def resizeEvent(self, event):
        QScrollArea.resizeEvent(self, event)
        size = event.size()
        w, h=self.width(), self.height()
        oldNumColumns = self.master.numColumns
        numColumns = w / (self.master.imageSize + self.master.gridLayout.horizontalSpacing()*2 + 20) or 1
        if numColumns != oldNumColumns:
            self.master.numColumns = numColumns
            self.master.rearrangeLayout()

class OWMoleculeVisualizer(OWWidget):
    settingsList=["colorFragmets","showFragments"]
    contextHandlers={"":DomainContextHandler("", [ContextField("moleculeTitleAttributeList",
                                                    DomainContextHandler.List + DomainContextHandler.SelectedRequired + DomainContextHandler.IncludeMetaAttributes,
                                                    selected="selectedMoleculeTitleAttrs"),
                                                  ContextField("moleculeSmilesAttr", DomainContextHandler.Required + DomainContextHandler.IncludeMetaAttributes)], maxAttributesToPickle=10000)} ##maxAttributesToPickle=10000 some bug in Context handler 
    def __init__(self, parent=None, signalManager=None, name="Molecule visualizer"):
        OWWidget.__init__(self, parent, signalManager, name)
        self.inputs=[("Molecules", ExampleTable, self.setMoleculeTable), ("Molecule subset", ExampleTable, self.setMoleculeSubset), ("Fragments", ExampleTable, self.setFragmentTable)]
        self.outputs=[("Selected Molecules", ExampleTable)]
        self.colorFragments=1
        self.showFragments=0
        self.selectedFragment=""
        self.moleculeSmiles=[]
        self.fragmentSmiles=[]
        self.defFragmentSmiles=[]
        self.moleculeSmilesAttr=0
        self.moleculeTitleAttr=0
        self.moleculeTitleAttributeList=[]
        self.selectedMoleculeTitleAttrs=[]
        self.fragmentSmilesAttr=0
        self.imageSize=200
        self.numColumns=4
        self.commitOnChange=0
        self.overRideCache=True
        ##GUI
        box=OWGUI.widgetBox(self.controlArea, "Info", addSpace = True)
        self.infoLabel=OWGUI.label(box, self, "Chemicals: ")
##        if not oasaLocal:
##            OWGUI.label(box, self, "OpenEye not installed, access to server required.")
##            self.serverInfo = OWGUI.label(box, self, "")
        box=OWGUI.radioButtonsInBox(self.controlArea, self, "showFragments", ["Show molecules", "Show fragments"], "Show", callback=self.showImages)
        self.showFragmentsRadioButton=box.buttons[-1]
        self.markFragmentsCheckBox=OWGUI.checkBox(box, self, "colorFragments", "Mark fragments", callback=self.redrawImages)
        box.setSizePolicy(QSizePolicy(QSizePolicy.Minimum, QSizePolicy.Maximum))
        OWGUI.separator(self.controlArea)
        self.moleculeSmilesCombo=OWGUI.comboBox(self.controlArea, self, "moleculeSmilesAttr", "Molecule SMILES Attribute",callback=self.showImages)
        self.moleculeSmilesCombo.box.setSizePolicy(QSizePolicy(QSizePolicy.Minimum, QSizePolicy.Maximum))
        OWGUI.separator(self.controlArea)
##        self.moleculeTitleCombo=OWGUI.comboBox(self.controlArea, self, "moleculeTitleAttr", "Molecule title attribute", callback=self.redrawImages)
        box=OWGUI.widgetBox(self.controlArea, "Molecule Title Attributes", addSpace = True)
##        self.moleculeTitleListBox=QListBox(box)
##        self.moleculeTitleListBox.setSelectionMode(QListBox.Extended)
##        self.moleculeTitleListBox.setMinimumHeight(100)
##        self.connect(self.moleculeTitleListBox, SIGNAL("selectionChanged()"), self.updateTitles)
        self.moleculeTitleListBox=OWGUI.listBox(box, self, "selectedMoleculeTitleAttrs", "moleculeTitleAttributeList", selectionMode = QListWidget.ExtendedSelection, callback=self.updateTitles)
        self.moleculeTitleListBox.setMinimumHeight(100)
##        OWGUI.separator(self.controlArea)
##        self.moleculeTitleCombo.box.setSizePolicy(QSizePolicy(QSizePolicy.Minimum, QSizePolicy.Maximum))
        OWGUI.separator(self.controlArea)
        self.fragmentSmilesCombo=OWGUI.comboBox(self.controlArea, self, "fragmentSmilesAttr", "Fragment SMILES Attribute", callback=self.updateFragmentsListBox)
        self.fragmentSmilesCombo.box.setSizePolicy(QSizePolicy(QSizePolicy.Minimum, QSizePolicy.Maximum))
        OWGUI.separator(self.controlArea)
        box=OWGUI.spin(self.controlArea, self, "imageSize", 50, 500, 50, box="Image Size", callback=self.redrawImages, callbackOnReturn = True)
        box.setSizePolicy(QSizePolicy(QSizePolicy.Minimum, QSizePolicy.Maximum))
        OWGUI.separator(self.controlArea)
        box=OWGUI.widgetBox(self.controlArea,"Selection", addSpace = True)
        OWGUI.checkBox(box, self, "commitOnChange", "Commit on change")
        self.selectMarkedMoleculesButton=OWGUI.button(box, self, "Select &matched molecules", self.selectMarked)
        OWGUI.button(box, self, "&Commit", callback=self.commit)
        OWGUI.separator(self.controlArea)
        OWGUI.button(self.controlArea, self, "&Save to HTML", self.saveToHTML, debuggingEnabled = 0)
        OWGUI.rubber(self.controlArea)
        
        spliter=QSplitter(Qt.Vertical)
        self.scrollArea=ScrollArea(self, spliter)
        
        self.molWidget=QWidget()
        self.scrollArea.setWidget(self.molWidget)
        self.mainArea.layout().addWidget(spliter)
        self.gridLayout=QGridLayout(self.molWidget)
        self.molWidget.setLayout(self.gridLayout)

        if pybel:
        	self.listBox=QListWidget(spliter)
        else:
        	self.listBox=QListWidget(None)
        	self.listBox.setHidden(True)

        self.connect(self.listBox, SIGNAL("itemClicked(QListWidgetItem *)"), self.fragmentSelection)
        
        self.fragmentSmilesCombo.box.setDisabled(not pybel)

        self.imageWidgets=[]
        self.candidateMolSmilesAttr=[]
        self.candidateMolTitleAttr=[None]
        self.candidateFragSmilesAttr=[None]
        self.molData=None
        self.molSubset=[]
        self.fragData=None
        self.ctrlPressed=FALSE
        self.resize(800,600)
        self.listBox.setMaximumHeight(150)
        self.fragmentSmilesCombo.setDisabled(True)
        self.selectMarkedMoleculesButton.setDisabled(True)
        self.markFragmentsCheckBox.setDisabled(True)
        self.showFragmentsRadioButton.setDisabled(True)
        self.loadSettings()
        self.failedCount=0
        self.fromCacheCount=0
        
        if not pybel:
        	self.showFragments = 0
        	self.warning(10, "Pybel module not installed. To view molecule fragments\nplease install openbabel python extension.")
        if not oasaLocal:
        	self.warning(11, "OASA module not installed. For faster local molecule\nrendering install OASA module (part of BKChem)")
        
    def setMoleculeTable(self, data):
        self.closeContext()
        self.molData=data
        if data:
            self.setMoleculeSmilesCombo()
            self.setMoleculeTitleListBox()
            self.setFragmentSmilesCombo()
            self.updateFragmentsListBox()
            if self.molSubset:
                try:
                    self.molSubset=self.molSubset.select(self.molData.domain)
                except:
                    self.molSubset=[]
            tmp=self.moleculeTitleAttributeList
            self.openContext("", data)
            if tmp and not self.moleculeTitleAttributeList: ##openContext somtimes crashes internaly and silently clears title list
                self.moleculeTitleAttributeList=tmp
            self.showImages()
        else:
            self.moleculeSmilesCombo.clear()
            self.moleculeTitleListBox.clear()
            self.moleculeSmilesAttr=0
            self.moleculeTitleAttributeList=[]
            self.selectedMoleculeTitleAttrs=[]
            self.defFragmentSmiles=[]
            if not self.fragmentSmilesAttr:
                self.listBox.clear()
            self.destroyImageWidgets()
            self.openContext("",data)
            self.send("Selected Molecules", None)
    
    def setMoleculeSubset(self, data):
        self.molSubset=data
        try:
            self.molSubset=self.molSubset.select(self.molData.domain)
        except:
            self.molSubset=[]
        self.showImages()
            
    def setFragmentTable(self, data):
        self.fragData=data
        if data:
            self.setFragmentSmilesCombo()
            self.updateFragmentsListBox()
            self.selectedFragment=""
            self.showImages()
        else:
            self.setFragmentSmilesCombo()
            #self.fragmentSmilesAttr=0
            self.updateFragmentsListBox()
            if self.showFragments:
                self.destroyImageWidgets()
        self.fragmentSmilesCombo.setDisabled(bool(data))

    def filterSmilesVariables(self, data):
#    	import pybel
        candidates = data.domain.variables+data.domain.getmetas().values()
        candidates = filter(lambda v:v.varType==orange.VarTypes.Discrete or v.varType==orange.VarTypes.String, candidates)
        if len(data)>20:
            data=data.select(orange.MakeRandomIndices2(data, 20))
        vars=[]
        def isValidSmiles(s):
            try:
                pybel.readstring("smi", s)
            except IOError:
                return False
            except ImportError:
            	return True
            return True
            
        import os
        tmpFd1=os.dup(1)
        tmpFd2=os.dup(2)
        fd=os.open(os.devnull, os.O_APPEND)
##        os.close(1)
        os.dup2(fd, 1)
        os.dup2(fd, 2)
##        os.close(fd)
        for var in candidates:
            count=0
            for e in data:
                if pybel and isValidSmiles(str(e[var])):
                    count+=1
            vars.append((count,var))
        names=[v.name for v in data.domain.variables+data.domain.getmetas().values()]
        names=filter(isValidSmiles, names)
##        os.close(1)
        os.dup2(tmpFd1, 1)
        os.dup2(tmpFd2, 2)
##        os.close(tmpFd)
        return vars, names
        
    def setMoleculeSmilesCombo(self):
        candidates, self.defFragmentSmiles=self.filterSmilesVariables(self.molData)
        self.candidateMolSmilesAttr=[c[1] for c in candidates]
        best = reduce(lambda best,current:best[0]<current[0] and current or best, candidates)
        self.moleculeSmilesCombo.clear()
        self.moleculeSmilesCombo.addItems([v.name for v in self.candidateMolSmilesAttr])
        self.moleculeSmilesAttr=candidates.index(best)

    def setMoleculeTitleListBox(self):
        self.icons=self.createAttributeIconDict()
        vars=self.molData.domain.variables+self.molData.domain.getmetas().values()
        self.moleculeTitleAttributeList=[attr.name for attr in vars]
        self.selectedMoleculeTitleAttrs=[]

    def updateTitles(self):
        if not self.molData:
            return
        smilesAttr=self.candidateMolSmilesAttr[min(self.moleculeSmilesAttr, len(self.candidateMolSmilesAttr)-1)]
        attrs = self.molData.domain.variables+self.molData.domain.getmetas().values()
        selected=[attrs[i] for i in self.selectedMoleculeTitleAttrs]
        for widget, example in zip(self.imageWidgets, filter(lambda e:not e[smilesAttr].isSpecial(),self.molData)):
            text=" / ".join(map(str, [example[attr] for attr in selected]))
            widget.label.setText(text)

    def setFragmentSmilesCombo(self):
        if self.fragData:
            candidates, names=self.filterSmilesVariables(self.fragData)
        else:
            candidates=[]
        self.candidateFragSmilesAttr=[None]+candidates
        self.fragmentSmilesCombo.clear()
        self.fragmentSmilesCombo.addItems(["Default"]+[v.name for v in candidates])
        if self.fragmentSmilesAttr>len(candidates):
            self.fragmentSmilesAttr=0

    def updateFragmentsListBox(self):
    	if not pybel:
    		return
        fAttr=self.candidateFragSmilesAttr[self.fragmentSmilesAttr]
        if fAttr:
            self.fragmentSmiles=[""]+[str(e[fAttr]) for e in self.fragData if not e[fAttr].isSpecial()]
        else:
            self.fragmentSmiles=[""]+self.defFragmentSmiles
        self.listBox.clear()
        self.listBox.addItems(self.fragmentSmiles)
        self.showFragmentsRadioButton.setDisabled(len(self.fragmentSmiles)==1)
        self.markFragmentsCheckBox.setDisabled(len(self.fragmentSmiles)==1)
        self.selectMarkedMoleculesButton.setDisabled(True)
        
    def fragmentSelection(self, item):
    	if not pybel:
    		return
        index = self.listBox.indexFromItem(item).row()
        if index == -1:
            index = 0
        self.selectedFragment=self.fragmentSmiles[index]
        self.selectMarkedMoleculesButton.setEnabled(bool(self.selectedFragment))
        self.markFragmentsCheckBox.setEnabled(bool(self.selectedFragment))
        if not self.showFragments and self.colorFragments:
            self.redrawImages()
        
    def renderImages(self,useCached=False):
##        def fixNumColumns(numItems, numColumns):
##            if (self.imageSize+4)*(numItems/numColumns+1)>30000:
##                return numItems/(30000/(self.imageSize+4))
##            else:
##                return numColumns
        
        self.numColumns=self.scrollArea.width()/(self.imageSize+4) or 1
        self.gridLayout = QGridLayout()
        self.scrollArea.takeWidget()
        self.molWidget = QWidget()
        self.scrollArea.setWidget(self.molWidget)
        self.molWidget.setLayout(self.gridLayout)
        self.imageWidgets=[]
##        self.imageCache.newEpoch()
        self.failedCount=0
        self.fromCacheCount=0
        if self.showFragments and self.fragmentSmiles:
            correctedNumColumns=self.numColumns #fixNumColumns(len(self.fragmentSmiles[1:]), self.numColumns)
            self.progressBarInit()
            for i,fragment in enumerate(self.fragmentSmiles[1:]):
                #imagename=self.imageprefix+str(i)+".bmp"
                #vis.molecule2BMP(fragment, imagename, self.imageSize)
##                image=MolImage(self,  self.molWidget, DrawContext(molecule=fragment, imagename=imagename, size=self.imageSize))
                image=MolWidget(self, self.molWidget, DrawContext(molecule=fragment, size=self.imageSize))
                self.gridLayout.addWidget(image, i/correctedNumColumns, i%correctedNumColumns)
                self.imageWidgets.append(image)
                self.progressBarSet(i*100/len(self.fragmentSmiles))
            self.progressBarFinished()
        elif self.molData and self.candidateMolSmilesAttr:
            sAttr=self.candidateMolSmilesAttr[min(self.moleculeSmilesAttr, len(self.candidateMolSmilesAttr)-1)]
            tAttr=self.candidateMolTitleAttr[min(self.moleculeTitleAttr, len(self.candidateMolTitleAttr)-1)]
            if self.moleculeTitleAttr:
                titleList=[str(e[tAttr]) for e in self.molData if not e[sAttr].isSpecial()]
            else:
                titleList=[]
                if not sAttr:
                    return
            molSmiles=[(str(e[sAttr]), e) for e in self.molData if not e[sAttr].isSpecial()]
            correctedNumColumns=self.numColumns #fixNumColumns(len(molSmiles), self.numColumns)
            self.progressBarInit()
            if self.colorFragments and self.selectedFragment:
                fMap=map_fragments([self.selectedFragment], [t[0] for t  in molSmiles])
            for i,((molecule, example), title) in enumerate(zip(molSmiles, titleList or [""]*len(molSmiles))):
                #imagename=self.imageprefix+str(i)+".bmp"
                if self.colorFragments and self.selectedFragment and fMap[molecule][self.selectedFragment]:
                    context=DrawContext(molecule=molecule, fragment=self.selectedFragment, size=self.imageSize, title=title, grayedBackground=example in self.molSubset, useCached=useCached)
                    #vis.moleculeFragment2BMP(molecule, self.selectedFragment, imagename, self.imageSize)
                else:
                    context=DrawContext(molecule=molecule, size=self.imageSize, title=title, grayedBackground=example in self.molSubset, useCached=useCached)
                    #vis.molecule2BMP(molecule, imagename, self.imageSize)
##                image=MolImage(self, self.molWidget, context)
                image=MolWidget(self, self.molWidget, context)
                self.gridLayout.addWidget(image, i/correctedNumColumns, i%correctedNumColumns)
                self.imageWidgets.append(image)
                self.progressBarSet(i*100/len(molSmiles))
            self.progressBarFinished()
            self.updateTitles()
        #print "done drawing"
        self.overRideCache=False
        self.fromCacheCount = sum(int(w.from_cache) for w in self.imageWidgets)
        self.failedCount = len(self.imageWidgets) - sum(int(w.state) for w in self.imageWidgets)
        
        self.molWidget.setMinimumSize(self.gridLayout.sizeHint())
        self.molWidget.show()

    def destroyImageWidgets(self):
        for w in self.imageWidgets:
            w.hide()
            self.gridLayout.removeWidget(w)
        self.imageWidgets=[]

    def rearrangeLayout(self):
        self.numColumns=self.scrollArea.width() / (self.imageSize + self.gridLayout.horizontalSpacing()*2 + 20) or 1
        self.molWidget = QWidget()
        self.gridLayout = QGridLayout()
        self.molWidget.setLayout(self.gridLayout)
        for i, w in enumerate(self.imageWidgets):
            self.gridLayout.addWidget(w, i/self.numColumns, i%self.numColumns)
        self.scrollArea.takeWidget()
        self.scrollArea.setWidget(self.molWidget)
        self.molWidget.setMinimumSize(self.gridLayout.sizeHint())
        self.molWidget.show()
            
    def showImages(self, useCached=False):
        self.destroyImageWidgets()
        self.warning(0)
        
        self.renderImages(useCached)
        if not oasaLocal:
            if self.failedCount>0:
                self.infoLabel.setText("%i chemicals\nFailed to retrieve %i images from server\n%i images from server\n%i images from local cache" \
									   % (len(self.imageWidgets), self.failedCount, len(self.imageWidgets)-self.failedCount-self.fromCacheCount, self.fromCacheCount))
                self.warning(0, "Failed to retrieve some images from server")
            elif self.fromCacheCount == len(self.imageWidgets):
                self.infoLabel.setText("%i chemicals\nAll images from local cache" % len(self.imageWidgets))
            elif self.fromCacheCount == 0:
                self.infoLabel.setText("%i chemicals\nAll images from server" % len(self.imageWidgets))
            else:
                self.infoLabel.setText("%i chemicals\n%i images from server\n%i from local cache" % (len(self.imageWidgets), len(self.imageWidgets)-self.fromCacheCount, self.fromCacheCount))
        else:
            self.infoLabel.setText("Chemicals %i" % len(self.imageWidgets))

    def redrawImages(self, useCached=False):
        selected=map(lambda i:self.imageWidgets.index(i), filter(lambda i:i.selected, self.imageWidgets))
        self.showImages(useCached=useCached)
        for i in selected:
            self.imageWidgets[i].selected=True
            self.imageWidgets[i].repaint()
            
    def mouseAction(self, image, event):
        if self.ctrlPressed:
            image.selected=not image.selected
        else:
            for i in self.imageWidgets:
                i.selected=False
                i.repaint()
            image.selected=True
        image.repaint()
        if self.commitOnChange:
            self.commit()

    def selectMarked(self):
    	if not pybel:
    		return
        if not self.showFragments:
            molecules=[i.context.molecule for i in self.imageWidgets]
            fMap=map_fragments([self.selectedFragment], molecules)
            for i in self.imageWidgets:
                if fMap[i.context.molecule][self.selectedFragment]:
                    i.selected=True
                else:
                    i.selected=False
                i.repaint()
        if self.commitOnChange:
            self.commit()
    
    def commit(self):
        if self.showFragments:
            sAttr=self.candidateMolSmilesAttr[self.moleculeSmilesAttr]
            molecules=[str(e[sAttr]) for e in self.molData]
            fragments=[i.context.molecule for i in self.imageWidgets if i.selected]
            fragmap=map_fragments(fragments, molecules)
            match=filter(lambda m:max(fragmap[m].values()), molecules)
            examples=[e for e in self.molData if str(e[sAttr]) in match]
            if examples:
                table=orange.ExampleTable(examples)
                self.send("Selected Molecules", table)
            else:
                self.send("Selected Molecules", None)                
        else:
            mols=[i.context.molecule for i in self.imageWidgets if i.selected]
            sAttr=self.candidateMolSmilesAttr[self.moleculeSmilesAttr]
            examples=[e for e in self.molData if str(e[sAttr]) in mols]
            if examples:
                table=orange.ExampleTable(examples)
                self.send("Selected Molecules", table)
            else:
                self.send("Selected Molecules", None)

    def keyPressEvent(self, key):
        if key.key()==Qt.Key_Control:
            self.ctrlPressed=TRUE
        else:
            OWWidget.keyPressEvent(self, key)

    def keyReleaseEvent(self, key):
        if key.key()==Qt.Key_Control:
            self.ctrlPressed=FALSE
        else:
            OWWidget.keyReleaseEvent(self, key)

    def saveToHTML(self):
        fileName=str(QFileDialog.getSaveFileName("index.html","HTML (.html)", None, "Save to.."))
        if not fileName:
            return
        else:
            file=open(fileName, "w")
        import os
        path, _ =os.path.split(fileName)
        if "molimages" not in os.listdir(path):
            os.mkdir(path+"/molimages")
        title="Molekule"
        file.write("<html><title>"+title+"</title>\n")
        file.write("<body> <table border=\"1\">\n")
        i=0
        try:
            import Image
        except:
            pass
        for row in range(len(self.imageWidgets)/self.numColumns+1):
            file.write("<tr>\n")
            for col in range(self.numColumns):
                if i>=len(self.imageWidgets):
                    break
                try:
                    im=Image.open(self.imageWidgets[i].context.imagename)
                    if im.mode!="RGB":
                        im=im.convert("RGB")
                    im.save(path+"/molimages/image"+str(i)+".gif", "GIF")
                    file.write("<td><p>"+str(self.imageWidgets[i].label.text())+"</p><img src=\"./molimages/image"+str(i)+".gif\"></td>\n")
                except:
                    from shutil import copy
                    copy(self.imageWidgets[i].context.imagename, path+"/molimages/")
                    file.write("<td><p>"+str(self.imageWidgets[i].label.text())+"</p><img src=\"./molimages/image"+str(i)+".bmp\"></td>\n")
                i+=1
            file.write("</tr>\n")
        file.write("</table></body></html>")
        file.close()

    def saveSettings(self, *args, **kw):
        OWWidget.saveSettings(self, *args, **kw)
        ImageCache().sync()

def molecule2BMP(molSmiles, filename, size=200, title="", grayedBackground=False):
    moleculeFragment2BMP(molSmiles, None, filename)

def moleculeFragment2BMP(molSmiles, fragSmiles, filename, size=200, title="", grayedBackground=False):
    obiChem.mol_to_svg(molSmiles, fragSmiles, filename)

def remoteMoleculeFragment2BMP(molSmiles, fragSmiles, filename, size=200, title="", grayedBackground=False):
    import cStringIO
    import urllib
    import Image, ImageDraw

    from openbabel import OBConversion, OBMol
    loader = OBConversion()
    loader.SetInAndOutFormats("smi", "smi")
    if not loader.ReadString(OBMol(), molSmiles):
        img = Image.new("RGB", (size, size), (255,255,255))
        draw = ImageDraw.Draw(img)
        draw.text((10, size/2-10), "wrong SMILES notation", fill=(0, 0, 0))
        img.save(filename)
        return
    
    params = urllib.urlencode({'password': OWMoleculeVisualizer.serverPwd, 'molSmiles': molSmiles, 'fragSmiles': fragSmiles, 'size': size, 'title': title, 'grayBack': grayedBackground})
    try:
        f = urllib.urlopen("http://212.235.189.53/openEye/drawMol.py", params)
        imgstr = f.read()
    except IOError, er:
        #print er
        raise er
    im = cStringIO.StringIO(imgstr)
    img = Image.open(im)
    #print img.format, img.size, img.mode
    img.save(filename)
    #del img

def remoteMolecule2BMP(molSmiles, filename, size=200, title="", grayedBackground=False):
    import cStringIO
    import urllib
    import Image, ImageDraw

    from openbabel import OBConversion, OBMol
    loader = OBConversion()
    loader.SetInAndOutFormats("smi", "smi")
    if not loader.ReadString(OBMol(), molSmiles):
        img = Image.new("RGB", (size, size), (255,255,255))
        draw = ImageDraw.Draw(img)
        draw.text((10, size/2-10), "wrong SMILES notation", fill=(0, 0, 0))
        img.save(filename)
        return    
    
    params = urllib.urlencode({'password': OWMoleculeVisualizer.serverPwd, 'molSmiles': molSmiles, 'size': size, 'title': title, 'grayBack': grayedBackground})
    try:
        f = urllib.urlopen("http://212.235.189.53/openEye/drawMol.py", params)
        imgstr = f.read()
    except IOError, er:
        #print er, "rm2BMP"
        raise er
    im = cStringIO.StringIO(imgstr)
    img = Image.open(im)
    #print img.format, img.size, img.mode
    img.save(filename)
    #del img        


def map_fragments(fragments, smiles, binary=True):
    ret = {}
    for s in smiles:
        mol = pybel.readstring("smi", s)
        d={}
        for f in fragments:
            pat = pybel.Smarts(f)
            count=0
            if pat.findall(mol):
                count+=1
            if binary:
                d[f]=count!=0 and 1 or 0
            else:
                d[f]=count
        ret[s]=d
    return ret

if not oasaLocal:
    moleculeFragment2BMP = remoteMoleculeFragment2BMP
    molecule2BMP = remoteMolecule2BMP

if __name__=="__main__":
    pass
##    app=QApplication(sys.argv)
##    from pywin.debugger import set_trace
####    set_trace()
##    w=OWMoleculeVisualizer()
####    app.setMainWidget(w)
##    w.show()
##    data=orange.ExampleTable("E://chem/chemdata/BCMData_growth_frag.tab")
##    
##    w.setMoleculeTable(data)
####    data=orange.ExampleTable("E://chem//new//sf.tab")
####    w.setFragmentTable(data)
##    app.exec_()
##    w.saveSettings()
