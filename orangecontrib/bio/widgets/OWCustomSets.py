"""
<name>Custom Gene Sets</name>
<description>Manage custom geneset files</description>
<icon>icons/customSets.svg</icon>
<contact>Vid Jelen (vid.jelen1@gmail.com)</contact>
"""
import os
from StringIO import StringIO

import Orange
import cPickle as pickle

from OWWidget import *
import OWGUI

from Orange.bio.obiGeneSets import loadGMT, list_local, register, local_path, remove_local, modification_date, getGenesetsStats

NAME = "Custom Gene Sets"
DESCRIPTION = "Manage custom geneset files"
ICON = "icons/customSets.svg"
PRIORITY = 5000

INPUTS = []
OUTPUTS = []

REPLACES = ["_bioinformatics.widgets.OWCustomSets.OWCustomSets"]


class OWCustomSets(OWWidget):

    def __init__(self, parent=None, signalManager=None,
                 title="Custom Gene Set Manager"):
        OWWidget.__init__(self, parent, signalManager, title,
                          wantMainArea=True)

        self.inputs = []
        self.outputs = []
        self.new_geneset = set()
        self.selected_file = "" 

        self.browse_button = OWGUI.button(self.controlArea, self, 'Import Gene Sets ...', callback = self.on_open_dialog)
        self.browse_button.setIcon(self.style().standardIcon(QStyle.SP_DirOpenIcon))

        # The preview field        
        box = OWGUI.widgetBox(self.controlArea, "Imported Gene Sets")
        self.preview_view = QTreeWidget()
        self.preview_view.setAllColumnsShowFocus(1)
        self.preview_view.setColumnCount(3)
        self.preview_view.setHeaderLabels(["Name", "# of Genes", "Genes"])

        self.preview_view.header().setStretchLastSection(True)
        self.preview_view.header().setClickable(True)
        self.preview_view.header().setSortIndicatorShown(True)
        self.preview_view.setSortingEnabled(True)

        # The geneset table
        ma = self.mainArea

        self.listView = QTreeWidget(ma)

        # Adding the widgets into separate layouts
        ma.layout().addWidget(self.preview_view)
        box.layout().addWidget(self.listView)

        self.listView.setAllColumnsShowFocus(1)
        self.listView.setColumnCount(2)
        self.listView.setHeaderLabels(["Name", "Import time"])

        self.listView.header().setStretchLastSection(True)

        self.listView.setSelectionMode(QAbstractItemView.SingleSelection)
        self.listView.setSelectionBehavior(QAbstractItemView.SelectRows)

        self.populate_table()

        self.resize(800, 500)

        #Data Set info bar
        info_box = OWGUI.widgetBox(self.controlArea, "Info")
        self.info = OWGUI.widgetLabel(info_box, "No gene set selected")
        self.connect(self.listView, SIGNAL("itemSelectionChanged()"), self.selection)
        self.connect(self.listView, SIGNAL("itemSelectionChanged()"), self.update_preview)

        info_box.layout().addWidget(self.info)       
        
        OWGUI.button(self.controlArea, self, "Delete", callback=self.delete_data)

    def selection(self): 
        if self.listView.selectedItems():
            self.info.clear()
            name = self.listView.selectedItems()[0].text(0).replace(" - ", "_._") 
            for geneset in os.listdir(local_path()):
                if geneset.__contains__(str(name)):
                    the_file = os.path.join(local_path(), geneset) 
                    sets = pickle.load(open(the_file, "rb"))
                    stats = getGenesetsStats(sets)
                    num_sets, uniq_genes, avg_genes = str(stats[0]), str(stats[1]), str(stats[2])
                    break
            self.info.setText("Gene Sets: %d\nUnique Genes: %d\nAverage Gene Set Size: %d" % (int(num_sets), int(uniq_genes), int(avg_genes)))
        else:
            self.info.setText("No gene set selected")


    def populate_table(self):
        self.listView.clear()
        for geneset in os.listdir(local_path()):
            item = QTreeWidgetItem(self.listView)
            name = geneset[geneset.index("_")+1:geneset.rfind("_._")]
            the_file = os.path.join(local_path(), geneset)
            mod_time = str(modification_date(the_file))
            item.setText(0, name.replace("_._", " - "))
            item.setText(1, mod_time[:mod_time.rfind(".")])

    def on_select_recent(self, recent):
        if isinstance(recent, int):
            recent = self.recent_files[recent]

        self.set_selected_file(recent)

    def on_open_dialog(self):
        last = os.path.expanduser("~/Documents")
        path = QFileDialog.getOpenFileName(self, "Open File","" ,"geneset file *.gmt (*.gmt);;All files (*)", last)
        path = unicode(path)
        if path:
            self.set_selected_file(path)

    def set_selected_file(self, filename):
        basedir, name = os.path.split(filename)
        self.selected_file = filename
        self.genesetname = name
            
        self.import_data()
                                                                    
    def update_preview(self):
        if self.listView.selectedItems():
            final_text = ""
            self.preview_view.clear()
            name = self.listView.selectedItems()[0].text(0).replace(" - ", "_._")            
            for geneset in os.listdir(local_path()):
                if geneset.__contains__(str(name)):
                    the_file = os.path.join(local_path(), geneset) 
                    sets = pickle.load(open(the_file, "rb"))
                    break
            for geneset in sets:
                item = QTreeWidgetItem(self.preview_view)
                item.setText(0, geneset.id)
                item.setData(1, Qt.DisplayRole, len(geneset.genes))
                item.setText(2, ", ".join(list(geneset.genes)[:5]) + ", ...")
        else:
            self.preview_view.clear()

    def import_data(self):    
        self.error(0)     
        if self.selected_file:
            try:
                self.new_geneset = loadGMT(open(self.selected_file, "rt").read(), self.genesetname)
                geneset_split = self.new_geneset.split_by_hierarchy()
                for geneset in geneset_split:
                    register(geneset)
                self.populate_table()
            except Exception, ex:
                self.error(0, "An error occurred while "
                              "loading the file:\n\t%r" % self.selected_file
                              )   

    def delete_data(self):   
        self.error(0)
        if self.listView.selectedItems():
            try:
                unwantedGeneset = str(self.listView.selectedItems()[0].text(0))
                remove_local(unwantedGeneset)
                self.populate_table()
            except Exception, ex:
                self.error(0, "An error occurred while "
                              "deleting the file:\n\t%r" % unwantedGeneset
                              )    
       
if __name__ == "__main__":
    import sys
    app = QApplication(sys.argv)
    w = OWCustomSets()
    w.show()
    app.exec_()
    w.saveSettings()
