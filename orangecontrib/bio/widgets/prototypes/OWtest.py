"""
<name>test</name>
<description>Manage custom geneset files</description>

"""
import os
from StringIO import StringIO

import Orange

from OWWidget import *
import OWGUI

from Orange.bio.obiGeneSets import loadGMT, list_local, register, local_path, remove_local, modification_date

class standard_icons(object):
    def __init__(self, qwidget=None, style=None):
        self.qwidget = qwidget
        if qwidget is None:
            self.style = QApplication.instance().style()
        else:
            self.style = qwidget.style()

    @property
    def dir_open_icon(self):
        return self.style.standardIcon(QStyle.SP_DirOpenIcon)

    @property
    def reload_icon(self):
        return self.style.standardIcon(QStyle.SP_BrowserReload)


class OWtest(OWWidget):
    settingsList = ["recent_files"]

    def __init__(self, parent=None, signalManager=None,
                 title="Custom Geneset File (*.gmt) Manager"):
        OWWidget.__init__(self, parent, signalManager, title,
                          wantMainArea=True)

        self.inputs = []
        self.outputs = []
        self.new_geneset = set()

        # List of recent opened files.
        self.recent_files = []
        self.loadSettings()
        self.recent_files = filter(os.path.exists, self.recent_files)

        layout = QHBoxLayout()
        box = OWGUI.widgetBox(self.controlArea, "File", orientation=layout)

        icons = standard_icons(self)

        self.recent_combo = QComboBox(self, objectName="recent_combo",
                                      toolTip="Recent files",
                                      activated=self.on_select_recent)
        self.recent_combo.addItems([os.path.basename(p) \
                                    for p in self.recent_files])

        self.browse_button = QPushButton("...", icon=icons.dir_open_icon,
                                         toolTip="Browse filesystem",
                                         clicked=self.on_open_dialog)

        layout.addWidget(self.recent_combo, 2)
        layout.addWidget(self.browse_button)
       
        # The preview field
        form = QFormLayout()
        
        box = OWGUI.widgetBox(self.controlArea, "Preview")
        self.preview_view = QTableView()

        box.layout().addWidget(self.preview_view)

        OWGUI.button(self.controlArea, self, "Import", callback=self.import_data)

        # The geneset table
        ma = self.mainArea

        self.listView = QTreeWidget(ma)
        ma.layout().addWidget(self.listView)
        self.listView.setAllColumnsShowFocus(1)
        self.listView.setColumnCount(2)
        self.listView.setHeaderLabels(["Genesets name", "Import time"])

        self.listView.header().setStretchLastSection(True)
        self.listView.header().setClickable(True)
        self.listView.header().setSortIndicatorShown(True)
        self.listView.setSortingEnabled(True)

        self.listView.setSelectionMode(QAbstractItemView.SingleSelection)
        self.listView.setSelectionBehavior(QAbstractItemView.SelectRows)

        self.populate_table()

        OWGUI.button(self.controlArea, self, "Delete", callback=self.delete_data)

        self.selected_file = None 

        self.resize(450, 500)
        if self.recent_files:
            QTimer.singleShot(1,
                    lambda: self.set_selected_file(self.recent_files[0])
                    )

    def populate_table(self):
        self.listView.clear()
        for geneset in os.listdir(local_path()):
            item = QTreeWidgetItem(self.listView)
            name = geneset[geneset.index("gs_")+3:geneset.index(".gmt")+4]
            the_file = os.path.join(local_path(), geneset)
            mod_time = str(modification_date(the_file))
            item.setText(0, name)
            item.setText(1, mod_time[:mod_time.rfind(".")])

        print list_local()

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

    def on_reload_file(self):
        if self.recent_files:
            self.set_selected_file(self.recent_files[0])

    def set_selected_file(self, filename):
        basedir, name = os.path.split(filename)
        self.selected_file = filename
        self.genesetname = name
        index_to_remove = None
        if filename in self.recent_files:
            index_to_remove = self.recent_files.index(filename)
        elif self.recent_combo.count() > 6:
            # Always keep 6 latest files in the list.
            index_to_remove = self.recent_combo.count() - 1
        self.recent_combo.insertItem(0, name)
        self.recent_combo.setCurrentIndex(0)
        self.recent_files.insert(0, filename)

        if index_to_remove is not None:
            self.recent_combo.removeItem(index_to_remove + 1)
            self.recent_files.pop(index_to_remove + 1)
 
        self.update_preview()
    
    def update_preview(self):
        pass
        """ Leftovers from copying another widget
        
        self.error(0)
        if self.selected_file:
            head = StringIO("".join(self.selected_file_head))
            hints = self.hints[self.selected_file]

            # Save hints for the selected file
            hints["quotechar"] = self.quote
            hints["delimiter"] = self.delimiter or self.other_delimiter
            hints["has_header"] = self.has_header
            hints["has_orange_header"] = self.has_orange_header
            hints["skipinitialspace"] = self.skipinitialspace
            hints["DK"] = self.missing or None
            try:
                data = Orange.data.io.load_csv(head, delimiter=self.delimiter,
                                   quotechar=self.quote,
                                   has_header=self.has_header,
                                   has_types=self.has_orange_header,
                                   has_annotations=self.has_orange_header,
                                   skipinitialspace=self.skipinitialspace,
                                   DK=self.missing or None,
                                   create_new_on=MakeStatus.OK)
            except Exception, ex:
                self.error(0, "Cannot parse (%r)" % ex)
                data = None

            if data is not None:
                model = ExampleTableModel(data, None, self)
            else:
                model = None
            self.preview_view.setModel(model)
        """
    
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
    w = OWCustomGeneSetUpload()
    w.show()
    app.exec_()
    w.saveSettings()
