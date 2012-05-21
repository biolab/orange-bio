"""
<name>Gene Info</name>
<description>Displays gene information from NCBI and other sources.</description>
<priority>2010</priority>
<contact>Ales Erjavec (ales.erjavec(@at@)fri.uni-lj.si)</contact>
<icon>icons/GeneInfo.png</icon>
"""

from __future__ import absolute_import, with_statement

from collections import defaultdict
from functools import partial

import orange
from Orange.orng import orngServerFiles
from Orange.orng.orngDataCaching import data_hints
from Orange.OrangeWidgets import OWGUI
from Orange.OrangeWidgets.OWWidget import *

from .. import obiGene, obiTaxonomy

class TreeModel(QAbstractItemModel):
    def __init__(self, data, header, parent):
        QAbstractItemModel.__init__(self, parent)
        self._data = [[QVariant(s) for s in row] for row in data]
        self._dataDict = {}
        self._header = header
        self._roleData = {Qt.DisplayRole:self._data}
        self._roleData = partial(defaultdict, partial(defaultdict, partial(defaultdict, QVariant)))(self._roleData)
    
    def setColumnLinks(self, column, links):
        font =QFont()
        font.setUnderline(True)
        font = QVariant(font)
        for i, link in enumerate(links):
            self._roleData[LinkRole][i][column] = QVariant(link)
            self._roleData[Qt.FontRole][i][column] = font
            self._roleData[Qt.ForegroundRole][i][column] = QVariant(QColor(Qt.blue))
    
    def setRoleData(self, role, row, col, data):
        self._roleData[role][row][col] = data
        
    def data(self, index, role=Qt.DisplayRole):
        row, col = index.row(), index.column()
        return self._roleData[role][row][col]
        
    def index(self, row, col, parent=QModelIndex()):
        return self.createIndex(row, col, 0)
    
    def parent(self, index):
        return QModelIndex()
    
    def rowCount(self, index=QModelIndex()):
        if index.isValid():
            return 0
        else:
            return len(self._data)
        
    def columnCount(self, index=QModelIndex()):
        return len(self._header)

    def headerData(self, section, orientation, role=Qt.DisplayRole):
        if role==Qt.DisplayRole:
            return QVariant(self._header[section])
        return QVariant()

from Orange.OrangeWidgets.OWGUI import LinkStyledItemDelegate, LinkRole

def lru_cache(maxsize=100):
    """ A least recently used cache function decorator.
    """
    
    def decorating_function(func):
        import functools
        cache = {}
        
        functools.wraps(func)
        def wrapped(*args, **kwargs):
            key = args + tuple(sorted(kwargs.items()))
            if key not in cache:
                res = func(*args, **kwargs)
                cache[key] = (time.time(), res)
                if len(cache) > maxsize:
                    key, (_, _) = min(cache.iteritems(), key=lambda item: item[1][0])
                    del cache[key]
            else:
                _, res = cache[key]
                cache[key] = (time.time(), res) # update the time
                
            return res
        
        def clear():
            cache.clear()
        
        wrapped.clear = clear
        
        return wrapped
    return decorating_function
                
class LinkFmt(object):
    def __init__(self, link_fmt, name):
        self.link_fmt = link_fmt
        self.name = name
        
    def format(self, *args, **kwargs):
        return Link(self.link_fmt.format(*args, **kwargs), **kwargs)
    
    def __repr__(self):
        return "<LinkFmt " + repr(self.name) + " >"
    
    def __str__(self):
        return self.name
    
class Link(object):
    def __init__(self, link, text=None, **kwargs):
        self.link = link
        self.text = text if text is not None else "link"
        self.__dict__.update(kwargs)
        
    def str(self):
        return link
    
    
@lru_cache(maxsize=2)
def get_ncbi_info(taxid):
    return obiGene.NCBIGeneInfo(taxid)

def ncbi_info(taxid, genes):
    info = get_ncbi_info(taxid)
    schema_link = LinkFmt("http://www.ncbi.nlm.nih.gov/sites/entrez?Db=gene&Cmd=ShowDetailView&TermToSearch={gene_id}", name="NCBI ID")
    schema = [schema_link, "Symbol", "Locus Tag", "Chromosome",
              "Description", "Synonyms", "Nomenclature"]
    ret = []
    for gene in genes:
        gi = info.get_info(gene)
        if gi:
            ret.append([schema_link.format(gene_id=gi.gene_id, text=gi.gene_id),
                        gi.symbol + " (%s)" % gene if gene != gi.symbol else gi.symbol,
                        gi.locus_tag or "",
                        gi.chromosome or "",
                        gi.description or "",
                        ", ".join(gi.synonyms),
                        gi.symbol_from_nomenclature_authority or ""
                        ])
        else:
            ret.append(None)
    return schema, ret
    
def dicty_info(taxid, genes):
    from .. import obiDicty 
    info = obiDicty.DictyBase()
    name_matcher = obiGene.GMDicty()
    name_matcher.set_targets(info.info.keys())
    schema_link = LinkFmt("http://dictybase.org/db/cgi-bin/gene_page.pl?dictybaseid={gene_id}", name="Dicty Base Id")
    schema = [schema_link, "Name", "Synonyms", "Gene Products"]
    
    ret = []
    for gene in genes:
        gene = name_matcher.umatch(gene)
        gi = info.info.get(gene, None)
        if gi:
            ret.append([schema_link.format(gene_id=gene, text=gene),
                        gi[0] + " (%s)" % gene if gene != gi[0] else gi[0], # Gene Name
                        ", ".join(gi[1]), # Synonyms
                        gi[2] or "", # Gene Products
                        ])
            
        else:
            ret.append(None)
    
    return schema, ret
    
    
INFO_SOURCES = {"default": [("NCBI Info", ncbi_info)],
                "352472": [("NCBI Info", ncbi_info),
                           ("Dicty Base", dicty_info)
                           ]
                }

class OWGeneInfo(OWWidget):
    settingsList = ["organismIndex", "geneAttr", "useAttr", "autoCommit"]
    contextHandlers = {"":DomainContextHandler("", ["organismIndex",
                                "geneAttr", "useAttr", "useAltSource"])}
    def __init__(self, parent=None, signalManager=None, name="Gene Info"):
        OWWidget.__init__(self, parent, signalManager, name)

        self.inputs = [("Examples", ExampleTable, self.setData)]
        self.outputs = [("Selected Examples", ExampleTable)]

        self.organismIndex = 0
        self.geneAttr = 0
        self.useAttr = False
        self.autoCommit = False
        self.searchString = ""
        self.selectionChangedFlag = False
        self.useAltSource = 0
        self.loadSettings()
        
        self.infoLabel = OWGUI.widgetLabel(OWGUI.widgetBox(self.controlArea,
                                                    "Info", addSpace=True),
                                           "No data on input\n")
        self.organisms = sorted(set([name.split(".")[-2] for name in \
                            orngServerFiles.listfiles("NCBI_geneinfo")] + \
                            obiGene.NCBIGeneInfo.essential_taxids()))
    
        self.organismBox = OWGUI.widgetBox(self.controlArea, "Organism",
                                           addSpace=True)
        self.organismComboBox = OWGUI.comboBox(self.organismBox, self,
                                "organismIndex", "Organism",
                                items=[obiTaxonomy.name(id) for id in self.organisms],
                                callback=self.setItems,
                                debuggingEnabled=0)
        
        # For now only support one alt source, with a checkbox
        # In the future this can be extended to multiple selections
        self.altSourceCheck = OWGUI.checkBox(self.organismBox, self,
                            "useAltSource", "Show information from dictyBase",
                            callback=self.onAltSourceChange,
#                            debuggingEnabled=0,
                            )
        self.altSourceCheck.hide()
        
        box = OWGUI.widgetBox(self.controlArea, "Gene names", addSpace=True)
        self.geneAttrComboBox = OWGUI.comboBox(box, self, "geneAttr",
                                "Gene atttibute", callback=self.setItems)
        
        c = OWGUI.checkBox(box, self, "useAttr", "Use attribute names",
                           callback=self.setItems,
                           disables=[(-1, self.geneAttrComboBox)])
        
        self.geneAttrComboBox.setDisabled(bool(self.useAttr))

        box = OWGUI.widgetBox(self.controlArea, "Commit", addSpace=True)
        b = OWGUI.button(box, self, "Commit", callback=self.commit)
        c = OWGUI.checkBox(box, self, "autoCommit", "Commit on change")
        OWGUI.setStopper(self, b, c, "selectionChangedFlag",
                         callback=self.commit)
        
        ## A label for dictyExpress link
        self.dictyExpressBox = OWGUI.widgetBox(self.controlArea, "Dicty Express")
        self.linkLabel = OWGUI.widgetLabel(self.dictyExpressBox, "")
        self.linkLabel.setOpenExternalLinks(False)
        self.connect(self.linkLabel, SIGNAL("linkActivated(QString)"),
                     self.onDictyExpressLink)
        self.dictyExpressBox.hide()
        
        OWGUI.rubber(self.controlArea)

        OWGUI.lineEdit(self.mainArea, self, "searchString", "Filter",
                       callbackOnType=True, callback=self.searchUpdate)
        
        self.treeWidget = QTreeView(self.mainArea)
        self.treeWidget.setRootIsDecorated(False)
        self.treeWidget.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.treeWidget.setItemDelegate(LinkStyledItemDelegate(self.treeWidget))
        #self.connect(self.treeWidget, SIGNAL("itemSelectionChanged()"), self.commitIf)
        self.treeWidget.viewport().setMouseTracking(True)
        self.treeWidget.setSortingEnabled(True)
        self.mainArea.layout().addWidget(self.treeWidget)
        
        box = OWGUI.widgetBox(self.mainArea, "",
                              orientation="horizontal")
        OWGUI.button(box, self, "Select Filtered",
                     callback=self.selectFiltered)
        OWGUI.button(box, self, "Clear Selection",
                     callback=self.treeWidget.clearSelection)
        
        self.resize(1000, 700)        

        self.geneinfo = []
        self.cells = []
        self.row2geneinfo = {}
        self.data = None
        self.currentLoaded = None, None
        self.selectionUpdateInProgress = False
        
    def setData(self, data=None):
        self.closeContext()
        self.data = data
        if data:
            self.geneAttrComboBox.clear()
            self.attributes = [attr for attr in self.data.domain.variables + \
                               self.data.domain.getmetas().values() \
                               if attr.varType in [orange.VarTypes.String,
                                                   orange.VarTypes.Discrete]]
            self.geneAttrComboBox.addItems([attr.name for attr in self.attributes])
            self.openContext("", data)
            self.geneAttr = min(self.geneAttr, len(self.attributes) - 1)
            
            taxid = data_hints.get_hint(self.data, "taxid", "")
            if taxid in self.organisms:
                self.organismIndex = self.organisms.index(taxid)
                
            self.useAttr = data_hints.get_hint(self.data, "genesinrows",  self.useAttr)
            
            self.setItems()
        else:
            self.clear()

    def infoSource(self):
        """ Return the current selected info source getter function from 
        INFO_SOURCES
        """
        org = self.organisms[min(self.organismIndex, len(self.organisms) - 1)]
        if org not in INFO_SOURCES:
            org = "default"
        sources = INFO_SOURCES[org]
        name, func =  sources[min(self.useAltSource, len(sources) - 1)]
        return name, func
        
    def setItems(self):
        self.warning(0)
        if not self.data:
            return
        if self.useAttr:
            genes = [attr.name for attr in self.data.domain.attributes]
        elif self.attributes:
            attr = self.attributes[self.geneAttr]
            genes = [str(ex[attr]) for ex in self.data if not ex[attr].isSpecial()]
        else:
            genes = []
        if not genes:
            self.warning(0, "Could not extract genes from input dataset.")
        self.warning(1)
        org = self.organisms[min(self.organismIndex, len(self.organisms) - 1)]
        source_name, info_getter = self.infoSource()
        info , currorg = self.currentLoaded
        self.error(0)
        
        self.updateDictyExpressLink(genes, show=org == "352472")
        self.altSourceCheck.setVisible(org == "352472")
        
        # get the info for the genes in a separate thread
        self.progressBarInit()
#        call = self.asyncCall(info_getter, (org, genes),
#                              name="Load NCBI Gene Info",
#                              blocking=False)
#        call.connect(call, SIGNAL("progressChanged(float)"), self.progressBarSet, Qt.QueuedConnection)
#        with orngServerFiles.DownloadProgress.setredirect(call.emitProgressChanged):
#            call.__call__()
#            schema, geneinfo = call.get_result()
#        call.__call__()
#        schema, geneinfo = call.get_result()
        with orngServerFiles.DownloadProgress.setredirect(self.progressBarSet):
            schema, geneinfo = info_getter(org, genes)
        self.progressBarFinished()
        # schema, geneinfo = info_getter(org, genes) 

        self.geneinfo = geneinfo = list(zip(genes, geneinfo))

        self.progressBarInit()
        milestones = set([i for i in range(0, len(geneinfo), max(len(geneinfo)/100, 1))])
        self.cells = cells = []
        self.row2geneinfo = {}
        links = []
        for i, (gene, gi) in enumerate(geneinfo):
            if gi:
                row = []
                for sch, item in zip(schema, gi):
                    if isinstance(item, Link): # TODO: This should be handled by delegates 
                        row.append(item.text)
                        links.append(item.link)
                    else:
                        row.append(item)
                cells.append(row)
                self.row2geneinfo[len(cells) - 1] = i
#                cells.append([gi.gene_id, gi.symbol + " (%s)" % gene if gene != gi.symbol else gi.symbol,
#                            gi.locus_tag or "", gi.chromosome or "", gi.description or "",
#                            ", ".join(gi.synonyms), gi.symbol_from_nomenclature_authority or ""])
#                links.append("http://www.ncbi.nlm.nih.gov/sites/entrez?Db=gene&Cmd=ShowDetailView&TermToSearch=%s" % gi.gene_id)
                

            if i in milestones:
                self.progressBarSet(100.0*i/len(geneinfo))
        model = TreeModel(cells, [str(col) for col in schema], self.treeWidget)
        
        model.setColumnLinks(0, links)
        proxyModel = QSortFilterProxyModel(self)
        proxyModel.setSourceModel(model)
        self.treeWidget.setModel(proxyModel)
        self.connect(self.treeWidget.selectionModel(), SIGNAL("selectionChanged(QItemSelection , QItemSelection )"), self.commitIf)
        for i in range(7):
            self.treeWidget.resizeColumnToContents(i)
            self.treeWidget.setColumnWidth(i, min(self.treeWidget.columnWidth(i), 200))
        self.treeWidget.update()
        self.progressBarFinished()

        self.infoLabel.setText("%i genes\n%i matched NCBI's IDs" % (len(genes), len(cells)))
        self.matchedInfo = len(genes), len(cells)

    def clear(self):
        self.infoLabel.setText("No data on input\n")
        self.treeWidget.setModel(TreeModel([], ["NCBI ID", "Symbol", "Locus Tag",
                                            "Chromosome", "Description", "Synonyms",
                                            "Nomenclature"], self.treeWidget))
        self.geneAttrComboBox.clear()
        self.send("Selected Examples", None)

    def commitIf(self, *args):
        if self.autoCommit and not self.selectionUpdateInProgress:
            self.commit()
        else:
            self.selectionChangedFlag = True

    def commit(self):
        if not self.data:
            return
        model = self.treeWidget.model()
        mapToSource = model.mapToSource
        selectedIds = [self.cells[mapToSource(index).row()][0] for index in self.treeWidget.selectedIndexes()]
        selectedRows = self.treeWidget.selectedIndexes()
        selectedRows = [mapToSource(index).row() for index in selectedRows]
        model = model.sourceModel()
        
        selectedGeneids = [self.row2geneinfo[row] for row in selectedRows]
        selectedIds = [self.geneinfo[i][0] for i in selectedGeneids]
        selectedIds = set(selectedIds)
        gene2row = dict((self.geneinfo[self.row2geneinfo[row]][0], row) \
                        for row in selectedRows)
        
        if self.useAttr:
            def is_selected(attr):
                return attr.name in selectedIds
            attrs = [attr for attr in self.data.domain.attributes if is_selected(attr)]
            domain = orange.Domain(attrs, self.data.domain.classVar)
            domain.addmetas(self.data.domain.getmetas())
            newdata = orange.ExampleTable(domain, self.data)
            self.send("Selected Examples", newdata)
        elif self.attributes:
            attr = self.attributes[self.geneAttr]
            geneinfo = dict(self.geneinfo)
            examples = [ex for ex in self.data if str(ex[attr]) in selectedIds]
            if True:  # Add gene info
                domain = orange.Domain(self.data.domain, self.data.domain.classVar)
                domain.addmetas(self.data.domain.getmetas())
                n_columns = model.columnCount()

                headers = [str(model.headerData(i, Qt.Horizontal, Qt.DisplayRole).toString()) \
                           for i in range(n_columns)]
                new_meta_attrs = [(orange.newmetaid(), orange.StringVariable(name)) \
                                  for name in headers]
                domain.addmetas(dict(new_meta_attrs))
                examples = [orange.Example(domain, ex) for ex in examples]
                for ex in examples:
                    for i, (_, meta) in enumerate(new_meta_attrs):
                        row = gene2row[str(ex[attr])]
                        ex[meta] = str(model.data(model.index(row, i), Qt.DisplayRole).toString())

            if examples:
                newdata = orange.ExampleTable(examples)
            else:
                newdata = None
            self.send("Selected Examples", newdata)
        else:
            self.send("Selected Examples", None)
            
    def rowFiltered(self, row):
        searchStrings = self.searchString.lower().split()
        row = unicode(" ".join(self.cells[row]).lower(), errors="ignore")
        return not all([s in row for s in searchStrings])
    
    def searchUpdate(self):
        if not self.data:
            return
        searchStrings = self.searchString.lower().split()
        index = self.treeWidget.model().sourceModel().index
        mapFromSource = self.treeWidget.model().mapFromSource
        for i, row in enumerate(self.cells):
            row = unicode(" ".join(row).lower(), errors="ignore")
            self.treeWidget.setRowHidden(mapFromSource(index(i, 0)).row(), QModelIndex(), not all([s in row for s in searchStrings]))
        #self.treeWidget.model().setFilterRegExp(QRegExp(self.searchString, Qt.CaseInsensitive, QRegExp.FixedString))
            
    def selectFiltered(self):
        if not self.data:
            return
        itemSelection = QItemSelection()
        
        index = self.treeWidget.model().sourceModel().index
        mapFromSource = self.treeWidget.model().mapFromSource
        for i, row in enumerate(self.cells):
            if not self.rowFiltered(i):
                itemSelection.select(mapFromSource(index(i, 0)), mapFromSource(index(i, 0)))
        self.treeWidget.selectionModel().select(itemSelection, QItemSelectionModel.Select | QItemSelectionModel.Rows)
        
    def sendReport(self):
        from Orange.OrangeWidgets import OWReport
        genes, matched = self.matchedInfo
        info, org = self.currentLoaded
        self.reportRaw("<p>Input: %i genes of which %i (%.1f%%) matched NCBI synonyms<br>Organism: %s<br>Filter: %s</p>" % (genes, matched, 100.0 * matched / genes, obiTaxonomy.name(org), self.searchString))
        self.reportSubsection("Gene list")
        self.reportRaw(reportItemView(self.treeWidget))
        
    def updateDictyExpressLink(self, genes, show=False):
        def fix(ddb):
            if ddb.startswith("DDB"): 
                if not ddb.startswith("DDB_G"):
                    ddb = ddb.replace("DDB", "DDB_G")
                return ddb
            return None 
        if show:
            genes = [fix(gene) for gene in genes if fix(gene)]
            link1 = '<a href="http://dictyexpress.biolab.si/run/index.php?gene=%s">Microarray profile</a>'
            link2 = '<a href="http://dictyexpress.biolab.si/run/index.php?gene=%s&db=rnaseq">RNA-Seq profile</a>'
            self.linkLabel.setText(link1 + "<br/>" + link2)
            
            show = any(genes)
                
        if show:
            self.dictyExpressBox.show()
        else:
            self.dictyExpressBox.hide()

    def onDictyExpressLink(self, link):
        if not self.data:
            return

        selectedIndexes = self.treeWidget.selectedIndexes()
        if not len(selectedIndexes):
            QMessageBox.information(self, "No gene ids selected", "Please select some genes and try again.")
            return
        model = self.treeWidget.model()
        mapToSource = model.mapToSource
        selectedIds = [self.cells[mapToSource(index).row()][0] for index in selectedIndexes]
        selectedRows = self.treeWidget.selectedIndexes()
        selectedRows = [mapToSource(index).row() for index in selectedRows]
        model = model.sourceModel()

        selectedGeneids = [self.row2geneinfo[row] for row in selectedRows]
        selectedIds = [self.geneinfo[i][0] for i in selectedGeneids]
        selectedIds = set(selectedIds)

        def fix(ddb):
            if ddb.startswith("DDB"):
                if not ddb.startswith("DDB_G"):
                    ddb = ddb.replace("DDB", "DDB_G")
                return ddb
            return None

        genes = [fix(gene) for gene in selectedIds if fix(gene)]
        url = str(link) % " ".join(genes)
        QDesktopServices.openUrl(QUrl(url))
            
    def onAltSourceChange(self):
        self.setItems()
        
def reportItemView(view):
    model = view.model()
    return reportItemModel(view, model)
    
def reportItemModel(view, model, index=QModelIndex()):
    if not index.isValid() or model.hasChildren(index):
        columnCount, rowCount = model.columnCount(index), model.rowCount(index)
        if not index.isValid():
            text = '<table>\n<tr>' + ''.join('<th>%s</th>' % model.headerData(i, Qt.Horizontal, Qt.DisplayRole).toString() for i in range(columnCount)) +'</tr>\n'
        else:
#            variant = model.data(index, Qt.DisplayRole)
#            text = '<table' + (' caption="%s"' % variant.toString() if variant.isValid() else '') + '>\n'
            pass
        text += ''.join('<tr>' + ''.join('<td>' + reportItemModel(view, model, model.index(row, column, index)) + '</td>' for column in range(columnCount)) + '</tr>\n' for row in range(rowCount) if not view.isRowHidden(row, index))
        text += '</table>'
        return text
    else:
        variant = model.data(index, Qt.DisplayRole)
        return str(variant.toString()) if variant.isValid() else ""
        
if __name__ == "__main__":
    app = QApplication(sys.argv)
    data = orange.ExampleTable("brown-selected.tab")
    w = OWGeneInfo()
    w.show()
    w.setData(data)
    app.exec_()
    w.saveSettings()
        
        
        
        
