"""
<name>Rules</name>
<description>Rules selects from input gene expression only those that match selected rule.</description>
<icon>icons/Default.png</icon>
<priority>600</priority>
"""

from qttable import *
from OWWidget import *
import OWGUI

##############################################################################

class OWRules(OWWidget):
    settingsList = ['ShowClosestOther']

    def __init__(self, parent=None, signalManager = None):
        OWWidget.__init__(self, parent, signalManager, "Rule Selection")

        self.inputs = [("Rules", ExampleTable, self.ruledataset), ("Expression", ExampleTable, self.expressiondataset)]
        self.outputs = [("Expression", ExampleTable), ("Genes", ExampleTable), ("Motifs", list)]
        
        self.ShowClosestOther = FALSE
        self.expression = None
        self.rules = None
        self.showMetas = 1

        # GUI
#        self.grid=QGridLayout(self,1,1,5)
#        self.grid.setColStretch(0,0)
        self.layout=QVBoxLayout(self.mainArea)
        self.ar = QVGroupBox(self.mainArea)
        OWGUI.checkBox(self.ar, self, 'ShowClosestOther', 'Show closest other expressions', tooltip='', callback=self.runRule)
        self.table = QTable(self.ar)
        self.table.setSelectionMode(QTable.Multi)

        #self.grid.addWidget(self.ar, 0, 0)
        self.layout.add(self.ar)
        self.grid.setColStretch(0,1)
        self.grid.setRowStretch(0,1)

    def ruledataset(self, data):
        self.rules = data
        if not(self.rules and self.expression):
            self.table.hide()
        else:
            self.progressBarInit()
            self.set_table()
            self.progressBarFinished()

    def expressiondataset(self, data):
        self.expression = data
        if not(self.rules and self.expression):
            self.table.hide()
        else:
            self.progressBarInit()
            self.set_table()
            self.progressBarFinished()

    def set_table(self):
        if self.rules == None or self.expression == None:
            return

        print self.rules.domain, self.rules.domain.getmetas()
        cols = len(self.rules.domain.attributes)
        if self.rules.domain.classVar:
            cols += 1
        if self.showMetas:
            m = self.rules.domain.getmetas() # getmetas returns a dictionary
            ml = [(k, m[k]) for k in m]
            ml.sort(lambda x,y: cmp(y[0], x[0]))
            metas = [x[1] for x in ml]
            cols += len(metas)
        self.table.setNumCols(cols)
        self.table.setNumRows(len(self.rules))

        # set the header (attribute names)
        self.header=self.table.horizontalHeader()
        for i in range(len(self.rules.domain.attributes)):
            self.header.setLabel(i, self.rules.domain.attributes[i].name)
        col = len(self.rules.domain.attributes)
        if self.rules.domain.classVar:
            self.header.setLabel(col, self.rules.domain.classVar.name)
            col += 1
        if self.showMetas:
            for (j,m) in enumerate(metas):
                self.header.setLabel(j+col, m.name)

        # set the contents of the table (values of attributes)
        formats = []
        for j in range(len(self.rules.domain.attributes)):
            if self.rules.domain.attributes[j].varType == orange.VarTypes.Continuous:
                formats.append( "%06.3f")
            else:
                formats.append( "%s")
        instances = len(self.rules)
        for i in range(instances):
            self.progressBarSet(i*50/instances)
            for j in range(len(self.rules.domain.attributes)):
                self.table.setText(i, j, formats[j] % (self.rules[i][j]))
        col = len(self.rules.domain.attributes)
        if self.rules.domain.classVar:
            self.progressBarSet(50+i*20/instances)
            for i in range(instances):
                OWGUI.tableItem(self.table, i, col, self.rules[i].getclass().native(), background=QColor(160,160,160))
            col += 1
        mlen = len(metas)
        for (j,m) in enumerate(metas):
            self.progressBarSet(70+j*30/mlen)
            for i in range(instances):
                #print m.name, m.varType, self.data[i][m].valueType, self.data[i][m].varType
                OWGUI.tableItem(self.table, i, j+col, str(self.rules[i][m]), background=QColor(220,220,220))

        # adjust the width of the table
        for i in range(cols):
            self.table.adjustColumn(i)

        # manage sorting (not correct, does not handle real values)
        self.connect(self.header,SIGNAL("clicked(int)"), self.sort)
        self.sortby = 0
        #self.connect(self.table,SIGNAL("clicked(int, int, int, const QPoint&)"), self.runRule)
        self.connect(self.table,SIGNAL("selectionChanged()"), self.runRule)
        #self.table.setColumnMovingEnabled(1)
        self.table.show()
##        self.layout.activate() # this is needed to scale the widget correctly

    def sort(self, col):
        print "sorts the table by column col", col, self.sortby
        if col == self.sortby-1:
            self.sortby = - self.sortby
        else:
            self.sortby = col+1
        self.table.sortColumn(col, self.sortby>=0, TRUE)

    def runRule(self):
        if self.table.numSelections() <= 0:
            self.send("Expression", None)
            self.send("Genes", None)
            self.send("Motifs", [])
            return None

        geneListColumn = self.table.numCols() - 1
        otherGeneListColumn = geneListColumn - 1
        rowsSelected = sum([self.table.isRowSelected(i) for i in range(self.table.numRows())])
        genesToSend = []
        motifsToSend = []
        if rowsSelected == 1: ## if only one row selected, return "original" domain
            multi = 0
            if self.ShowClosestOther:
                newcls = ['other'] + [str(self.table.text(i, 0)) for i in range(self.table.numRows()) if self.table.isRowSelected(i)]
                newclvar = orange.EnumVariable('rule', values = newcls)
                ndom = orange.Domain(self.expression.domain.attributes + [newclvar])
            else:
                ndom = self.expression.domain
        else: ## otherwise create a new domain, where rule is class
            multi = 1
            newcls = [str(self.table.text(i, 0)) for i in range(self.table.numRows()) if self.table.isRowSelected(i)]
            newclvar = orange.EnumVariable('rule', values = newcls)
            ndom = orange.Domain(self.expression.domain.attributes + [newclvar])

        nt = orange.ExampleTable(ndom)
        for cn in range(self.table.numRows()):
            if self.table.isRowSelected(cn):
                geneList = str(self.table.text(cn, geneListColumn)) ##self.rules[row]['geneList'])
                geneList = eval(geneList)
                if self.ShowClosestOther:
                    otherGeneList = str(self.table.text(cn, otherGeneListColumn))
                    otherGeneList = eval(otherGeneList)
                for e in self.expression:
                    DDB = str(e['DDB'])
                    if DDB in geneList:
                       if DDB not in genesToSend: genesToSend.append( DDB)
                       ne = orange.Example(ndom, e)
                       rule = str(self.table.text(cn, 0))
                       motsInRule = rule.split(' and ')
                       if multi:
                           ne.setclass(str(rule))
                       elif self.ShowClosestOther:
                           ne.setclass(str(rule))

                       for mot in motsInRule:
                           if mot not in motifsToSend: motifsToSend.append( mot)
                       nt.append( ne)
                    if not(multi) and self.ShowClosestOther and DDB in otherGeneList:
                       ne = orange.Example(ndom, e)
                       ne.setclass("other")
                       nt.append( ne)

	## make an ExampleTable containing the gene list
	if len(genesToSend) > 0:
		gvar = orange.EnumVariable("geneList", values = genesToSend)
		gdom = orange.Domain([gvar])
		gtable = orange.ExampleTable(orange.Domain(gdom, 0))
		for g in genesToSend:
			e = orange.Example(gdom)
			e[0] = g
			gtable.append( e)
	else:
		gtable = None

        self.send("Expression", nt)
        self.send("Genes", gtable)
        self.send("Motifs", motifsToSend)
        print genesToSend
        print motifsToSend
##############################################################################
# Test the widget, run from DOS prompt
# > python OWDataTable.py)
# Make sure that a sample data set (adult_sample.tab) is in the directory

if __name__=="__main__":
    a = QApplication(sys.argv)
    ow = OWRules()
    a.setMainWidget(ow)

#    d = orange.ExampleTable('adult_sample')
    d = orange.ExampleTable('wtclassed')
    ow.show()
    ow.dataset(d)
    a.exec_loop()
