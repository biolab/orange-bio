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
    settingsList = []

    def __init__(self, parent=None):
        OWWidget.__init__(self, parent,"Rule Selection")

        self.inputs = [("Rules", ExampleTable, self.ruledataset), ("Expression", ExampleTable, self.expressiondataset)]
        self.outputs = [("Expression", ExampleTable)]
        
        self.expression = None
        self.rules = None
        self.showMetas = 1

        # GUI
        self.layout=QVBoxLayout(self.mainArea)
        self.table=QTable(self.mainArea)
        self.table.setSelectionMode(QTable.NoSelection)
        self.layout.add(self.table)
        self.table.hide()
        self.grid.setColStretch(0,0)

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
        if self.rules==None or self.expression==None:
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
        instances = len(self.rules)
        for i in range(instances):
            self.progressBarSet(i*50/instances)
            for j in range(len(self.rules.domain.attributes)):
                self.table.setText(i, j, str(self.rules[i][j]))
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
        self.connect(self.table,SIGNAL("clicked(int, int, int, const QPoint&)"), self.runRule)
        #self.table.setColumnMovingEnabled(1)
        self.table.show()
        self.layout.activate() # this is needed to scale the widget correctly

    def sort(self, col):
        "sorts the table by column col"
        if col == self.sortby-1:
            self.sortby = - self.sortby
        else:
            self.sortby = col+1
        self.table.sortColumn(col, self.sortby>=0, TRUE)

    def runRule(self, row, col, button, point):
        nt = orange.ExampleTable(self.expression.domain)
        if row >= 0 and row < len(self.rules):
            rule = str(self.rules[row]['geneList'])
            rule = eval(rule)
            for e in self.expression:
                DDB = str(e['DDB'])
                if DDB in rule:
                    nt.append( e)
        self.send("Expression", nt)
		
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
