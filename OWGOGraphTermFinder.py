"""
<name>GO Graph</name>
<description>GO Term Finder</description>
<icon>icons/GOGraph.png</icon>
<priority>150</priority>
"""
from __future__ import division                     #replaces the / sign with true division
import orange, math, glob,string
import GOlib ## function needed to handle the GO and annotation
import OWGUI
from qt import *
from qtcanvas import *
from OWWidget import *
from OWGOTermFinder import *
from OWOptions import *
from qttable import *
from qwt import *
import copy, time
from operator import setitem

evidenceColors = [QPen(QColor(0, 200, 0), 2), QPen(QColor(0, 0, 0), 2)]
geneNodeColor = QColor(159,247,228)
nodeGeneLine = QColor(0,0,200)
BodyColor_Default = QColor(255, 225, 10)
BodyCasesColor_Default = QColor(100, 100, 210)       
TabColor_Default = QColor(255, 255, 255)
TabShadow_Default = (Qt.darkGray)
TabEdgeColor_Default = QColor(248, 250, 188)

nodeWidth = 80     #Witdth of node
nodeHeight = 60     #Height of node
nodeGeneWidth = 80     #Witdth of Gene node
nodeGeneHeight = 20     #Height of Gene node
infoTabWidth = 120
infoTabHeight = 120
edgeInfoWidth = 100
stripePos = (20,30) #Position of the legend stripe
stripePosOnCanvas = (20,20) #Position of the legend stripe on canvas
stripeSize = (404,30)
stripeSizeOnCanvas = (204,30)
nudge = 8    #Nudge genes a bit left (or right)from the center of node (where needed)

def dict_povezave(file):                             #funkcija, ki zgradi povezave iz datoteke
    f = open(file,"r")
    povezave = {}                                    #dictionary kamor shranjujem povezave
    buf = f.readline()
    while len(buf)!=0:                               #berem po vrsticah dokler je kaj v datoteki
        loceno = buf.split()                         #split razdeli string v seznam posameznih 'besed' - whitespace doloci sam - tab presledek...
        if povezave.has_key(loceno[2]) == 0:         #Ce kljuc se ne obstaja v slovarju ga dodam in mu priredim vrednost
            povezave[loceno[2]] = [loceno[3]]
        else:
            povezave[loceno[2]] = povezave[loceno[2]]+[loceno[3]] #Ce pa kljuc ze obstaja, mu samo prilepim se eno vrednost
        buf = f.readline()
    f.close()
    return povezave


def hasKey(key,list = []):
    for i in list:
        if i == key:
            return 1
    return 0
    
def mono_ovojnica(node,list=[],dag = {}):         #Funkcija, ki vrne tranzitivno ovojnico dolocenega noda
    if dag.has_key(node):
        for x in dag[node]:
            if hasKey(x[0],list) == 0: #** dodal x[0] - prej x
                list.append(x[0]) #**
            mono_ovojnica(x[0],list,dag) #**
    return list

def dict_ovojnica(povezave = {}):                      #Funkcija, ki vrne slovar vseh tranzitivnih ovojnic
    slovar = {}
    for node in povezave.keys():
        lista = []
        slovar[node] = mono_ovojnica(node,lista,povezave)
    return slovar
    
def dict_node(file):                                    #zgradi slovar opisov posameznih GOtermov
    f = open(file,'r')
    desc_node = {}
    buf = f.readline()
    while len(buf)!=0:
        loceno = buf.split(chr(9))                      #chr(9) = horizontal tab   
        desc_node[loceno[3]] = [loceno[1],loceno[2]]   
        buf = f.readline()
    f.close()
    return desc_node

def dict_GO_gene(file):                             #funkcija, ki zgradi slovar tipa GOterm -> [seznam genov]
    f = open(file,"r")
    GO_gene = {}                                    #dictionary kamor shranjujem povezave
    for i in range(1,8):
        buf = f.readline()
    while len(buf)!=0:                               #berem po vrsticah dokler je kaj v datoteki
        loceno = buf.split(chr(9))                   #split razdeli string v seznam posameznih 'besed' - whitespace doloci sam - tab presledek...
        if GO_gene.has_key(loceno[4]) == 0:          #Ce kljuc se ne obstaja v slovarju ga dodam in mu priredim vrednost
           GO_gene[loceno[4]] = [loceno[9]]
        else:
            GO_gene[loceno[4]] = GO_gene[loceno[4]]+[loceno[9]] #Ce pa kljuc ze obstaja, mu samo prilepim se eno vrednost
        buf = f.readline()
    f.close()
    return GO_gene

def popNext(ready = {}):                                 #Returns the node from priority queue (ofcourse that's the node with longest path)
    if ready <> {}:
      maxLenNode =  ready.keys()                           #maxLength stores the name of node with greatest path length
      maxLenNode = maxLenNode[0]                           #at the beginning the node with the max length is the first one   
      maxLength = ready[maxLenNode]                        #maximum length
      maxLength = maxLength[1]
      for i in ready.keys():
        if ready[i][1] > maxLength:
            maxLength = ready[i][1]
            maxLenNode = i
      return maxLenNode
                 
def rangNodes(start,links = {}):                          #Function finds the apropriate rang for each node
    ready = {start:([start],0)}                           #priority queue
    rang = {}
    while ready <> {}:
       current=popNext(ready)
       if rang.has_key(current):
          if rang[current]<ready[current][1]:
            rang[current] = ready[current][1]
       if rang.has_key(current) == 0:
          rang[current] = ready[current][1]
       if links.has_key(current):
         for i in links[current]:
           ready[i[0]] = (ready[current][0]+[i[0]],ready[current][1]+1)
       del ready[current]
    return rang

def insertDummy(povezave = {}, rang = {}):               #procedure for inserting dummy nodes
    def toLong(povezave = {},rang = {}):                 #returns connections that are spanned over more than one layer
        long = []
        len = 0
        for n in povezave.keys():
          for i in povezave[n]:
            len = (rang[i[0]]-rang[n]) #**
            if len > 1:
               long = long + [(n,i[0],len,i[1])] #**
        return long
     
    def merge(povezave = {}, nr = {}):                    #merge the old LONG edges with new dummy-ones
        for i in nr.keys():
           if povezave.has_key(i) <> 0:
               povezave[i] = povezave[i] + [nr[i]]
           if povezave.has_key(i) == 0:
               povezave[i] = [nr[i]]

    longs = toLong(povezave,rang)                         
    cnt = 0
    
    #print povezave
    for i in longs:
       #print i[0]
       #print 'Dolga: ', i
       #print povezave[i[0]]
       for n in povezave[i[0]]:
          cnt = cnt + 1
          if n[0] == i[1]: #**
             #print 'P:', povezave[i[0]][cnt-1]
             del povezave[i[0]][cnt-1]      #delete old LONG edges
       #print '-----------------'
       cnt=0
    #print povezave  
    #print dummyRel
    cnt = 0                                                #add dummy nodes to rangs
    newRoute = {}
    #print '__________'
    for i in longs:  
        #print i[3]
        newRoute[i[0]] = ('d'+str(cnt),i[3]) #**
        for l in range(rang[i[0]]+1,rang[i[1]]):
          prefix = 'd'+str(cnt)
          rang[prefix] = l
          cnt = cnt + 1
          newRoute[prefix] = ('d'+str(cnt),i[3]) #**
       
        newRoute[prefix] = (i[1],i[3]) #**
        #print newRoute
        #print rang
        merge(povezave,newRoute)
        newRoute = {}
    #print povezave

def barryCenter(povezave = {}, rang = {}):             #Crossing minimization proposed by Sugiyama,Tagawa and Toda
    def edgeReverse(povezave = {}):                    #Create reverse DAG, for faster "father of the node" search
        reverse = {}
        for n in povezave.keys():
            for i in povezave[n]:
                if reverse.has_key(i[0]): #** prej i sedaj i[0]
                    reverse[i[0]] = reverse[i[0]] + [n] #**
                if reverse.has_key(i[0]) == 0: #**
                    reverse[i[0]] = [n] #**
        return reverse     
               
    def rangPack(povezave = {},reverse = {}, rang = {}):    #pack nodes into distinct rangs (type: {rang: {node:(fathers,children,order no.)} }
        pack = {}
        for i in rang.keys():                           #make rang dictionaryes
            pack[rang[i]] = {}
        for i in rang.keys():
           pack[rang[i]][i] = []
        
        #print 'PACK :',pack
        children = [] #** 
        for i in pack.keys():
            for n in pack[i]:
               if reverse.has_key(n):
                 father = reverse[n]
               else: 
                 father = []
               if povezave.has_key(n):
                 for c in povezave[n]: #** dodana for zanka
                    children.append(c) #**staro glej komentar -SPREMEMBA
                 #children = povezave[n]
                 #print 'children: ',children
               else:
                 children = []
               pack[i][n] = [father,children,1]
               children = [] #** treba ga je pocistiti vsak korak
        #print 'PACK :',pack
        return pack
        
    def has_key(key,i=[]):  #pogleda ce imamo kljuc ze v seznamu
        r = 0
        for n in i:
            if n == key:
                r=1
                return r
        return r   
    
    
    reverse = edgeReverse(povezave)       #reverses the DAG for bottom up traversal
    pack = rangPack(povezave,reverse,rang)
    #print pack
    a = []    
    nr = 0
    for i in pack.keys():
        if i <> 0:
            for n in pack[i].keys():
                a.append([len(pack[i][n][0]),n])
            a.sort()
            a.reverse()
            pos = 1
            for n in a:
                n[0] = pos
                pos = pos + 1
            pozicije = [] #zasedene pozicije
            for z in a:
                #print z[1]  #node za katerega racunam barry
                #print pack[i][z[1]][0] #vsi predhodniki noda
                deljivec = len(pack[i][z[1]][0])
                sest = 0
                for m in pack[i][z[1]][0]: #vsak predhodnik posebej
                    sest = sest + pack[i][z[1]][2] #pozicija predhodnika            
                #print "bary:"+ str(sest/deljivec) # to je pa zdej bary center
                if has_key(int(sest/deljivec),pozicije) == 0:
                    pozicije.append(int(sest/deljivec))
                    pack[i][z[1]][2] = int(sest/deljivec)
                    #print pozicije
                else:
                    t = 1 #pomozna spremenljivka ki oscilira
                    stevec = 0
                    odmik = 0
                    original = int(sest/deljivec) #originalna - zeljena pozicija, ki pa jo bomo spreminjali
                    while has_key(original+odmik,pozicije) != 0:
                        stevec = stevec + 1
                        if t == 1:
                            t = 0
                            odmik = stevec * 1
                        if t == 0:
                            t = 1
                            odmik = stevec * (-1)
                        #print odmik
                    pozicije.append(original+odmik)
                    pack[i][z[1]][2] = original+odmik
            #print "Pozicije_"
            #print pozicije
            #print has_key(3,pozicije)
            pozicije = []
            a = []
            #print '__________'
        nr = i #just to know how many levels the graph has
    
    for i in range(nr-1,-1,-1):  # barrycentering from bottom up
        for n in pack[i].keys():
            #print len(pack[i][n][1])
            a.append([len(pack[i][n][1]),n])
        a.sort()
        a.reverse()   
        pos = 1
        for n in a:
            n[0] = pos
            pos = pos + 1
        pozicije = []
        for z in a:
            #print "node:"
            #print z[1]  #node za katerega racunam barry
            #print pack[i][z[1]][0] #vsi predhodniki noda
            deljivec = len(pack[i][z[1]][1])
            sest = 0
            for m in pack[i][z[1]][1]: #vsak predhodnik posebej
                sest = sest + pack[i][z[1]][2] #pozicija predhodnika            
            #print "bary:"+ str(sest/deljivec) # to je pa zdej bary center
            if (sest != 0)&(deljivec != 0):
                if has_key(int(sest/deljivec),pozicije) == 0:
                    pozicije.append(int(sest/deljivec))
                    pack[i][z[1]][2] = int(sest/deljivec)
                    #print pozicije
                else:
                    t = 1 #pomozna spremenljivka ki oscilira
                    stevec = 0
                    odmik = 0
                    original = int(sest/deljivec) #originalna - zeljena pozicija, ki pa jo bomo spreminjali
                    while has_key(original+odmik,pozicije) != 0:
                        stevec = stevec + 1
                        if t == 1:
                            t = 0
                            odmik = stevec * 1
                        if t == 0:
                            t = 1
                            odmik = stevec * (-1)
                        #print odmik
                    pozicije.append(original+odmik)
                    pack[i][z[1]][2] = original+odmik
        #print "Pozicije_"
        #print pozicije
        #print has_key(3,pozicije)
        pozicije = []
        a = []    
    
    #print pack
    return pack
                
def setX(minCross = {}):                   # calculates x positions of nodes in each layer
    def calcBarryF(node,rang,minCross = {}):
        i = rang
        fatherSum = 0                           #the sum of father positions
        fatherDiv = 0                           #number of the fathers
        if minCross.has_key(i-1):
           for n in minCross[i][node][0]:       #all the fathers of the node
              fatherDiv = fatherDiv + 1
              fatherSum = fatherSum + minCross[i-1][n][2]
              #print n
        if fatherDiv == 0:
           fatherDiv = 1
        #print 'father: '+str(fatherSum/fatherDiv)+' FDIV: '+str(fatherDiv)
        fBarry = round(fatherSum/fatherDiv)
        return fBarry
    
    def calcBarryS(node,rang,minCross = {}):
        sonSum = 0                           #the sum of son positions
        sonDiv = 0                           #number of the sons
        i = rang
        if minCross.has_key(i+1):
           for n in minCross[i][node][1]:       #all the fathers of the node
              sonDiv = sonDiv + 1
              sonSum = sonSum + minCross[i+1][n[0]][2]
              #print n
           #print '='
           #print 'son: '+str(fatherSum/fatherDiv)
        #print 'father: '+str(fatherSum/fatherDiv)+' FDIV: '+str(fatherDiv)
        if minCross[i][node][1] == []:
            sBarry = minCross[i][node][2]
        else:
            sBarry = round(sonSum/sonDiv)
        return sBarry   
    
    def priority(minCross = {}):
        prior = {}
        max = 0
        currPrior = 0
        for i in minCross.keys():
            for n in minCross[i]:
                currPrior = len(minCross[i][n][0])+len(minCross[i][n][1]) #priority of the node is determined by the number of it's neigbours
                if max < currPrior:
                    max = currPrior                                       #maximum priority is stored for reassignment of the dummy nodes
                prior[n] = currPrior
            
            #for n in prior.keys():                                        #assign the highest possible priority to dummy nodes
            #   if n[0] == 'd':       
            #      prior[n] = max + 1
            #      max = max + 1
        return prior
        
        
    def nodeOnPos(node,i = []):
        for n in i:
            if n[0] == node:
                return n[1]
        return 0

    def has_key(key,i=[]):  #pogleda ce imamo kljuc ze v seznamu
        r = 0
        #print i
        for n in i:
            #print n
            if n[0] == key:
                r=1
                return r
        return r

    def vriniL(dp,node,i = []):  #vrivanje s prestavljanjem nodov v desno
        prvi = 0
        cnt = 0  
        vrni = []
        for n in i:
            if n[0] == dp:
                break
            cnt = cnt + 1
        free = 0
        for n in range(cnt,prvi,-1):
            if i[n][0]-i[n-1][0] <> 1:
                free = n
                break
        for n in range(free,cnt+1):
            i[n][0] = i[n][0] - 1
        vrni = i + [[dp,node]]
        return vrni

    def vriniD(dp,node,i = []):  #vrivanje s prestavljanjem nodov v levo
        zadnji = len (i)
        cnt = 0  
        vrni = []
        for n in i:
            if n[0] == dp:
                break
            cnt = cnt + 1
        free = 0
        for n in range(cnt,zadnji-2):
            if i[n+1][0]-i[n][0] <> 1:
                free = n
                break
        
        if free == 0:
            free = zadnji-1
        
        for n in range(cnt,free+1):
            i[n][0] = i[n][0] + 1
        
        vrni = i + [[dp,node]]
        return vrni

    prior=priority(minCross)   #Define priority number to each node
    
    order = {}                 #define the order of nodes in layer
    for n in minCross.keys():
        for i in minCross[n].keys():
            if order.has_key(n) == 0:
               order[n] = [[minCross[n][i][2],i]]
            else:
               order[n] = order[n] + [[minCross[n][i][2],i]]
        order[n].sort()
    #print 'Order:'
    #print order
    
    balance = {}                            #balance is built like: {Level:[[priority,node],...]}
    for n in minCross.keys():
        for i in minCross[n].keys():
            if balance.has_key(n):
               balance[n]=balance[n]+[[prior[i],i]]
            else:
               balance[n]=[[prior[i],i]]
               
    for i in balance.keys():   #sort the balance in the way that the nodes with higher prioritys are to the left
        balance[i].sort()
        balance[i].reverse()
    #print 'Balance:'    
    #print balance    
    #print '-------------------------------------'
    #b=1
    #while balance.has_key(b):
    #    for n in balance[b]:
    #        print n
    #    b = b + 1
    for p in range(1,10):
        i=1     #we begin with down-ordering on the second layer (0 is the first!)
        newPos = []
        while balance.has_key(i):
            for n in balance[i]:
                dp = calcBarryF(n[1],i,minCross)    #calculate the desired position
                if has_key(dp,newPos)==0:
                    newPos = newPos + [[dp,n[1]]]
                    newPos.sort()
                    #print newPos
                else:
                    #print minCross[i][n[1]][2]
                    #print minCross[i][nodeOnPos(dp,newPos)][2]
                    #print 'hocem: '+str(dp)
                    if minCross[i][n[1]][2] - minCross[i][nodeOnPos(dp,newPos)][2] > 0:
                        while has_key(dp,newPos):
                            dp = dp + 1
                            if has_key(dp,newPos):
                                if minCross[i][nodeOnPos(dp,newPos)][2] > minCross[i][n[1]][2]:
                                    vriniD(dp,n[1],newPos);
                                    newPos.sort()
                        if has_key(dp,newPos) == 0:
                            newPos = newPos + [[dp,n[1]]]
                            newPos.sort()
                        #print 'dej ga na desno!'
                        #print newPos
                    else:
                        while has_key(dp,newPos):
                            dp = dp - 1
                            if has_key(dp,newPos):
                                if minCross[i][nodeOnPos(dp,newPos)][2] < minCross[i][n[1]][2]:
                                    vriniL(dp,n[1],newPos);
                                    newPos.sort()
                        if has_key(dp,newPos) == 0:
                            newPos = newPos + [[dp,n[1]]]
                            newPos.sort()
                        #print 'dej ga na levo!'
                        #print newPos
            #print newPos
            for t in newPos:
                #print t
                minCross[i][t[1]][2] = t[0]
            newPos = []
            #print '---sep---'
            i = i + 1
        #print i
        maxI = i
         
        i = i - 2
        newPos = []   #Up ordering!
        while balance.has_key(i):
            for n in balance[i]:
                dp = calcBarryS(n[1],i,minCross)    #calculate the desired position
                if has_key(dp,newPos)==0:
                    newPos = newPos + [[dp,n[1]]]
                    newPos.sort()
                    #print newPos
                else:
                    #print minCross[i][n[1]][2]
                    #print minCross[i][nodeOnPos(dp,newPos)][2]
                    #print 'hocem: '+str(dp)
                    if minCross[i][n[1]][2] - minCross[i][nodeOnPos(dp,newPos)][2] > 0:
                        while has_key(dp,newPos):
                            dp = dp + 1
                            if has_key(dp,newPos):
                                if minCross[i][nodeOnPos(dp,newPos)][2] > minCross[i][n[1]][2]:
                                    vriniD(dp,n[1],newPos);
                                    newPos.sort()
                        if has_key(dp,newPos) == 0:
                            newPos = newPos + [[dp,n[1]]]
                            newPos.sort()
                        #print 'dej ga na desno!'
                        #print newPos
                    else:
                        while has_key(dp,newPos):
                            dp = dp - 1
                            if has_key(dp,newPos):
                                if minCross[i][nodeOnPos(dp,newPos)][2] < minCross[i][n[1]][2]:
                                    vriniL(dp,n[1],newPos);
                                    newPos.sort()
                        if has_key(dp,newPos) == 0:
                            newPos = newPos + [[dp,n[1]]]
                            newPos.sort()
                        #print 'dej ga na levo!'
                        #print newPos
            #print newPos
            for t in newPos:
                #print t
                minCross[i][t[1]][2] = t[0]
            newPos = []
            #print '---sep---'
            i = i - 1
    
    i=maxI/2 - 1    #At the end of repositioning we reposition just 3 layers from the bottom - better results
    newPos = []
    while balance.has_key(i):
        for n in balance[i]:
            dp = calcBarryF(n[1],i,minCross)    #calculate the desired position
            if has_key(dp,newPos)==0:
                newPos = newPos + [[dp,n[1]]]
                newPos.sort()
                #print newPos
            else:
                #print minCross[i][n[1]][2]
                #print minCross[i][nodeOnPos(dp,newPos)][2]
                #print 'hocem: '+str(dp)
                if minCross[i][n[1]][2] - minCross[i][nodeOnPos(dp,newPos)][2] > 0:
                    while has_key(dp,newPos):
                        dp = dp + 1
                        if has_key(dp,newPos):
                            if minCross[i][nodeOnPos(dp,newPos)][2] > minCross[i][n[1]][2]:
                                vriniD(dp,n[1],newPos);
                                newPos.sort()
                    if has_key(dp,newPos) == 0:
                        newPos = newPos + [[dp,n[1]]]
                        newPos.sort()
                    #print 'dej ga na desno!'
                    #print newPos
                else:
                    while has_key(dp,newPos):
                        dp = dp - 1
                        if has_key(dp,newPos):
                            if minCross[i][nodeOnPos(dp,newPos)][2] < minCross[i][n[1]][2]:
                                vriniL(dp,n[1],newPos);
                                newPos.sort()
                    if has_key(dp,newPos) == 0:
                        newPos = newPos + [[dp,n[1]]]
                        newPos.sort()
                    #print 'dej ga na levo!'
                    #print newPos
        #print newPos
        for t in newPos:
            #print t
            minCross[i][t[1]][2] = t[0]
        newPos = []
        #print '---sep---'
        i = i + 1
    #print i    
    

def normalize(minCross = {}):   #normalisation of x position of the nodes - makes all x coord's positive, but not too positive (align graph to the left side of canvas)
    min = 1000                  # nodes probably will not have position greater then 1000
    for i in minCross.keys():
        for n in minCross[i].keys():
            if minCross[i][n][2] < min:
                min = minCross[i][n][2]
    for i in minCross.keys():
        for n in minCross[i].keys():
            minCross[i][n][2] = minCross[i][n][2] + min*(-1) 
    
def text_graph(minCross = {}):
    pos = ''
    i=0
    for n in minCross.keys():
       pos = ''
       for l in range(-10,10):
          for i in minCross[n]:
             if minCross[n][i][2] == l:
                 #if i[0] <> 'd':   
                    pos = pos[0:len(pos)-len(str(i))]
                    pos = pos + str(i)
          pos = pos + '     |';
       
       #print pos

#<---------------------------------- UP LOGIC | DOWN GRAPHICS -------------------------------------->

def setTextLine(container=None, canvas=None, text="", width=1,height=1,z=None, color=None,TextOnly=False):
    if TextOnly == False:
        line = QCanvasText(canvas)    
        line.setFont(QFont('Cyrillic',8))
        line.setZ(z)
        if color:
            line.setColor(color)
    test = QCanvasText(canvas)    
    test.setFont(QFont('Cyrillic',8))
    subText = ''   
    whole = text
    r = 0
    breaks = 0
    for i in range(1,len(whole)+1):
        subText =  whole[r:i]
        test.setText(subText)
        if test.boundingRect().right() > width-15:
            if whole[i:] != '':
                whole = whole[0:i]+'-\n'+whole[i:]
                r = i+1
                test.setText(whole[0:len(whole)-2])
                if test.boundingRect().bottom() - test.boundingRect().top() > height - 13:
                    whole = whole[0:r]
                    break
    if TextOnly == False: #return object QCanvasText with the text in it
        line.setText(whole)
        container.append(line)
        return line
    else:
        return whole  #return just formated text
    
class GeneNode(QCanvasRectangle):
    def __init__(self, *args):
        apply(QCanvasRectangle.__init__, (self,) +args)
        self.canvas = args[0]
        self.setZ(6)
        self.setBrush(QBrush(geneNodeColor))
        self.geneLine = QCanvasLine(self.canvas)
        self.geneLine.setPen(QPen(nodeGeneLine, 2))
        self.geneLine.setZ(-2)
        self.geneText = QCanvasText(self.canvas) 
        self.bounding = 0

class Node(QCanvasRectangle):
    def __init__(self, *args):
        apply(QCanvasRectangle.__init__, (self,) +args)
        self.setZ(0)
        self.canvas = args[0]
        self.line = QCanvasLine(self.canvas)
        self.line.setPen(QPen(QColor(0, 0, 0), 1))
        self.line.setZ(5)
        self.text = []
        self.textAnnotation = []
        self.dropplet = QCanvasEllipse(self.canvas)
        self.dropplet.setBrush(QBrush(Qt.gray))
        self.dropplet.setZ(-1) 
        self.selection = []
        for i in range(0,8):
            self.selection.append(QCanvasLine(self.canvas))
        for i in self.selection:
            i.setPen(QPen(QColor(0, 0, 200), 4))
        self.selection[0].setPoints(2,2,2,10)           
        self.selection[1].setPoints(0,0,12,0)
        self.selection[2].setPoints(12,0,0,0)
        self.selection[3].setPoints(10,0,10,10)
        self.selection[4].setPoints(10,12,10,0)
        self.selection[5].setPoints(10,10,0,10)
        self.selection[6].setPoints(1,12,1,1)
        self.selection[7].setPoints(2,10,10,10)
        self.selection[0].setZ(5)           
        self.selection[1].setZ(5)
        self.selection[2].setZ(5)
        self.selection[3].setZ(5)
        self.selection[4].setZ(5)
        self.selection[5].setZ(5)
        self.selection[6].setZ(5)
        self.selection[7].setZ(5)
        self.selected = False       #rectangle around node ?
        self.visible = True         #Node visible ?
        self.droppletClicked = False       #children nodes showed or hided
        self.genes = ''
    
    def setNodeText(self,text,val,gene):
        setTextLine(self.text,self.canvas,text,nodeWidth,nodeHeight,10)
        setTextLine(self.textAnnotation,self.canvas,'p='+val[0:5],nodeWidth,nodeHeight,10)
        if gene: #If there is any direct genes I create direct gene tab to display it/them
            self.genes = GeneNode(self.canvas)
            test = QCanvasText(self.canvas)    #text for boundry testing purposses
            test.setFont(QFont('Cyrillic',8))
            wholeLine = ''
            oldLine = ''
            actualLines = 0
            for geneItem in gene:
                oldLine = wholeLine
                wholeLine = wholeLine + geneItem
                test.setText(wholeLine)
                wholeLine = wholeLine + '\n'
                actualLines = actualLines + 1
                if (test.boundingRect().bottom()-test.boundingRect().top()) > nodeGeneHeight:
                    wholeLine = oldLine + '...('+str(len(gene)-(actualLines-1))+' more) ' + '\n'
                    break
            test.setText(wholeLine[0:len(wholeLine)-1]) #last character is the newline character '\n'
            self.genes.bounding = (test.boundingRect().bottom()-test.boundingRect().top())
            self.genes.geneText.setText(wholeLine[0:len(wholeLine)-2])    
            self.genes.geneText.setZ(7)

    def showNode(self):
        self.show()
        self.text[0].show() #name
        self.textAnnotation[0].show()
        self.line.show()
        if self.childrenL:
            self.dropplet.show()
          
    def resizeNode(self,factor):
        factor = factor * 0.1
        self.setSize(nodeWidth*factor,nodeHeight*factor)    
        self.move(self.posx*factor,self.posy*factor)
        if self.childrenL:
            if self.genes:
                self.genes.setSize(nodeGeneWidth*factor,self.genes.bounding*factor)
                self.genes.geneLine.setPoints(0*factor,(nodeHeight)*factor,0*factor-nudge*factor,(nodeHeight+21)*factor)
                self.genes.geneText.setFont(QFont('Cyrillic',8*factor))
                self.genes.move(self.x()+0*factor-nudge*factor,self.y()+(nodeHeight+21)*factor)
                self.genes.geneText.move(2+self.x()-nudge*factor,self.y()+(nodeHeight+21)*factor)                
                self.genes.geneLine.move(self.x()+(nodeWidth/2)*factor,self.y()+0*factor)                
        else:
            if self.genes: 
                self.genes.setSize(nodeGeneWidth*factor,self.genes.bounding*factor)
                self.genes.geneLine.setPoints(0*factor,(nodeHeight)*factor,0*factor,(nodeHeight+21)*factor)
                self.genes.geneText.setFont(QFont('Cyrillic',8*factor))
                self.genes.move(self.x()+0*factor,self.y()+(nodeHeight+21)*factor)
                self.genes.geneText.move(2+self.x()+0*factor,self.y()+(nodeHeight+21)*factor)
                self.genes.geneLine.move(self.x()+(nodeWidth/2)*factor,self.y()+0*factor)
        
        self.line.setPoints(0*factor,(nodeHeight-15)*factor,nodeWidth*factor,(nodeHeight-15)*factor)
        self.line.move(self.x()+0*factor,self.y()+0*factor) #0,0 are the relative coordinates of the line in the box
        self.text[0].setFont(QFont('Cyrillic',8*factor))
        self.text[0].move(self.x()+2*factor,self.y()+1*factor) #2,1 are the relative coordinates of text in the box
        self.textAnnotation[0].setFont(QFont('Cyrillic',8*factor))
        self.textAnnotation[0].move(self.x()+2,self.y()+(nodeHeight-15)*factor) #2,1 are the relative coordinates of text in the box
        self.dropplet.setSize(15*factor,15*factor)
        self.dropplet.move(self.x()+(nodeWidth/2)*factor,self.y()+nodeHeight*factor)  #30,40 are the relative coordinates of the dropplet
        i = 0
        for p in self.parentObj:
            self.edgesI[i].setPoints(p.x()+((nodeWidth/2)*factor),p.y()+(nodeHeight*factor),self.x()+((nodeWidth/2)*factor),self.y())
            i = i+1
        if self.dummyLine != []:
            self.dummyLine[0].setPoints(self.x()+((nodeWidth/2)*factor),self.y(),self.x()+((nodeWidth/2)*factor),self.y()+(nodeHeight*factor))
        self.selection[0].move(self.x()-5,self.y()-3)
        self.selection[1].move(self.x()-5,self.y()-3)
        self.selection[2].move(self.x()+self.width()-7,self.y()-3)
        self.selection[3].move(self.x()+self.width()-7,self.y()-3)
        self.selection[4].move(self.x()+self.width()-7,self.y()+self.height()-7)
        self.selection[5].move(self.x()+self.width()-7,self.y()+self.height()-7)
        self.selection[6].move(self.x()-4,self.y()+self.height()-7)
        self.selection[7].move(self.x()-4,self.y()+self.height()-7)
        
        
    def verticalResize(self,factor):        #Works fine for both vertical and horizontal space corrections
        factor = factor * 0.1
        self.move(self.posx*factor,self.posy*factor)
        if self.childrenL:
            if self.genes:
                self.genes.setSize(nodeGeneWidth*factor,self.genes.bounding*factor)
                self.genes.geneLine.setPoints(0*factor,(nodeHeight)*factor,0*factor-nudge*factor,(nodeHeight+21)*factor)
                self.genes.geneText.setFont(QFont('Cyrillic',8*factor))
                self.genes.move(self.x()+0*factor-nudge*factor,self.y()+(nodeHeight+21)*factor)
                self.genes.geneText.move(2+self.x()-nudge*factor,self.y()+(nodeHeight+21)*factor)                
                self.genes.geneLine.move(self.x()+(nodeWidth/2)*factor,self.y()+0*factor)                
        else:
            if self.genes: 
                self.genes.setSize(nodeGeneWidth*factor,self.genes.bounding*factor)
                self.genes.geneLine.setPoints(0*factor,(nodeHeight)*factor,0*factor,(nodeHeight+21)*factor)
                self.genes.geneText.setFont(QFont('Cyrillic',8*factor))
                self.genes.move(self.x()+0*factor,self.y()+(nodeHeight+21)*factor)
                self.genes.geneText.move(2+self.x()+0*factor,self.y()+(nodeHeight+21)*factor)
                self.genes.geneLine.move(self.x()+(nodeWidth/2)*factor,self.y()+0*factor)
                
        self.setSize(nodeWidth*factor,nodeHeight*factor)    
        self.line.setPoints(0*factor,(nodeHeight-15)*factor,nodeWidth*factor,(nodeHeight-15)*factor)
        self.line.move(self.x()+0*factor,self.y()+0*factor) #0,0 are the relative coordinates of the line in the box
        self.text[0].setFont(QFont('Cyrillic',8*factor))
        self.text[0].move(self.x()+2*factor,self.y()+1*factor) #2,1 are the relative coordinates of text in the box
        self.textAnnotation[0].setFont(QFont('Cyrillic',8*factor))
        self.textAnnotation[0].move(self.x()+2,self.y()+(nodeHeight-15)*factor) #2,1 are the relative coordinates of text in the box
        self.dropplet.setSize(15*factor,15*factor)
        self.dropplet.move(self.x()+(nodeWidth/2)*factor,self.y()+nodeHeight*factor)  #30,40 are the relative coordinates of the dropplet
        i = 0
        for p in self.parentObj:
            self.edgesI[i].setPoints(p.x()+((nodeWidth/2)*factor),p.y()+(nodeHeight*factor),self.x()+((nodeWidth/2)*factor),self.y())
            i = i+1
        if self.dummyLine != []:
            self.dummyLine[0].setPoints(self.x()+((nodeWidth/2)*factor),self.y(),self.x()+((nodeWidth/2)*factor),self.y()+(nodeHeight*factor))
        self.selection[0].move(self.x()-5,self.y()-3)
        self.selection[1].move(self.x()-5,self.y()-3)
        self.selection[2].move(self.x()+self.width()-7,self.y()-3)
        self.selection[3].move(self.x()+self.width()-7,self.y()-3)
        self.selection[4].move(self.x()+self.width()-7,self.y()+self.height()-7)
        self.selection[5].move(self.x()+self.width()-7,self.y()+self.height()-7)
        self.selection[6].move(self.x()-4,self.y()+self.height()-7)
        self.selection[7].move(self.x()-4,self.y()+self.height()-7)          
    


class Edge(QCanvasLine):
    def __init__(self, *args):
        apply(QCanvasLine.__init__, (self,) +args)
        self.setZ(-3)        
        self.father = ()   #the node, that an edge is connected to from bellow (used for lines highlighting)
        self.son = ()   #the node, that an edge is connected to from above (used for lines highlighting)

class InfoTab(QCanvasRectangle):
    def __init__(self, *args):
        apply(QCanvasRectangle.__init__, (self,) +args)
        self.canvas = args[0]
        self.setSize(infoTabWidth,infoTabHeight)
        self.setZ(100)
        self.setBrush(QBrush(TabColor_Default))
        self.sepLine1 = QCanvasLine(self.canvas)
        self.sepLine1.setPoints(0,0,infoTabWidth,0)
        self.sepLine1.setZ(101)
        self.sepLine2 = QCanvasLine(self.canvas)
        self.sepLine2.setPoints(0,0,infoTabWidth,0)
        self.sepLine2.setZ(101)
        self.titleDG = QCanvasText(self.canvas)
        self.titleDG.setFont(QFont('Cyrillic',8,75))
        self.titleDG.setZ(101)
        self.titleDIG = QCanvasText(self.canvas)
        self.titleDIG.setFont(QFont('Cyrillic',8,75))
        self.titleDIG.setZ(101)
        self.dotsD = QCanvasText(self.canvas)
        self.dotsD.setFont(QFont('Cyrillic',8,87))
        self.dotsD.setZ(101)
        self.dotsDI = QCanvasText(self.canvas)
        self.dotsDI.setFont(QFont('Cyrillic',8,87))
        self.dotsDI.setZ(101)
        self.text = []
        self.direct = []
        self.directIndirect = []   
        self.shadow = QCanvasRectangle(self.canvas)
        self.shadow.setSize(infoTabWidth,infoTabHeight)
        self.shadow.setZ(99)
        self.shadow.setBrush(QBrush(TabShadow_Default))
        self.visible = True
    
    def setText(self,text,textD,textDI):
        if self.text == []: #fill it only if it is empty
            setTextLine(self.text,self.canvas,text,infoTabWidth-20,infoTabHeight,110)
            self.text[0].setFont(QFont('Cyrillic',8,75))
        if self.direct == []:
            for TD in textD:
                setTextLine(self.direct,self.canvas,TD,infoTabWidth,infoTabHeight,110)
        if self.directIndirect == []:
            for TDI in textDI:
                setTextLine(self.directIndirect,self.canvas,TDI,infoTabWidth,infoTabHeight,110)
    
    def showTab(self):
        self.show()
        self.sepLine1.show()
        self.text[0].show()
        if self.direct:
            self.titleDG.show()
        if self.directIndirect:
            self.sepLine2.show()
            self.titleDIG.show()
        self.shadow.show()
    
    def hideTab(self):
        self.hide()
        self.sepLine1.hide()
        self.text[0].hide()
        self.titleDG.hide()
        self.sepLine2.hide()
        self.titleDIG.hide()
        self.dotsD.hide()
        self.dotsDI.hide()
        for item in self.direct:
            item.hide()
        for item in self.directIndirect:
            item.hide()
        self.shadow.hide()
    
    def moveTab(self,x,y):
        showedDI = 0    #how many direct&indirect genes is actually showed in tab
        showedD = 0     #how many direct genes is actually showed in tab
        actualNr = 0    #number of all actually displayed lines
        line = 0        #number of all lines (displayd and not displayed)
        moreD = False   #is "...(X more) " line already displayed ?
        moreDI = False  #is "...(X more) " line already displayed ?
        self.move(x,y)
        self.text[0].move(x+2,y+(line*13))
        s = ''
        s = self.text[0].text()
        nr = s.contains('\n')
        line = line + nr + 1.3
        actualNr = actualNr + nr + 1.3
        self.sepLine1.move(x,y+(line*13))
        if self.direct:
            self.titleDG.setText('Direct Genes ('+str(len(self.direct))+')')
            self.titleDG.move(x+2,y+(line*13))
            line = line + 1
            actualNr = actualNr + 1
        for item in self.direct:
            if ((line+1)*13) < infoTabHeight - 26:
                item.move(x+2,y+(actualNr*13))
                item.show()
                showedD = showedD + 1
                actualNr = actualNr + 1
            else:
                if moreD == False:
                    self.dotsD.setText('...('+str(len(self.direct)-showedD)+' more)')
                    self.dotsD.move(x+2,y+actualNr*13)
                    self.dotsD.show()
                    actualNr = actualNr + 1
                    moreD = True
            line = line + 1
        
        if self.directIndirect:
            self.sepLine2.move(x,y+(actualNr*13))
            self.titleDIG.setText('Indirect ('+str(len(self.directIndirect))+')')
            self.titleDIG.move(x+2,y+(actualNr*13))        
            line = line + 1
            actualNr = actualNr + 1
        for item in self.directIndirect:
            if ((line+1)*13) < infoTabHeight:
                item.move(x+2,y+(actualNr*13))
                item.show()
                showedDI = showedDI + 1
                actualNr = actualNr + 1
            else:
                if moreDI == False:
                    self.dotsDI.setText('...('+str(len(self.directIndirect)-showedDI)+' more)')
                    self.dotsDI.move(x+2,y+actualNr*13)
                    self.dotsDI.show()                
                    actualNr = actualNr + 1
                    moreDI = True
            line = line + 1
        
        if ((line)*13) > infoTabHeight:
            self.setSize(infoTabWidth,actualNr*13)
            self.shadow.setSize(infoTabWidth,actualNr*13)
        else:
            self.setSize(infoTabWidth,line*13+5)
            self.shadow.setSize(infoTabWidth,line*13+5)

        self.shadow.move(x+3,y+3)
        

class edgeInfoTab(QCanvasRectangle):
    def __init__(self, *args):
        apply(QCanvasRectangle.__init__, (self,) +args)
        self.canvas = args[0]
        #self.setSize(edgeInfoWidth,edgeInfoHeight)
        self.setZ(100)
        self.setBrush(QBrush(TabEdgeColor_Default))
        self.edgeInfoText = QCanvasText(self.canvas)
        self.edgeInfoText.setZ(101)
    def edgeInfoShow(self):
        self.show()
        self.edgeInfoText.show()
    def edgeInfoHide(self):
        self.hide()
        self.edgeInfoText.hide()
    def edgeInfoMove(self,x,y):
        self.move(x,y)
        self.edgeInfoText.move(x+2,y)

class MyCanvasView(QCanvasView):
    def __init__(self, *args):
        apply(QCanvasView.__init__,(self,)+args)
        self.canvas = args[0]
        self.viewport().setMouseTracking(True)
        self.last = []
        self.oldLines = []
        self.oldGenes = []
        self.oldColor = {}
    
    def contentsMouseMoveEvent(self, ev):
        try:
            #Part for mouse-over-droopplet detection
            items = filter(lambda ci: ci.z()==-1, self.canvas.collisions(ev.pos())) #makes dropplets dark gray when the mouse is over them
            if len(items) != 0:
                items[0].setBrush(QBrush(Qt.darkGray))
                self.last = items
                self.canvas.update()
            else:
                for i in self.last:
                    i.setBrush(QBrush(Qt.gray))
                    self.canvas.update()
            
            #Part for mouse-over-edge detection
            for i in self.oldLines:    #delete all old red lines
                    if i not in items:
                        if self.oldColor.has_key(i):
                            i.setPen(self.oldColor[i])#QPen(QColor(0, 0, 0), 2))
                            self.canvas.update()
                            i.setZ(-3)
            
            items = []
            items = filter(lambda ci: ci.z()==-3, self.canvas.collisions(ev.pos())) #makes lines red when the mouse is over them
            oldLines = []
            if len(items) != 0:
                items = [items[0]] #only one edge is colored! #Here is some abrakadabra for interconnection of long edges!
            for n in items:                                 #Hmmm...quite a code...I'am not sure I understand it ;)
                if n.fs[1].name[0] == 'd':                  #...yes I do...I think
                    items.append(n.fs[1].edgesO[0])         #get serious! ...Here I append all the lines, that are downwards from the one I am over with the mouse
            for n in items:                                 #..Here I append all the lines, that are upwards from the one I am over with the mouse
                if n.fs[0].name[0] == 'd':
                    items.append(n.fs[0].edgesI[0])            
            dummies = []                                     #dummy nodes are also seen as lines, aren't they?
            for n in items:                                  #...so let's include them!
                if n.fs[0].name[0] == 'd':
                    dummies.append(n.fs[0].dummyLine[0])
            items = items + dummies        
            self.oldColor = {}
            if len(items) != 0:
                self.oldLines = items
                for o in self.oldLines:
                    self.oldColor[o] = o.pen() 
                for l in items:
                    l.setPen(QPen(QColor(255, 0, 0), 2))
                    l.setZ(-2)
            if len(items) != 0:
                self.emit(PYSIGNAL("sigMouseMoveOverEdge"), (ev.pos,items[0]))
            else:
                self.emit(PYSIGNAL("sigMouseNotOverEdge"), ())
            
            
            #Part for mouse-over-node detection
            items = []
            items = filter(lambda ci: ci.z()==0, self.canvas.collisions(ev.pos()))
            if len(items) != 0:
                self.emit(PYSIGNAL("sigMouseMoveOverNode"), (ev.pos,items[0]))
            else:
                self.emit(PYSIGNAL("sigMouseNotOverNode"), ())
            
            self.canvas.update()
        except IndexError:
            pass
        
    def contentsMousePressEvent(self, ev):
        try:
            items = filter(lambda ci: ci.z()==0 or ci.z==1, self.canvas.collisions(ev.pos()))
            if len(items) != 0:
                self.emit(PYSIGNAL("sigMouseButtonPressed"), (ev.pos,items[0]))
                self.emit(PYSIGNAL("sigMouseSelectPressed"), ())
            items = []
            items = filter(lambda ci: ci.z()==-1, self.canvas.collisions(ev.pos()))
            if len(items) != 0:
                self.emit(PYSIGNAL("sigMouseDroppletPressed"), (ev.pos,items[0].father))

        except IndexError:
            pass



DEBUG = 0

class OWGOGraphTermFinder(OWGOTermFinder):
    settingsList = ["subNodes","CBLegend","CBGenes","FilterNumEnabled",
                    "FilterPvalEnabled","FilterDepthEnabled","RBNodeColor","RBNodeAnnotation"]
    def __init__(self, parent=None, name='OWGraphicGoTermFinder'):
        self.callbackDeposit = [] # deposit for OWGUI callback functions
        OWGOTermFinder.__init__(self, parent, name)
        self.inputs = [("Cluster Examples", ExampleTable, self.clusterDataset, 0), ("Reference Examples", ExampleTable, self.referenceDataset, 0)]
        self.outputs = [("Examples", ExampleTable), ("Classified Examples", ExampleTableWithClass)]
        self.goLV.reparent(None,0,QPoint(0,0),False) #remove the textual viewer of DAG
        self.sigTermsTable.reparent(None,0,QPoint(0,0),False) #remove table of significant GOids
        self.tabs.removePage(self.filterTab)
        self.zoom = 10
        self.sizeV = 10
        self.sizeH = 10
        self.subNodes = False  #CheckBox for selecting subnodes(or not)
        self.CBLegend  = True  #CheckBox for show/hide legend        
        self.CBGenes = False
        self.FilterNumEnabled = True        
        self.FilterNumValue = 1
        self.FilterPvalEnabled = True
        self.FilterPvalue = 1.00       
        self.FilterDepthEnabled = True
        self.FilterDepthValue = 8
        self.RBNodeColor = 0
        self.RBNodeAnnotation = 0
        self.maxReference = 0
        self.maxDirectIndirect = 0
        self.maxDirect = 0
        self.setColorPalette()        
        
        self.loadSettings()
      
        #VISUAL SETTINGS TAB
        Visual = QVGroupBox(self)
        OWGUI.hSlider(Visual, self, 'sizeH', box='Horizontal Space', minValue=1, maxValue=20, step=1, callback=self.toggleHorizontalSpaceSlider, ticks=1)
        OWGUI.hSlider(Visual, self, 'sizeV', box='Vertical Space', minValue=1, maxValue=50, step=50, callback=self.toggleVerticalSpaceSlider, ticks=0)
        OWGUI.hSlider(Visual, self, 'zoom', box='Zoom', minValue=1, maxValue=20, step=1, callback=self.toggleZoomSlider, ticks=1)
        OWGUI.checkBox(Visual,self,"CBLegend","Show/Hide legend", callback=self.showHideLegend)
        OWGUI.checkBox(Visual,self,"CBGenes","Show/Hide direct genes", callback=self.showHideGenes)
        OWGUI.radioButtonsInBox(Visual,self,'RBNodeColor',['P value','Direct map','Direct&Indirect map','Reference number','None'],box = 'Node color',callback=self.changeColoration)
        OWGUI.radioButtonsInBox(Visual,self,'RBNodeAnnotation',['P value','Direct map','Direct&Indirect map','Reference number'],box = 'Node annotation',callback=self.changeNodeAnnotation)
        self.tabs.insertTab(Visual, "Visual sett.")    
    
        # FILTER TAB
        filterTab1 = QVGroupBox(self)
        box = QVButtonGroup("Filter GO Term Nodes", filterTab1)
        
        OWGUI.checkBox(box, self, 'FilterNumEnabled', "Cluster frequency", callback=self.setFilterNumEnabled)
        self.sliderFilterNumValue = OWGUI.qwtHSlider(box, self, 'FilterNumValue', label='#:', labelWidth=33, minValue=1, maxValue=1000, step=1.0, precision=1, ticks=0, maxWidth=80, callback=self.runMyFilters)
        if not self.FilterNumEnabled:
            self.sliderFilterNumValue.box.setDisabled(1)
        #
        OWGUI.checkBox(box, self, 'FilterPvalEnabled', "p. value", callback=self.setFilterPvalEnabled)
        self.sliderFilterPvalue = OWGUI.qwtHSlider(box, self, 'FilterPvalue', label='#:', labelWidth=33, minValue=0.0, maxValue=1.0, step=0.001, precision=3.0, ticks=0, maxWidth=80, callback=self.runMyFilters)
        if not self.FilterPvalEnabled:
            self.sliderFilterPvalue.box.setDisabled(1)
        #
        OWGUI.checkBox(box, self, 'FilterDepthEnabled', "GO depth", callback=self.setFilterDepthEnabled)
        self.sliderFilterDepthValue = OWGUI.qwtHSlider(box, self, 'FilterDepthValue', label='#:', labelWidth=33, minValue=0.0, maxValue=100, step=1.0, precision=1.0, ticks=0, maxWidth=80, callback=self.runMyFilters)
        if not self.FilterDepthEnabled:
            self.sliderFilterDepthValue.box.setDisabled(1)
        self.tabs.insertTab(filterTab1, "Filter")

        #self.layout=QVBoxLayout(self.mainArea) #spodnji widget ze ima layout (GOTermFinder)
        self.splitter.destroy() #destroy the splitter from GOTermFinder
        self.splitter = QSplitter(QSplitter.Vertical, self.mainArea)
        self.canvas = QCanvas(1,1)
        self.layout.add(self.splitter)
        self.view = MyCanvasView(self.canvas,self.splitter)

        self.info = InfoTab(self.canvas)
        self.edgeInfo = edgeInfoTab(self.canvas)
        
        self.canvasLegend = QCanvas(0, 0)
        self.canvasLegend.resize(700,90)
        self.legend = MyCanvasView(self.canvasLegend,self.splitter)        
        
        self.colorItems = [self.createColorStripe(i,stripeSize[0],stripeSize[1]) for i in range(len(self.ColorPalettes))]
        self.colorItemsOnCanvas = [self.createColorStripe(i,stripeSizeOnCanvas[0],stripeSizeOnCanvas[1]) for i in range(len(self.ColorPalettes))]
        self.drawLegend()

        self.nodes = []
        self.connect(self.view,PYSIGNAL("sigMouseButtonPressed"),self.select)
        self.connect(self.view,PYSIGNAL("sigMouseSelectPressed"), self.viewSelectionChanged)
        self.connect(self.view,PYSIGNAL("sigMouseDroppletPressed"),self.showHide)
        self.connect(self.view,PYSIGNAL("sigMouseMoveOverNode"),self.showInfo)
        self.connect(self.view,PYSIGNAL("sigMouseNotOverNode"),self.hideInfo)        
        self.connect(self.view,PYSIGNAL("sigMouseMoveOverEdge"),self.showEdgeInfo)
        self.connect(self.view,PYSIGNAL("sigMouseNotOverEdge"),self.hideEdgeInfo)

        s = sum(self.splitter.sizes())
        self.splitter.setSizes( [s-10, 10] )
        self.resize(1020,700)
        

    def viewSelectionChanged(self):
        geneToGOterm = {}
        allGOterms = []
        for li in self.nodes:
            if li.selected:
                GOID = li.name #self.goLVitem2GOID.get(li, None)
                GOterm, x, G, pval, genesInGOID, genesInGOIDdirect = self.GOtermValues.get(GOID, (GOID+'?', '', '', '', [], []))
                #print GOterm, x, G, pval, genesInGOID, genesInGOIDdirect
                if GOID == 'root': ## put real aspect instead of 'root'
                    GOterm = self.GO.get('aspect', GOID+'?')
                if GOterm not in allGOterms:
                    allGOterms.append( GOterm)

                ## make gene -> GOterm annotations only for some genes; depending on the selection type
                if self.SelectMode == 1: 
                    geneList = genesInGOIDdirect # node specific: use just the direct annotations
                else:
                    geneList = genesInGOID # subgraph: use both directly and indirectly annotated genes

                for gene in geneList:
                    tmpl = geneToGOterm.get(gene, [])
                    if GOterm not in tmpl:
                        tmpl.append(GOterm)
                        geneToGOterm[gene] = tmpl
        self.sendSelectedGenesData(geneToGOterm, allGOterms)        

    def changeNodeAnnotation(self):
        for i in self.nodes:
            x= 0
            value = self.RBNodeAnnotation  #simulating switch statement
            result = {
                0: lambda x: 'p='+str(i.PVal)[0:5],
                1: lambda x: 'Dir.='+str(i.direct),
                2: lambda x: 'Dir&Indir.='+str(i.directIndirect),
                3: lambda x: 'Ref.='+str(i.reference),
                }[value](x)
            i.textAnnotation[0].setText(result)
        self.canvas.update()
        
    def showEdgeInfo(self,arg,item):
        father =  item.fs[0].name
        son = item.fs[1].name
        while father[0] == 'd':
            father = self.connects[father].parentsL[0]
        while son[0] == 'd':
            son = self.connects[son].childrenL[0][0]
        if father == 'root':
            father = self.GO.get('aspect', 'root'+'?')
            father = setTextLine(text = father,width = edgeInfoWidth,height = 200,TextOnly = True)
        else:
            father = setTextLine(text = self.GOtermValues[father][0],width = edgeInfoWidth,height = 200,TextOnly = True)
        son = setTextLine(text = self.GOtermValues[son][0],width = edgeInfoWidth,height = 200,TextOnly = True)
        if item.type == 0:
            self.edgeInfo.edgeInfoText.setText(son+'\nis a\n'+father)
        else:
            self.edgeInfo.edgeInfoText.setText(son+'\nis part of\n'+father)
        self.edgeInfo.setSize(edgeInfoWidth,(self.edgeInfo.edgeInfoText.boundingRect().bottom()-self.edgeInfo.edgeInfoText.boundingRect().top())+4)
        self.edgeInfo.edgeInfoMove(arg().x()+20,arg().y()+20)
        self.edgeInfo.edgeInfoShow()
        self.canvas.update()
    
    def hideEdgeInfo(self):
        self.edgeInfo.edgeInfoHide()
        self.canvas.update()
    
    def createColorStripe(self, palette,dx,dy):
        bmp = chr(252)*dx*2 + reduce(lambda x,y:x+y, [chr(int(i*250/dx)) for i in range(dx)] * (dy-4)) + chr(252)*dx*2 
        image = QImage(bmp, dx, dy, 8, self.ColorPalettes[palette], 256, QImage.LittleEndian)
        pm = QPixmap()
        pm.convertFromImage(image, QPixmap.Color);
        return pm
    
    def setColorPalette(self):
        white = qRgb(255,255,255)
        gray = qRgb(200,200,200)
        self.ColorPalettes = \
          ([qRgb(255-(255.*i), 150, 255.*i) for i in range(250)] + [white]*3 + [qRgb(0., 0., 255.), qRgb(255., 255., 0.), gray],
           [qRgb(0, 255.*i*2/250., 0) for i in range(125, 0, -1)] + [qRgb(255.*i*2/250., 0, 0) for i in range(125)] + [white]*3 + [qRgb(0, 255., 0), qRgb(255., 0, 0), gray],
           [qRgb(255.*i/250., 0, 0) for i in range(250)] + [white]*3 + [qRgb(0., 0, 0), qRgb(255., 0, 0), gray])

    def drawLegend(self): #draws range stripe and legend in canvasLegend
        item = []
        point = []
        point.append(QPoint(1,1))
        item.append(self.colorItems[0])
        self.arrayPM = QCanvasPixmapArray(item,point)
        self.sprite = QCanvasSprite(self.arrayPM,self.canvasLegend)
        self.sprite.show()
        self.sprite.setX(stripePos[0])
        self.sprite.setY(stripePos[1])
        
        self.titleLegend = QCanvasText(self.canvasLegend)
        self.ticks = []
        for x in range(1,6):
            self.ticks.append(QCanvasLine(self.canvasLegend))
        pos = 0
        for tick in self.ticks:
            if (pos==0)|((pos==4)):
                tick.setPoints(stripePos[0]+pos*(stripeSize[0]/4),stripePos[1]+1,stripePos[0]+pos*(stripeSize[0]/4),stripePos[1]+stripeSize[1]-3)        
            else:
                tick.setPoints(stripePos[0]+pos*(stripeSize[0]/4),stripePos[1]+(stripeSize[1]/1.3),stripePos[0]+pos*(stripeSize[0]/4),stripePos[1]+stripeSize[1]-3)        
            tick.setPen(QPen(QColor(0, 0, 40*pos), 2))
            tick.show()
            tick.setZ(100)
            pos = pos + 1
        self.rngText = []
        for x in range(0,5):
            self.rngText.append(QCanvasText(self.canvasLegend))

        self.markPolyArray = QPointArray(3)
        self.markPolyArray.setPoint(0, 0,0)
        self.markPolyArray.setPoint(1, 8,0)
        self.markPolyArray.setPoint(2, 4,8)
        self.activeMark = QCanvasPolygon(self.canvasLegend)
        self.activeMark.setPoints(self.markPolyArray)
        self.activeMark.setBrush(QBrush(Qt.black))
        self.activeMark.setZ(111) 
        self.activeMarkText = QCanvasText(self.canvasLegend)
        
    def createLegendOnCanvas(self): #draws range stripe and legend in main canvas
        self.frameOnCanvas = QCanvasRectangle(self.canvas)
        self.frameOnCanvas.setZ(199)
        self.frameOnCanvas.setBrush(QBrush(TabColor_Default))
        self.frameOnCanvas.move(stripePosOnCanvas[0]-15,stripePosOnCanvas[1]-15)
        self.frameOnCanvas.setSize(stripeSizeOnCanvas[0]+30,stripeSizeOnCanvas[1]+85)
        self.frameOnCanvas.show()
        
        self.frameShadowOnCanvas = QCanvasRectangle(self.canvas)
        self.frameShadowOnCanvas.setZ(198)
        self.frameShadowOnCanvas.setBrush(QBrush(TabShadow_Default))
        self.frameShadowOnCanvas.move(stripePosOnCanvas[0]-10,stripePosOnCanvas[1]-10)
        self.frameShadowOnCanvas.setSize(stripeSizeOnCanvas[0]+30,stripeSizeOnCanvas[1]+85)
        self.frameShadowOnCanvas.show()        
        
        itemOnCanvas = []
        pointOnCanvas = []
        pointOnCanvas.append(QPoint(1,1))
        itemOnCanvas.append(self.colorItemsOnCanvas[0])
        self.arrayPMOnCanvas = QCanvasPixmapArray(itemOnCanvas,pointOnCanvas)
        self.spriteOnCanvas = QCanvasSprite(self.arrayPMOnCanvas,self.canvas)
        self.spriteOnCanvas.show()
        self.spriteOnCanvas.setX(stripePosOnCanvas[0])
        self.spriteOnCanvas.setY(stripePosOnCanvas[1])
        self.spriteOnCanvas.setZ(200)
        
        self.titleLegendOnCanvas = QCanvasText(self.canvas)
        self.titleLegendOnCanvas.setZ(200)
        
        self.ticksOnCanvas = []
        for x in range(1,6):
            self.ticksOnCanvas.append(QCanvasLine(self.canvas))
        for x in self.ticksOnCanvas:
            x.setZ(200)
        
        pos = 0
        for tickOnCanvas in self.ticksOnCanvas:
            if (pos==0)|((pos==4)):
                tickOnCanvas.setPoints(stripePosOnCanvas[0]+pos*(stripeSizeOnCanvas[0]/4),stripePosOnCanvas[1]+1,stripePosOnCanvas[0]+pos*(stripeSizeOnCanvas[0]/4),stripePosOnCanvas[1]+stripeSizeOnCanvas[1]-3)        
            else:
                tickOnCanvas.setPoints(stripePosOnCanvas[0]+pos*(stripeSizeOnCanvas[0]/4),stripePosOnCanvas[1]+(stripeSizeOnCanvas[1]/1.3),stripePosOnCanvas[0]+pos*(stripeSizeOnCanvas[0]/4),stripePosOnCanvas[1]+stripeSizeOnCanvas[1]-3)        
            tickOnCanvas.setPen(QPen(QColor(0, 0, 40*pos), 2))
            tickOnCanvas.show()
            tickOnCanvas.setZ(201)
            pos = pos + 1

        self.rngTextOnCanvas = []
        for x in range(0,5):
            self.rngTextOnCanvas.append(QCanvasText(self.canvas))

        self.edgeLegendTextOnCanvas = QCanvasText(self.canvas)
        self.edgeLegendTextOnCanvas.setText('Edge legend:')
        self.edgeLegendTextOnCanvas.setFont(QFont('Cyrillic',8,75))    
        self.edgeLegendTextOnCanvas.move(stripePosOnCanvas[0],stripePosOnCanvas[1]+stripeSizeOnCanvas[1]+15)
        self.edgeLegendTextOnCanvas.show()
        self.edgeLegendTextOnCanvas.setZ(200)
        
        self.isAOnCanvas = QCanvasText(self.canvas)
        self.isAOnCanvas.setText('Is a:')
        self.isAOnCanvas.setFont(QFont('Cyrillic',8,75))    
        self.isAOnCanvas.move(stripePosOnCanvas[0],stripePosOnCanvas[1]+stripeSizeOnCanvas[1]+30)
        self.isAOnCanvas.show()
        self.isAOnCanvas.setZ(200)
        self.isALineOnCanvas = QCanvasLine(self.canvas)    
        self.isALineOnCanvas.setPen(evidenceColors[0])
        self.isALineOnCanvas.setPoints(stripePosOnCanvas[0]+40,stripePosOnCanvas[1]+stripeSizeOnCanvas[1]+38,stripePosOnCanvas[0]+stripeSizeOnCanvas[0],stripePosOnCanvas[1]+stripeSizeOnCanvas[1]+38)
        self.isALineOnCanvas.show()
        self.isALineOnCanvas.setZ(200)
        
        self.partOfOnCanvas = QCanvasText(self.canvas)
        self.partOfOnCanvas.setText('Part of:')
        self.partOfOnCanvas.setFont(QFont('Cyrillic',8,75))    
        self.partOfOnCanvas.move(stripePosOnCanvas[0],stripePosOnCanvas[1]+stripeSizeOnCanvas[1]+50)
        self.partOfOnCanvas.show()
        self.partOfOnCanvas.setZ(200)
        self.partOfLineOnCanvas = QCanvasLine(self.canvas)    
        self.partOfLineOnCanvas.setPen(evidenceColors[1])
        self.partOfLineOnCanvas.setPoints(stripePosOnCanvas[0]+40,stripePosOnCanvas[1]+stripeSizeOnCanvas[1]+58,stripePosOnCanvas[0]+stripeSizeOnCanvas[0],stripePosOnCanvas[1]+stripeSizeOnCanvas[1]+58)
        self.partOfLineOnCanvas.show()
        self.partOfLineOnCanvas.setZ(200)

    def setMarks(self,end,exp,discrete): 
        f = (end / 4)
        lo =''
        hi =''
        mark = ''
        self.rngText[0].setText("0")
        self.rngTextOnCanvas[0].setText("0")
        for x in range(1,5):
            if discrete:
                hi = str(int(f*x))
            else:
                hi = str(f*x)
            self.rngText[x].setText(hi)
            self.rngTextOnCanvas[x].setText(hi)

        pos = 0
        for x in self.rngText:
            x.setZ(110)
            x.setX(stripePos[0]+pos*(stripeSize[0]/4)-((x.boundingRect().right()-x.boundingRect().left())/2))
            x.setY(stripePos[1]+stripeSize[1])
            x.show()
            pos = pos +1        
        
        pos = 0
        
        for x in self.rngTextOnCanvas:
            x.setZ(201)
            x.setX(stripePosOnCanvas[0]+pos*(stripeSizeOnCanvas[0]/4)-((x.boundingRect().right()-x.boundingRect().left())/2))
            x.setY(stripePosOnCanvas[1]+stripeSizeOnCanvas[1])
            x.show()
            pos = pos +1        
        
        self.titleLegend.setText(exp)
        self.titleLegend.setFont(QFont('Cyrillic',8,75))
        self.titleLegend.move(stripePos[0],0)
        self.titleLegend.show()
        
        self.titleLegendOnCanvas.setText(exp)
        self.titleLegendOnCanvas.setFont(QFont('Cyrillic',8,75))
        self.titleLegendOnCanvas.move(stripePosOnCanvas[0], stripePosOnCanvas[1]-14)
        self.titleLegendOnCanvas.setZ(201)
        self.titleLegendOnCanvas.show()

        self.canvasLegend.update()
        self.canvas.update()
        
    def createNodes(self):
        x=0
        y=0
        k = 0.1
        self.connects = {} #connection between logical and physical objects
        name = ''
        if self.dag['root'] != []:
            for i in self.minCross.keys():
                y = i                       #level 
                for n in self.minCross[i].keys():
                    #print 'node ',n,'pos: ', minCross[i][n][2]    
                    x=self.minCross[i][n][2]
                    name = n
                    tmp = Node(self.canvas)
                    self.connects[name] = tmp #keeping connection between name and the graphical object that ascends from it
                    self.nodes.append(tmp)
                    if (name[0] != 'd'):
                        if name == 'root':
                            #print self.GO.get('aspect', 'root'+'?')
                            self.nodes[len(self.nodes)-1].setNodeText(self.GO.get('aspect', 'root'+'?'),str(self.GOtermValues[name][3]),self.GOtermValues[name][5]) #[5] are the direct genes
                        else:    
                            self.nodes[len(self.nodes)-1].setNodeText(self.GOtermValues[name][0],str(self.GOtermValues[name][3]),self.GOtermValues[name][5]) #[5] are the direct genes
                        self.nodes[len(self.nodes)-1].setBrush(QBrush(QColor(self.GOtermValues[name][3]*255, 150,255-self.GOtermValues[name][3]*255 ))) #SPREMENI-CASE odvisen od izbire barvanja!
                        self.nodes[len(self.nodes)-1].cluster = len(self.GOtermValues[name][4]) #number in reference
                        self.nodes[len(self.nodes)-1].level = i #nodes level
                        self.nodes[len(self.nodes)-1].PVal = self.GOtermValues[name][3] #P value of the node
                        self.nodes[len(self.nodes)-1].directIndirect = len(self.GOtermValues[name][5])+len(self.GOtermValues[name][4]) #number of direct and indirect number
                        self.nodes[len(self.nodes)-1].direct = len(self.GOtermValues[name][5]) #direct genes in node                    
                        self.nodes[len(self.nodes)-1].reference = self.GOtermValues[name][2] #number in reference
                        if self.nodes[len(self.nodes)-1].directIndirect > self.maxDirectIndirect:
                            self.maxDirectIndirect = self.nodes[len(self.nodes)-1].directIndirect
                        if self.nodes[len(self.nodes)-1].direct > self.maxDirect:
                            self.maxDirect = self.nodes[len(self.nodes)-1].direct                    
                        if self.nodes[len(self.nodes)-1].reference > self.maxReference:
                            self.maxReference = self.nodes[len(self.nodes)-1].reference
                    else:
                        self.nodes[len(self.nodes)-1].cluster = 0
                        self.nodes[len(self.nodes)-1].setNodeText('dummy','dummy',[])
                        self.nodes[len(self.nodes)-1].level = 0
                        self.nodes[len(self.nodes)-1].PVal = 0
                        self.nodes[len(self.nodes)-1].directIndirect = 0
                        self.nodes[len(self.nodes)-1].direct = 0
                        self.nodes[len(self.nodes)-1].reference = 0
    
                    self.nodes[len(self.nodes)-1].posx = (x*90)+20
                    self.nodes[len(self.nodes)-1].posy = (y*12*self.sizeV)+20
                    self.nodes[len(self.nodes)-1].tmpx = (x*90)+20   #used for resizing
                    self.nodes[len(self.nodes)-1].tmpy = (y*12*self.sizeV)+20  
                    self.nodes[len(self.nodes)-1].name = name
                    self.nodes[len(self.nodes)-1].parentsL = self.minCross[i][n][0] #logical parents
                    self.nodes[len(self.nodes)-1].childrenL = self.minCross[i][n][1] #logical children
                    self.nodes[len(self.nodes)-1].visible = True
                    self.nodes[len(self.nodes)-1].dropplet.father = tmp #used for determing the association of dropplet with node
                    self.nodes[len(self.nodes)-1].parentObj = []    #list of parent graphical objects (node class instances)
                    self.nodes[len(self.nodes)-1].edgeLine = []     #list of edges stars and ends coordinates to the parent nodes
                    self.nodes[len(self.nodes)-1].edgesI = []       #graphical elements that represent edges(that come to the node)
                    self.nodes[len(self.nodes)-1].edgesO = []       #graphical elements that represent edges(that go out of the node)
                    self.nodes[len(self.nodes)-1].dummyLine = []    #line that replaces dummy node
                    self.nodes[len(self.nodes)-1].viewNr = 0        #number tels us how many times the node is selected (if we select parent nodes which are connected to the same children)
                    if self.maxX < self.nodes[len(self.nodes)-1].posx*(self.zoom*k):    #set the canvas size
                        self.maxX = self.nodes[len(self.nodes)-1].posx*(self.zoom*k)
                    if self.maxY < self.nodes[len(self.nodes)-1].posy*(self.zoom*k):
                        self.maxY = self.nodes[len(self.nodes)-1].posy*(self.zoom*k)
            
                    if self.nodes[len(self.nodes)-1].childrenL == []: #genes
                        #print  self.GOtermValues[name][5]
                        self.nodes[len(self.nodes)-1].geneName = self.GOtermValues[name][5]
                    
            for i in self.connects:
                for p in self.connects[i].parentsL:
                    self.connects[i].parentObj.append(self.connects[p])

    def createEdges(self):
        for i in self.nodes:
            for p in self.connects[i.name].parentsL:
                tmp = Edge(self.canvas)
                tmp.setPoints(self.connects[p].posx+30,self.connects[p].posy+45,i.posx+30,i.posy)
                tmp.show()
                tmp.fs = (self.connects[p],i)   #tuple which tells us which node begins and which ends the edge
                for r in self.connects[p].childrenL:
                    if i.name == r[0]:
                        try:
                            eind = self.GO['relationTypes'].keys().index(r[1])
                        except:
                            eind = 0
                        tmp.type = eind
                        tmp.setPen(evidenceColors[eind])
                i.edgesI.append(tmp)
                if i.name[0] == "d":
                    tmp1 =Edge(self.canvas)
                    tmp1.setPoints(i.x()+30,i.y()+46,i.x()+30,i.y())
                    tmp1.setPen(self.connects[i.name].edgesI[0].pen()) #the dummy edge is of fathers color
                    tmp1.type =self.connects[i.name].edgesI[0].type
                    tmp1.show()
                    tmp1.fs = (self.connects[p],i)
                    i.dummyLine.append(tmp1)
        
        for i in self.nodes:
            for n in i.edgesI:
                n.fs[0].edgesO.append(n)
                
    def showHide(self,arg,node):        #Hide if showed, show if hided (click on dropplet)
        lista = []
        ovojnica = mono_ovojnica(node.name,lista,self.dag)              
        ovojnica.sort()
        for i in ovojnica:                      #include dummy node edges and nodes hided by other dropplets
            for p in self.connects[i].parentsL:
                if (p[0] == 'd')or(self.connects[p].visible == False):
                    ovojnica.append(p)
        if node.droppletClicked:
            node.droppletClicked = False
        else:
            node.droppletClicked = True
        
        if node.droppletClicked:    
            for i in ovojnica:
                self.connects[i].visible = False
                for n in self.connects[i].parentObj:
                    n.droppletClicked = True
                self.connects[i].hide()
                self.connects[i].line.hide()
                self.connects[i].text[0].hide()
                self.connects[i].textAnnotation[0].hide()
                self.connects[i].dropplet.hide()
                if self.CBGenes:
                    if self.connects[i].genes:
                        self.connects[i].genes.hide()
                        self.connects[i].genes.geneLine.hide()
                        self.connects[i].genes.geneText.hide()
                for p in self.connects[i].selection:
                    p.hide()
                for p in self.connects[i].edgesI:
                    p.hide()
                for p in self.connects[i].dummyLine:
                    p.hide()
        else:
            for i in ovojnica:        
                self.connects[i].visible = True
                for n in self.connects[i].parentObj:
                    n.droppletClicked = False
                if self.connects[i].name[0] != "d":
                    self.connects[i].show()
                    self.connects[i].line.show()
                    self.connects[i].text[0].show()
                    self.connects[i].textAnnotation[0].show()
                    if self.connects[i].childrenL:      
                        self.connects[i].dropplet.show()
                if self.CBGenes:
                    if self.connects[i].genes:
                        self.connects[i].genes.show()
                        self.connects[i].genes.geneLine.show()
                        self.connects[i].genes.geneText.show()
                if self.connects[i].selected:
                    for p in self.connects[i].selection:
                        p.show()
                for p in self.connects[i].edgesI:
                    p.show()
                for p in self.connects[i].dummyLine:
                    p.show()
        self.canvas.update()

    def showHideLegend(self): # show or hide the legend
        if self.CBLegend:
            self.legend.show()
        else:
           self.legend.hide()
           
    def select(self,arg,node):
        if self.SelectMode == 1:
            node.selection[0].move(node.x()-5,node.y()-3)
            node.selection[1].move(node.x()-5,node.y()-3)
            node.selection[2].move(node.x()+node.width()-7,node.y()-3)
            node.selection[3].move(node.x()+node.width()-7,node.y()-3)
            node.selection[4].move(node.x()+node.width()-7,node.y()+node.height()-7)
            node.selection[5].move(node.x()+node.width()-7,node.y()+node.height()-7)
            node.selection[6].move(node.x()-4,node.y()+node.height()-7)
            node.selection[7].move(node.x()-4,node.y()+node.height()-7)        
            if node.selected != True:
                for i in range(0,8):
                    node.selection[i].show()
                node.viewNr = node.viewNr + 1
                node.selected = True
            else:
                for i in range(0,8):
                    node.selection[i].hide()    
                node.viewNr = node.viewNr - 1    
                node.selected = False
            self.canvas.update()
        else:
            lista = []
            ovojnica = mono_ovojnica(node.name,lista,self.dag)              
            ovojnica.append(node.name)
            ovojnica.sort()
            subRoot = node.selected
            for i in ovojnica:
                tmp = self.connects[i]
                if tmp.name[0] != "d":     #select the node if it isn't dummy
                    tmp.selection[0].move(tmp.x()-5,tmp.y()-3)
                    tmp.selection[1].move(tmp.x()-5,tmp.y()-3)
                    tmp.selection[2].move(tmp.x()+tmp.width()-7,tmp.y()-3)
                    tmp.selection[3].move(tmp.x()+tmp.width()-7,tmp.y()-3)
                    tmp.selection[4].move(tmp.x()+tmp.width()-7,tmp.y()+tmp.height()-7)
                    tmp.selection[5].move(tmp.x()+tmp.width()-7,tmp.y()+tmp.height()-7)
                    tmp.selection[6].move(tmp.x()-4,tmp.y()+tmp.height()-7)
                    tmp.selection[7].move(tmp.x()-4,tmp.y()+tmp.height()-7)        
                    if subRoot != True:  #All the sub-nodes have to act like the root node of the sub-tree
                        tmp.viewNr = tmp.viewNr + 1
                        if tmp.visible:
                            for r in range(0,8):
                                tmp.selection[r].show()
                        tmp.selected = True
                    else:
                        if node.viewNr < 2:        #If node has viewNr bigger than one that means it is subNode-selected more than 1 time-DON'T BREAK THAT
                            tmp.viewNr = tmp.viewNr - 1
                            if tmp.viewNr == 0:
                                for r in range(0,8):
                                    tmp.selection[r].hide()    
                                self.connects[i].selected = False
                    
                    if tmp.viewNr < 0:
                        tmp.viewNr = 0 
            self.canvas.update()
    
    def showInfo(self,arg,node):
        self.info.visible = True
        indirect = [g for g in self.GOtermValues[node.name][4] if g not in self.GOtermValues[node.name][5]]
        if self.GOtermValues[node.name][0] == '?':
            self.info.setText(self.GO.get('aspect', 'root'+'?'),self.GOtermValues[node.name][5],indirect)    
        else:    
            self.info.setText(self.GOtermValues[node.name][0],self.GOtermValues[node.name][5],indirect)
        self.info.moveTab(arg().x()+20,arg().y()+20)
        self.info.showTab()
        r = node.brush().color().red()
        self.activeMark.move(stripePos[0]+r*(stripeSize[0]/255)-4,stripePos[1])
        self.activeMark.show()
        x=''
        value = self.RBNodeColor  #simulating switch statement
        result = {
            0: lambda x: node.PVal,
            1: lambda x: node.direct,
            2: lambda x: node.directIndirect,
            3: lambda x: node.reference,
            4: lambda x: '0'
            }[value](x)
        self.activeMarkText.setText(str(result)[0:8])
        self.activeMarkText.move((stripePos[0]+r*(stripeSize[0]/255))-(((self.activeMarkText.boundingRect().right()-self.activeMarkText.boundingRect().left())/2)),stripePos[1]-13)
        self.activeMarkText.show()
        self.canvas.update()
        self.canvasLegend.update()
        
    def hideInfo(self):
        self.activeMark.hide()
        self.activeMarkText.hide()
        if self.info.visible:
            self.info.hideTab()
            self.info.visible = False
            self.info.text = []
            self.info.direct = []
            self.info.directIndirect = []
        self.canvas.update()
        self.canvasLegend.update()
    
    def toggleZoomSlider(self):
        End = False
        colItems = filter(lambda ci: ci.z()<197, self.frameShadowOnCanvas.collisions(True)) #move DAG if it collides with legend on Canvas
        while (colItems == [])&(not End):
            for i in self.nodes:
                i.posx = i.posx-(10*(20-self.zoom))
                self.maxX = self.maxX - (10*(20-self.zoom))
                i.resizeNode(self.zoom)   
                if i.posx < 40:
                    End = True
                    for i in self.nodes:
                        i.posx = i.posx+(10*(20-self.zoom))
                        self.maxX = self.maxX + (10*(20-self.zoom))
                        i.resizeNode(self.zoom)
            colItems = filter(lambda ci: ci.z()<197, self.frameShadowOnCanvas.collisions(True))

        colItems = filter(lambda ci: ci.z()<197, self.frameShadowOnCanvas.collisions(True)) #move DAG if it collides with legend on Canvas
        while colItems != []:
            for i in self.nodes:
                i.posx = i.posx+(10*(20-self.zoom))
                self.maxX = self.maxX + (10*(20-self.zoom))
                i.resizeNode(self.zoom)   
            colItems = filter(lambda ci: ci.z()<197, self.frameShadowOnCanvas.collisions(True))
        
        k = 0.1
        self.maxX=0;self.maxY=0;
        for i in self.nodes:
            if self.maxX < i.posx*(self.zoom*k):
                self.maxX = i.posx*(self.zoom*k)
            if self.maxY < i.posy*(self.zoom*k):
                self.maxY = i.posy*(self.zoom*k)
            i.resizeNode(self.zoom)
        self.canvas.resize(self.maxX+250,self.maxY+200)
        self.canvas.update()

    def toggleVerticalSpaceSlider(self):
        nodeGeneHeight = (13*self.sizeV)-110 #dynamic gene-node height adjustment
        #print nodeGeneHeight
        for i in self.nodes:
            if i.name[0] != 'd':
                gene = self.GOtermValues[i.name][5]
                if self.GOtermValues[i.name][5]: #If there is any direct genes I create direct gene tab to display it/them
                    test = QCanvasText(self.canvas)    #text for boundry testing purposses
                    test.setFont(QFont('Cyrillic',8))
                    wholeLine = ''
                    oldLine = ''
                    actualLines = 0
                    for geneItem in gene:
                        oldLine = wholeLine
                        wholeLine = wholeLine + geneItem
                        test.setText(wholeLine)
                        wholeLine = wholeLine + '\n'
                        actualLines = actualLines + 1
                        if (test.boundingRect().bottom()-test.boundingRect().top()) > nodeGeneHeight:
                            wholeLine = oldLine + '...('+str(len(gene)-(actualLines-1))+' more) ' + '\n'
                            break
                    test.setText(wholeLine[0:len(wholeLine)-1]) #last character is the newline character '\n'
                    i.genes.bounding = (test.boundingRect().bottom()-test.boundingRect().top())
                    i.genes.geneText.setText(wholeLine[0:len(wholeLine)-2])    
                    i.genes.geneText.setZ(7)        

        k = 0.1
        for i in self.nodes:
            i.posy = i.tmpy * self.sizeV * k # self.sizeV 
        End = False
        colItems = filter(lambda ci: ci.z()<197, self.frameShadowOnCanvas.collisions(True)) #move DAG if it collides with legend on Canvas
        while (colItems == [])&(not End):
            for i in self.nodes:
                i.posy = i.posy-40
                self.maxY = self.maxY - 40
                i.verticalResize(self.zoom)   
                if i.posy < 40:
                    End = True
                    for i in self.nodes:
                        i.posy = i.posy+40
                        self.maxY = self.maxY + 40
                        i.verticalResize(self.zoom)
            colItems = filter(lambda ci: ci.z()<197, self.frameShadowOnCanvas.collisions(True))

        colItems = filter(lambda ci: ci.z()<197, self.frameShadowOnCanvas.collisions(True)) #move DAG if it collides with legend on Canvas
        while colItems != []:
            for i in self.nodes:
                i.posy = i.posy+40
                self.maxY = self.maxY + 40
                i.verticalResize(self.zoom)   
            colItems = filter(lambda ci: ci.z()<197, self.frameShadowOnCanvas.collisions(True))


        self.maxX=0;self.maxY=0;
        for i in self.nodes:
            if self.maxX < i.posx*(self.zoom*k):
                self.maxX = i.posx*(self.zoom*k)
            if self.maxY < i.posy*(self.zoom*k):
                self.maxY = i.posy*(self.zoom*k)
            i.verticalResize(self.zoom)
        self.canvas.resize(self.maxX+250,self.maxY+200)
        self.canvas.update()

    def toggleHorizontalSpaceSlider(self):
        k = 0.1
        for i in self.nodes:
            i.posx = i.tmpx * self.sizeH * k 
        End = False
        colItems = filter(lambda ci: ci.z()<197, self.frameShadowOnCanvas.collisions(True)) #move DAG if it collides with legend on Canvas
        while (colItems == [])&(not End):
            for i in self.nodes:
                i.posx = i.posx-40
                self.maxX = self.maxX - 40
                i.verticalResize(self.zoom)   
                if i.posx < 40:
                    End = True
                    for i in self.nodes:
                        i.posx = i.posx+40
                        self.maxX = self.maxX + 40
                        i.verticalResize(self.zoom)
            colItems = filter(lambda ci: ci.z()<197, self.frameShadowOnCanvas.collisions(True))

        colItems = filter(lambda ci: ci.z()<197, self.frameShadowOnCanvas.collisions(True)) #move DAG if it collides with legend on Canvas
        while colItems != []:
            for i in self.nodes:
                i.posx = i.posx+40
                self.maxX = self.maxX + 40
                i.verticalResize(self.zoom)   
            colItems = filter(lambda ci: ci.z()<197, self.frameShadowOnCanvas.collisions(True))

        self.maxX=0;self.maxY=0;
        for i in self.nodes:
            if self.maxX < i.posx*(self.zoom*k):
                self.maxX = i.posx*(self.zoom*k)
            if self.maxY < i.posy*(self.zoom*k):
                self.maxY = i.posy*(self.zoom*k)
            i.verticalResize(self.zoom)
        self.canvas.resize(self.maxX+250,self.maxY+200)
        self.canvas.update()
                        
    def setFilterNumEnabled(self):
        self.sliderFilterNumValue.box.setDisabled(not self.FilterNumEnabled)
        if not self.FilterNumEnabled:
            self.oldNumVal = self.FilterNumValue
            self.FilterNumValue = 1
        else:
            self.FilterNumValue = 1#self.oldNumVal
        self.runFilters()

    def setFilterPvalEnabled(self):
        self.sliderFilterPvalue.box.setDisabled(not self.FilterPvalEnabled)
        if not self.FilterPvalEnabled:
            self.oldPVal = self.FilterPvalue
            self.FilterPvalNumValue = 1
        else:
            self.FilterPvalue = 1#self.oldPVal
        self.runFilters()

    def setFilterDepthEnabled(self):
        self.sliderFilterDepthValue.box.setDisabled(not self.FilterDepthEnabled)
        if not self.FilterDepthEnabled:
            self.oldDepth = self.FilterDepthValue
            self.FilterDepthValue = self.maxDepth
        else:
            self.FilterDepthValue = 8#self.oldDepth
        self.runFilters()
            
    def hideNode(self,node):
        allNodes = [node.name]
        for i in allNodes:
            for p in self.connects[i].parentsL:
                if (p[0] == 'd'):
                    allNodes.append(p)
        for i in allNodes:
            self.connects[i].visible = False
            for n in self.connects[i].parentObj:
                n.droppletClicked = True
            self.connects[i].hide()
            self.connects[i].line.hide()
            self.connects[i].text[0].hide()
            self.connects[i].textAnnotation[0].hide()
            self.connects[i].dropplet.hide()
            if self.connects[i].genes:
                self.connects[i].genes.geneText.hide()
                self.connects[i].genes.hide()
                self.connects[i].genes.geneLine.hide()            
            for p in self.connects[i].selection:
                p.hide()
            for p in self.connects[i].edgesI:
                p.hide()
            for p in self.connects[i].dummyLine:
                p.hide()        
    
    def showNode(self,node):
        allNodes = [node.name]
        for i in allNodes:
            for p in self.connects[i].parentsL:
                if (p[0] == 'd'):
                    allNodes.append(p)
        for i in allNodes:        
            self.connects[i].visible = True
            for n in self.connects[i].parentObj:
                n.droppletClicked = False
            if self.connects[i].name[0] != "d":
                self.connects[i].show()
                self.connects[i].line.show()
                self.connects[i].text[0].show()
                self.connects[i].textAnnotation[0].show()
                if self.connects[i].genes:
                    if self.CBGenes:
                        self.connects[i].genes.show()
                        self.connects[i].genes.geneLine.show()
                        self.connects[i].genes.geneText.show()
                if self.connects[i].childrenL:      
                    self.connects[i].dropplet.show()
            if self.connects[i].selected:
                for p in self.connects[i].selection:
                    p.show()
            for p in self.connects[i].edgesI:
                p.show()
            for p in self.connects[i].dummyLine:
                p.show()

    def updateDAG(self, updateonly=0):
        if self.dag['root'] == []:
            self.goLV.clear()
            return
        else:
            for i in self.canvas.allItems():  #clean the canvas
                i.setCanvas(None)
            self.minCross = {}
            self.nodes = []
            self.connects = {}
            self.maxDirectIndirect = 0
            self.maxDirect = 0
            self.maxReference = 0
            self.maxX = 0
            self.maxY = 0
            
            rang = rangNodes('root',self.dag)
            insertDummy(self.dag,rang)
            self.minCross=barryCenter(self.dag,rang)
            setX(self.minCross)
            normalize(self.minCross)
    
            self.maxNumInstances = max( [1] + [x for (GOterm, x, G, pval, genesInGOID, genesInGOIDdirect) in self.GOtermValues.values()]) #calculate the maximum number of instances in node
            self.maxDepth = max([x for x in self.minCross.keys()])
            
            self.sliderFilterDepthValue.setValue(self.maxDepth+1) 
            self.sliderFilterDepthValue.setRange(1, self.maxDepth+1, 1)
            
            self.sliderFilterPvalue.setValue(1)
            self.sliderFilterPvalue.setRange(0.0, 1.0, 0.001)
            
            self.sliderFilterNumValue.setValue(1)
            self.sliderFilterNumValue.setRange(1, self.maxNumInstances, 1)
            
            self.createLegendOnCanvas()
            
            nodeGeneHeight = (13*self.sizeV)-110 #dynamic gene-node height adjustment
            self.createNodes()  #Create graphical objects from the program constructed description
            self.createEdges()
            for i in self.nodes:    #Initial resizing of nodes
                i.resizeNode(self.zoom)
            for i in self.nodes:
                if i.name[0] != 'd':    #ce je node dummy ga ne kazem
                    i.showNode()        
            
            if self.CBLegend:
                self.legend.show()
            else:
               self.legend.hide()
                     
            self.showHideGenes()
            self.changeNodeAnnotation()
            self.changeColoration()
            self.toggleZoomSlider()
            self.toggleVerticalSpaceSlider()
    
    def distinct(self,lst): #return distinct values from list
        dict = {}
        set = setitem
        map(set, (dict,)*len(lst), lst, [])
        return dict.keys()
        
    def runMyFilters(self):
        leaves = [] 
        for i in self.nodes:
            if (not i.childrenL)&(i.PVal < self.FilterPvalue):
                leaves.append(i)
        for i in leaves:
            for p in i.parentsL:
                leaves.append(self.connects[p])
        PNodes = {}
        for i in self.distinct(leaves):
            PNodes[i] = '0'
        
        for i in self.nodes:
            if ((i.level >= self.FilterDepthValue)|(i.cluster < self.FilterNumValue))|(not PNodes.has_key(i)):
                self.hideNode(i)
            else:
                self.showNode(i)
        self.canvas.update()
    
    def showHideGenes(self): #show/hide directly annotated genes
        if self.CBGenes:
            for i in self.nodes:
                if i.genes:
                    if i.visible:
                        i.genes.show()
                        i.genes.geneLine.show()
                        i.genes.geneText.show()
        else:
            for i in self.nodes:
                if i.genes:
                    if i.visible:
                        i.genes.hide()
                        i.genes.geneLine.hide()
                        i.genes.geneText.hide()
        self.canvas.update()
    
    def changeColoration(self):
        for i in self.minCross.keys():
            for n in self.minCross[i].keys():
                if self.RBNodeColor == 0:
                    self.connects[n].setBrush(QBrush(QColor(self.connects[n].PVal*255, 150,255-self.connects[n].PVal*255 )))
                    self.setMarks(1,'p. value',False)
                if self.RBNodeColor == 1:
                    self.connects[n].setBrush(QBrush(QColor((self.connects[n].direct/self.maxDirect)*255, 150,255-(self.connects[n].direct/self.maxDirect)*255 )))
                    self.setMarks(self.maxDirect,'Direct genes',True)
                if self.RBNodeColor == 2:
                    self.connects[n].setBrush(QBrush(QColor((self.connects[n].directIndirect/self.maxDirectIndirect)*255, 150,255-(self.connects[n].directIndirect/self.maxDirectIndirect)*255 )))
                    self.setMarks(self.maxDirectIndirect,'Direct&Indirect genes',True)
                if self.RBNodeColor == 3:
                    self.connects[n].setBrush(QBrush(QColor((self.connects[n].reference/self.maxReference)*255, 150,255-(self.connects[n].reference/self.maxReference)*255 )))
                    self.setMarks(self.maxReference,'Number in reference',True)
                if self.RBNodeColor == 4:
                    self.connects[n].setBrush(QBrush(QColor(254,206,29)))
                    self.setMarks(0,'No color',True)
        self.canvas.update()

if __name__=="__main__":
    import orange
    a = QApplication(sys.argv)
    ow = OWGOGraphTermFinder()
    a.setMainWidget(ow)
    d = orange.ExampleTable('dicty2.tab', dontCheckStored=1)
    ow.clusterDataset(d, 0)
    ow.show()
    a.exec_loop()
    ow.saveSettings()
