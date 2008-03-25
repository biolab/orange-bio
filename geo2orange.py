import string
import re
import sys
import os.path

def usage():
   print "%s geo-file-name orange-file-name" % os.path.basename(sys.argv[0])

if len(sys.argv)<>3:
   usage()
   sys.exit(0)
   
inname = sys.argv[1]
outname = sys.argv[2]

#outname = "out.tab"
#inname = "GDS1962.soft.tab"

attvalue = re.compile("(.*) = (.*)")

def geo_process_file(inname, choice=None):
   print "Processing: %s" % inname
   f = open(inname)
   subsets = {}
   data = {}
   genes = []
   while True:
      line = f.readline()
      if "!dataset_table_begin" in line:
         # which sample to include
         samples = []
         for v in subsets.values():
            samples.extend(v)
         for s in samples:
            data[s] = []

         # classification class(sample)
         classification = {}
         for c, ss in subsets.items():
            for s in ss:
               classification[s] = c

         # header with sample names          
         line = f.readline() 
         line = line[:-1]
         header = string.split(line, "\t")
         accept = dict([(h,i) for i,h in enumerate(header[2:]) if h in samples])

         # read the data
         while True:
            line = f.readline() # read data line
            if "!dataset_table_end" in line:
               f.close()
               return data, genes, classification
            line = line[:-1]
            items = string.split(line, "\t")
            genes.append(items[1]) # 1=gene name, 0=spot id
            items = items[2:] # only the data values
            for k,i in accept.items():
               data[k].append(items[i])

         
      if "!subset_description" in line:
         m = attvalue.search(line)
         sname = m.group(2)
         line = f.readline()
         m = attvalue.search(line)
         if not choice or sname in choice:
            subsets[sname] = string.split(m.group(2), ",")
            
      if "!dataset_table_end" in line:
         break
      
   return subsets

#rawdata, geneids = geo_process_file(inname, ['tumor grade IV', 'tumor grade II', 'tumor grade III'])
#rawdata, geneids, classification = geo_process_file("GDS330.soft.txt")
rawdata, geneids, classification = geo_process_file(inname)


genepos = {}
for i, g in enumerate(geneids):
   genepos[g] = genepos.get(g, []) + [i]

fout = open(outname, "w")
print "Saving to:  %s" % outname
fout.write("id\toutcome\t" + "\t".join(geneids) + "\n")
fout.write("string\td" + "\tc"*len(geneids) + "\n")
fout.write("meta\tclass" + "\t"*len(geneids) + "\n")
for s, items in rawdata.items():
   fout.write("%s\t%s\t" % (s, classification[s]) + "\t".join(items) + "\n")
fout.close()

