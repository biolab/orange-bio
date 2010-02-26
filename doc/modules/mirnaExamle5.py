import random
import obimiRNA

miRNAs = random.sample(obimiRNA.ids('hsa'),10)

mirPath_all= obimiRNA.get_pathways(miRNAs,enrichment=False, pathSwitch=False)
mirPath_enr = obimiRNA.get_pathways(miRNAs,enrichment=True, pathSwitch=False)

print 'miRNA_name\t# of pathways\t# of enriched pathways\n'
for m in miRNAs:
    print '%s\t\t%d\t\t%d' % (m,len(mirPath_all[m]),len(mirPath_enr[m]))