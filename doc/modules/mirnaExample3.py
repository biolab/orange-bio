import random
import obiGO
import obimiRNA

annotations = obiGO.Annotations('hsa',obiGO.Ontology())
miRNAs = random.sample(obimiRNA.ids('hsa'),10)

print 'miRNA\tNumber of annotations\tGO_IDs\n'
for mi,goList in obimiRNA.get_GO(miRNAs, annotations, goSwitch=False).items():
    if goList:
        print '%s\t%d\t%s' % (mi, len(goList), ','.join(goList[0:4])+'...')