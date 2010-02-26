import random
import obiGO
import obimiRNA

annotations = obiGO.Annotations('hsa',obiGO.Ontology())

miRNAs = random.sample(obimiRNA.ids('hsa'),10)

dict_all = obimiRNA.get_GO(miRNAs, annotations, goSwitch=False)
dict_enr = obimiRNA.get_GO(miRNAs, annotations, enrichment=True, goSwitch=False)

dict_tfidf = obimiRNA.filter_GO(dict_all, annotations, reverse=False)

print '#\tmiRNA name\t# All GO terms\t# Enriched GO terms\t# Filtred GO terms\n'
for n,m in enumerate(miRNAs):
    print '%d\t%s\t\t%d\t\t%d\t\t%d' % (n+1,m,len(dict_all[m]),len(dict_enr[m]),len(dict_tfidf[m]))