import random
import obimiRNA

mirnaHSA = obimiRNA.ids('hsa')

for pm in reduce(lambda x,y: x+y, [obimiRNA.get_info(m).pre_forms.split(',') for m in random.sample(mirnaHSA,3)]):                                    
    pre_miRNA = obimiRNA.get_info(pm,type='pre')
    print
    print 'Pre-miRNA name: %s' % pm
    print 'Accession Number: %s' % pre_miRNA.preACC
    print 'Accession Number of mature form(s): %s' % pre_miRNA.matACCs
    print 'PubMed accession number(s): %s' % pre_miRNA.pubIDs
    print 'Pre-miRNAs clustered together with %s: %s' % (pm, pre_miRNA.clusters)
    print 'Link to miRBase: %s' % pre_miRNA.web_addr