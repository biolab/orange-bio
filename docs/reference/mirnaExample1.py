import random
import obimiRNA

miRNAs = obimiRNA.ids()

print 'miRNA name\tAccession_Number\t\tSequence\t\tPre-forms\n'
for m in random.sample(miRNAs, 10):
    accession = obimiRNA.get_info(m).matACC
    sequence = obimiRNA.get_info(m).matSQ
    preForms = obimiRNA.get_info(m).pre_forms
    print '%s\t%s\t\t%s\t\t%s' % (m, accession, sequence, preForms)