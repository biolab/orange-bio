import gzip
import re
import zipfile
from server_update import *
from server_update.tests.test_miRNA import MiRNATest

from io import BytesIO
from urllib.request import urlopen
from orangecontrib.bio.obimiRNA import toTaxo
import orangecontrib.bio.obiTaxonomy as tax

IDpat = re.compile('ID\s*(\S*)\s*standard;')
ACpat = re.compile('AC\s*(\S*);')
RXpat = re.compile('RX\s*PUBMED;\s(\d*).')
FT1pat = re.compile('FT\s*miRNA\s*(\d{1,}\.\.\d{1,})')
FT2pat = re.compile('FT\s*/accession="(MIMAT[0-9]*)"')
FT3pat = re.compile('FT\s*/product="(\S*)"')
SQpat = re.compile('SQ\s*(.*other;)')
seqpat = re.compile('\s*([a-z\s]*)\s*\d*')


def fastprint(filename, mode, what):
    with open(filename, mode) as file:
        file.write(what)


def format_checker(content):
    if len(re.findall('(ID.*?)ID', content.replace('\n',''))):
        return True
    else:
        print('Uncorrect format of miRBase data-file.')
        return False

    
def get_intoFiles(path, data_webPage):
    
    sections = data_webPage.split('//\n')
    sections.pop()
    
    files = []
    
    for s in sections:
        org = str(re.findall('ID\s*(\S*)\s*standard;', s.splitlines()[0])[0]).split('-')[0]
        fastprint(os.path.join(path, '%s_sections.txt' % org), 'a', s+'//\n')
        
        if not('%s_sections.txt' % org) in files:
            files.append('%s_sections.txt' % org)
            
    content = '\n'.join(list(set(files)))    
    fastprint(os.path.join(path, 'fileList.txt'), 'w', content)
            
    return os.path.join(path, 'fileList.txt')
    

def miRNA_info(path, object, org_name):
    
    address = os.path.join(path, '%s' % object)
    print(address)
    prefix = str(re.findall('(\S*)_sections\.txt', object)[0])

    data_webPage = open(address, 'r+t').read()
    if not data_webPage:
        print('Cannot read %s ' % address)
    else:
        format_checker(data_webPage)
            
        print('I have read: %s' % address)
        sections = data_webPage.split('//\n')
        sections.pop()
        print('Sections found: ', str(len(sections)))
            
        num_s = 0
        # files to write
        fastprint(os.path.join(path, '%s_premiRNA.txt' % prefix), 'w',
                  'preID'+'\t'+'preACC'+'\t'+'preSQ'+'\t'+'matACCs'+'\t'+'pubIDs'+'\t'+'clusters'+'\t'+'web_addr'+'\n')
        fastprint(os.path.join(path, '%s_matmiRNA.txt' % prefix), 'w',
                  'matID'+'\t'+'matACC'+'\t'+'matSQ'+'\t'+'pre_forms'+'\t'+'targets'+'\n')
        dictG = {}
        dictP = {}
            
        for s in sections:
            num_s += 1
            print('section: ', num_s, '/', str(len(sections)))
                            
            pubIDs = []
            matIDs = ''
            matACCs = ''
            preSQ = []
            
            my_ids = []
            my_accs = []
            my_locs = []  # if it's [61..81] you have to take from 60 to 81.
            
            rows = s.splitlines()
                
            for r in rows:
                
                if r[0:2] == 'ID':
                    preID = str(IDpat.findall(r)[0])
                    print(preID)
                        
                elif r[0:2] == 'AC':
                    preACC = str(ACpat.findall(r)[0])
                    web_addr = 'http://www.mirbase.org/cgi-bin/mirna_entry.pl?acc=%s' % preACC
                        
                elif r[0:2] == 'RX' and not(RXpat.findall(r) == []):
                    pubIDs.append(str(RXpat.findall(r)[0]))

                elif r[0:2] == 'FT' and not(FT1pat.findall(r) == []):
                    loc_mat = str(FT1pat.findall(r)[0])
                        
                    if not(loc_mat==[]):
                         my_locs.append(loc_mat)
                
                elif r[0:2] == 'FT' and not(FT2pat.findall(r) == []):
                     mat_acc = str(FT2pat.findall(r)[0])
                        
                     if matACCs == '':
                         matACCs = mat_acc
                     else:
                         matACCs = matACCs + ',' + mat_acc
                            
                     if not(mat_acc == []):
                         my_accs.append(mat_acc)    
                                
                elif r[0:2] == 'FT' and not(FT3pat.findall(r) == []):
                    mat_id = str(FT3pat.findall(r)[0])
                        
                    if matIDs == '':
                         matIDs = mat_id
                    else:
                         matIDs = matIDs + ',' + mat_id     
                        
                    if not(mat_id == []):
                         my_ids.append(mat_id)
                                          
                elif r[0:2] == 'SQ':
                     preSQ_INFO = str(SQpat.findall(r)[0])
                     seq = 'on'
            
                elif r[0:2] == '  ' and seq == 'on':
                     preSQ.append(str(seqpat.findall(r)[0]).replace(' ', ''))
                     
            # cluster search
            clusters = ''
            try:
                mirna_page = urlopen('http://www.mirbase.org/cgi-bin/mirna_entry.pl?acc=%s' % preACC).read().decode()
            except IOError:
                print('miRNA_info Error: Check the address for the miRNA page.')
                pass
            
            clust_check = re.findall('<td class="\S*">(Clustered miRNAs)</td>', mirna_page)
                
            if clust_check != [] and str(clust_check[0]) == 'Clustered miRNAs':    
                 clusters = ','.join(re.findall('<td><a href="/cgi-bin/mirna_entry.pl\?acc=MI\d*">(\S*?)</a></td>',mirna_page))
                      
            if clusters == '':
                clusters = 'None'
            
            # before printing:
            if not pubIDs:
                 pubIDs = 'None'
            else:
                pubIDs = ','.join(pubIDs)
            
            preSQ = ''.join(preSQ)
            
            fastprint(os.path.join(path, '%s_premiRNA.txt' % prefix), 'a',
                      preID+'\t'+preACC+'\t'+preSQ+'\t'+matACCs+'\t'+pubIDs+'\t'+clusters+'\t'+web_addr+'\n')
                
            for tup in zip(my_ids, my_accs, my_locs):
                
                [start, stop] = tup[2].split('..')
                
                if not(tup[0] in dictG):
                    dictG[tup[0]] = []
                
                dictG[tup[0]] = [tup[1], preSQ[int(start)-1:int(stop)]]
                
                if not(tup[0] in dictP):
                    dictP[tup[0]] = []

                dictP[tup[0]].append(preID)
                
        for k, v in dictG.items():
            pre_forms = ','.join(dictP[k]) 
            
            # targets
            targets = 'None'
            if k in TargetScanLib:
                targets = ','.join(TargetScanLib[k])
           
            fastprint(os.path.join(path, '%s_matmiRNA.txt' % prefix), 'a',
                      k+'\t'+v[0]+'\t'+v[1]+'\t'+pre_forms+'\t'+targets+'\n')

        return [os.path.join(path, '%s_matmiRNA.txt' % prefix), os.path.join(path, '%s_premiRNA.txt' % prefix)]


"""
Update files
"""
DOMAIN = 'miRNA'

download_path = sf_local.localpath('downloaded_files')
domain_path = sf_local.localpath(DOMAIN)
create_folder(domain_path)
create_folder(download_path)

org_taxo = [tax.name(id) for id in tax.common_taxids()]

# targets library from TargetScan
tarscan_url = 'http://www.targetscan.org//vert_50//vert_50_data_download/Conserved_Site_Context_Scores.txt.zip'

with urlopen(tarscan_url) as stream:
    content = stream.read()

zip_f = zipfile.ZipFile(BytesIO(content))
arch = zip_f.read(zip_f.namelist()[0]).splitlines()[1:]
arch.pop()  # why this pop?
mirnas = [a.decode().split('\t')[3] for a in arch]
gene_ids = [a.decode().split('\t')[1] for a in arch]
    
TargetScanLib = {}
for m, t in zip(mirnas, gene_ids):
    if not(m in TargetScanLib):
        TargetScanLib[m] = []
    if not(t in TargetScanLib[m]):
        TargetScanLib[m].append(t)


# miRNA library form miRBase
print('\nBuilding miRNA library...')
address = 'ftp://mirbase.org/pub/mirbase/CURRENT/miRNA.dat.gz'

data_webPage = gzip.GzipFile(fileobj=BytesIO(urlopen(address).read())).read()
# decode from bytes to strings
data_webPage = data_webPage.decode()


orgs = [re.findall('ID\s*(\S+?)-\S*\s*standard;',l)[0] for l in data_webPage.splitlines() if l[:2]=='ID']
des = [re.findall('DE\s*(.*)\s\S*.*\sstem[\s|-]loop',l)[0] for l in data_webPage.splitlines() if l[:2]=='DE']

assert len(orgs) == len(des)
orgs_des = dict(zip(orgs, des))

file_org = get_intoFiles(download_path, data_webPage)

miRNA_path = os.path.join(download_path, 'miRNA.txt')
print('miRNA file path: %s' % miRNA_path)
premiRNA_path = os.path.join(download_path, 'premiRNA.txt')
print('pre-miRNA file path: %s' % premiRNA_path)
    
fastprint(miRNA_path, 'w', 'matID'+'\t'+'matACC'+'\t'+'matSQ'+'\t'+'pre_forms'+'\t'+'targets'+'\n')
fastprint(premiRNA_path, 'w', 'preID'+'\t'+'preACC'+'\t'+'preSQ'+'\t'+'matACCs'+'\t'+'pubIDs'+'\t'+'clusters'+'\t'+'web_addr'+'\n')
    
for fx in [l.rstrip() for l in open(file_org).readlines()]:
    if orgs_des[fx.split('_')[0]] in org_taxo:
        end_files = miRNA_info(download_path, fx, orgs_des[fx.split('_')[0]])
            
        for filename in end_files:
            print("Now reading %s..." % filename)
            base_name = os.path.basename(filename)
            org = re.findall('/(\S{3,4})_\S{3}miRNA\.txt', filename)[0]
            type_file = re.findall(org+'_(\S*)miRNA\.txt', filename)[0]
            label = re.findall('/(\S{3,4}_\S{3}miRNA?)\.txt', filename)[0]
    
            org_taxid = str(toTaxo.get(org))
            org = tax.name(str(toTaxo.get(org)))

            if type_file == 'mat':
                TITLE = 'miRNA: {} mature form'.format(org)
                TAGS = ['miRNA'] + tax.shortname(org_taxid)
                fastprint(os.path.join(domain_path, base_name), 'w', open(filename, 'r').read())
                create_info_file(os.path.join(domain_path, base_name), title=TITLE, tags=TAGS)
                print('%s mat created' % org)

                for file_line in open(filename).readlines()[1:]:
                    fastprint(miRNA_path, 'a', file_line)

            elif type_file == 'pre':
                TITLE = 'miRNA: {} pre-form'.format(org)
                TAGS = ['miRNA'] + tax.shortname(org_taxid)
                fastprint(os.path.join(domain_path, base_name), 'w', open(filename, 'r').read())
                create_info_file(os.path.join(domain_path, base_name), title=TITLE, tags=TAGS)

                for file_line in open(filename).readlines()[1:]:
                    fastprint(premiRNA_path, 'a', file_line)
            else:
                print('Check the label.')

fastprint(os.path.join(domain_path, 'miRNA.txt'), 'w', open(miRNA_path, 'r').read())
create_info_file(os.path.join(domain_path, 'miRNA.txt'), title='miRNA: miRNA library', tags=['miRNA'])
print('\nmiRNA.txt uploaded')

fastprint(os.path.join(domain_path, 'premiRNA.txt'), 'w', open(premiRNA_path, 'r').read())
create_info_file(os.path.join(domain_path, 'premiRNA.txt'), title='miRNA: pre-form library', tags=['miRNA'])
print('premiRNA.txt uploaded\n')


helper = SyncHelper(DOMAIN, MiRNATest)
helper.run_tests()
helper.sync_files()

helper.remove_update_folder()
