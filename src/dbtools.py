# coding=utf-8
import io
import json
import os

def read_cnr_DBase(csvFilename):
    dic = {}
    print 'Open %s' % csvFilename
    f = open(csvFilename,'r')
    for n,line in enumerate(f):
        line = line[:-1]
        if n == 0:
            colnames = line.split('\t')
        elif n != 0 and line != '' and line.startswith('#') == False:
            row = line.split('\t')
            dic[row[0]]=dict(zip(colnames[1:], row[1:]))
    return dic


def write_csv_DBase(cnrDict, csvFilename):
    colnames = ['key','source','entry_name','alternative_names','nomenclature_trouble','mol_type',
                'dna_accession','blastn_evalue','dna_snp', 'dna_sequence',
                'prot_accession','blastp_evalue','prot_snp', 'prot_sequence',
                'function_grp_names','mechanism_names','operon_grp_name','cluster90_grp_name',
                'comment','taxonomy','reference','curated_by','db_names']

    if os.path.exists(csvFilename):
        rename_file(csvFilename, 'previous_version')
    f = open(csvFilename, 'w')
    f.write('\t'.join(colnames)+'\n')
    keys = cnrDict.keys()
    keys.sort()
    for key in keys:
        line = key
        record = cnrDict[key]
        for colname in colnames[1:]:
            line = line + '\t' + record[colname]
        line = line + '\n'
        try :
            f.write(line)
        except UnicodeEncodeError:
            line = line.encode('ascii','ignore')
            print line
            f.write(line)
    f.close()


def read_json_cnrDBase(jsonDBase_filename):
    jsndbTxt = open(jsonDBase_filename, 'r').readlines()[0]
    jsndbDic = json.loads(jsndbTxt)
    return jsndbDic


def write_json_DBase(jsonDBase, jsonFilename):
    with io.open(jsonFilename, 'w', encoding='utf-8') as f:
        f.write(unicode(json.dumps(jsonDBase, ensure_ascii=False, sort_keys=True)))


def generate_DBkey():
    time.sleep(2)
    y = time.localtime(time.time()).tm_year
    m = time.localtime(time.time()).tm_mon
    j = time.localtime(time.time()).tm_mday
    h = time.localtime(time.time()).tm_hour
    min = time.localtime(time.time()).tm_min
    s = time.localtime(time.time()).tm_sec
    key = '%i%i%i.%i%i%i' % (y,m,j,h,min,s)
    return key


def print_record(dict):
    print 'entry_name: %s' % dict['entry_name'], '\talternative names: %s' % dict['alternative_names'], '\tsource: %s' % dict['source'],'\tcutared by: %s' % dict['curated_by']
    print 'nomenclature_trouble: %s' % dict['nomenclature_trouble']
    print 'comment: %s' % dict['comment']
    print 'dna_accession: %s' % dict['dna_accession'], '\tblastn_evalue: %s' % dict['blastn_evalue'], '\tdna_snp: %s' % dict['dna_snp']
    print 'dna sequence: %s' % dict['dna_sequence']
    print 'protein_accession: %s' % dict['prot_accession'], '\tblastp_evalue: %s' % dict['blastp_evalue'], '\tprot_snp: %s' % dict['prot_snp']
    print 'protein sequence: %s' % dict['prot_sequence']
    print 'function_grp_names: %s' % dict['function_grp_names'], '\tmechanism_names: %s' % dict['mechanism_names']
    print 'cluster90_grp_name: %s' % dict['cluster90_grp_name'], '\toperon_grp_name: %s' % dict['operon_grp_name'], '\ttaxonomy: %s' % dict['taxonomy']
    print 'db_names: %s' % dict['db_names']
    print  ''


def select_record(Dict, selectDict):
    if selectDict.keys() != []:
        selectedDict = {}
        for key in Dict.keys():
            record = Dict[key]
            for item in selectDict:
                comparison = selectDict[item][0]
                target_value = selectDict[item][1]
                value = Dict[item]
                if comparison == '=':
                    if value == target_value:
                        option = True
                    else:
                        option = False
                        break
                elif comparison == '!=':
                    if value != target_value:
                        option = True
                    else:
                        option = False
                        break
                elif comparison == 'in':
                    if target_value in value:
                        option = True
                    else:
                        option = False
                        break
                elif comparison == 'not in':
                    if target_value not in value:
                        option = True
                    else:
                        option = False
                        break
                elif comparison == '>=':
                    if value >= target_value:
                        option = True
                    else:
                        option = False
                        break
                elif comparison == '<=':
                    if value <= target_value:
                        option = True
                    else:
                        option = False
                        break
                elif comparison == '>':
                    if value > target_value:
                        option = True
                    else:
                        option = False
                        break
                elif comparison == '<':
                    if value < target_value:
                        option = True
                    else:
                        option = False
                        break
            if option == 'True':
                selectedDict[key] = record
    else:
        return Dict
    print '\n%i records have been selected\n'
    return selectedDict


def write_DB_faa(cnrDict, filename='cnrDBase.faa', selectDict={}):
    write_cnrDict = select_record(cnrDict, selectDict)
    f = open(filename, 'w')
    keyList = write_cnrDict.keys()
    keyList.sort()
    for key in keyList:
        rec = write_cnrDict[key]
        if rec['prot_snp'] == '':
            search_type = 'blast'
        else:
            search_type = 'snp'
        name = '>' + key + '_' + rec['entry_name'] + ' | search_type: %s' % search_type + " | pass_evalue: " + rec['blastp_evalue'] + ' | snp: ' + rec['prot_snp'] + '\n'
        seq = rec['prot_sequence'] + '\n'
        f.write(name+seq)
    f.close()


def rename_file(filename, tag):
    cmd = 'mv %s %s_%s' % (filename, filename, tag)
    os.system(cmd)


def pairwisealn(sbjt_seq, qury_seq):
    from Bio import pairwise2
    from Bio.SubsMat import MatrixInfo as matlist
    if 'X' in sbjt_seq or '*' in sbjt_seq:
        print '### Probleme in database! Please check:', sbjt_seq
        print exit(1)
    if 'X' in qury_seq:
        print '#X# query: ', qury_seq
        qury_seq = qury_seq.replace('X','')
    if '*' in qury_seq:
        index = qury_seq.index('*')
        print '#*# query:', qury_seq
        qury_seq = qury_seq[:index]
        print '#*# query:', qury_seq

    matrix = matlist.ident
    gap_open = -10
    gap_extend = -1.0
    alns = pairwise2.align.globalds(sbjt_seq, qury_seq, matrix, gap_open, gap_extend)
    #print alns
    aln_sbjt, aln_qury, score, begin, end = alns[0]
    subdicList = []
    snpdicList = []
    subj, subi = '', ''
    difList = [i for i in xrange(len(aln_sbjt)) if aln_sbjt[i] != aln_qury[i]]

    idper = round(100 * ((len(aln_sbjt) - len(difList))/float(len(aln_sbjt))), 2)
    sbjt_cov   = round(100 * (len(aln_sbjt.replace('-',''))/float(len(aln_sbjt))), 2)
    qury_cov   = round(100 * (len(aln_qury.replace('-',''))/float(len(aln_qury))), 2)

    return aln_sbjt, aln_qury, idper, sbjt_cov, qury_cov
