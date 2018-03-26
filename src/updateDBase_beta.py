#!/usr/bin/python
# coding=utf-8
# Write by Richard Bonnet
# Date: 20/12/2015

import os
import argparse
import filepaths
import logging
import glob
from Bio import SeqIO
from Bio.Seq import Seq
import Bio
import re
import json
import time
from dbtools import print_record
from dbtools import write_csv_DBase
from dbtools import write_json_DBase
from dbtools import read_cnr_DBase
from dbtools import read_json_cnrDBase
from dbtools import write_DB_faa
from dbtools import pairwisealn
import pandas as pd

# key = dmy.hmins
# {key:{'entry_name': , 'alternative_names': , 'nomenclature_trouble': , 'comment: ', 'curated_by': ,
# 'source': , 'mol_type':'cds'/'dna',
# 'blastn_evalue': , 'dna_snp': , 'dna_seq': , 'dna_accession': , 'dna_length': ,
# 'blastp_evalue': , 'prot_snp': , 'prot_seq': , 'prot_accession': , 'prot_length': ,
# 'reference': , 'cluster90_grp_name: ', 'operon_grp_name': , 'function_grp_names': , 'mechanism_names': , 'taxonomy': code|species, 'db_names'}


path = filepaths.determine_path()
working_directory = os.getcwd()

def read_card_DBase(jsonCARD_filename):
    jsndbTxt = open(jsonCARD_filename, 'r').readlines()[0]
    jsndbDic = json.loads(jsndbTxt)

    cardDict = {}
    excluded_seqDict = {}

    model_id_list = jsndbDic.keys()
    model_id_list.sort()
    modelList = []
    paramList = []
    ARO_categoryList = []
    n = 0
    for d in model_id_list:
        if '_' not in d and 'model_sequences' in jsndbDic[d].keys():
            #try :
            #print jsndbDic[d]['model_sequences']['sequence']
            for seq_nub in jsndbDic[d]['model_sequences']['sequence'].keys():
                rec_id = str(d) + '_' + str(seq_nub)
                entry_name = jsndbDic[d]['ARO_name']
                source = 'CARD_%s' % rec_id
                print ''
                print rec_id
                print 'entry_name:', entry_name, '\tsource: %s' % source
                dna_sequence = jsndbDic[d]['model_sequences']['sequence'][seq_nub]['dna_sequence']['sequence']
                if dna_sequence[0:3] in ['CAT','CAA','CAC'] and dna_sequence[-3:] in ['TTA','CTA']:
                    dna_sequence = Seq(dna_sequence)
                    dna_sequence = dna_sequence.reverse_complement()
                print dna_sequence
                dna_accession = jsndbDic[d]['model_sequences']['sequence'][seq_nub]['dna_sequence']['accession']
                dna_length = len(dna_sequence)
                try:
                    blastn_evalue = jsndbDic[d]['model_param']['blastn_evalue']['param_value']
                except KeyError:
                    blastn_evalue = '1e-30'
                print 'dna_accession:', dna_accession
                print 'blastn_evalue:', blastn_evalue
                print dna_sequence

                prot_sequence = jsndbDic[d]['model_sequences']['sequence'][seq_nub]['protein_sequence']['sequence']
                if prot_sequence == '':
                    dna_sequence = Seq(dna_sequence)
                    dna_sequence, prot_sequence = translate(rec_id, 'card.json', dna_sequence)
                    prot_accession = ''
                else:
                    prot_accession = jsndbDic[d]['model_sequences']['sequence'][seq_nub]['protein_sequence']['GI']
                prot_length = len(prot_sequence)
                try:
                    blastp_evalue = jsndbDic[d]['model_param']['blastp_evalue']['param_value']
                except KeyError:
                    blastp_evalue = '1e-50'
                try:
                    prot_snp = ''
                    keyList = jsndbDic[d]['model_param']['snp']['param_value'].keys()
                    keyList.sort()
                    for key in keyList:
                        prot_snp = prot_snp + ',' + jsndbDic[d]['model_param']['snp']['param_value'][key]
                    prot_snp = prot_snp[1:]
                except KeyError:
                    prot_snp = ''
                print 'prot_accession', prot_accession
                print 'blastp-evalue: %s' % blastp_evalue, '\tblastn-evalue: %s' % blastn_evalue
                print 'snp: %s' % prot_snp
                print prot_sequence

                comment = jsndbDic[d]['ARO_description']
                mechanism_names = []
                try :
                    for ARO_category_cvterm_id in jsndbDic[d]['ARO_category'].keys():
                        mechanism_names.append(jsndbDic[d]['ARO_category'][ARO_category_cvterm_id]['category_aro_name'])
                except AttributeError:
                    x = ''
                mechanism_names.sort()
                mechanism_names = '|'.join(mechanism_names)
                print 'mechanism_names:', mechanism_names
                print 'comment: %s\n' % comment


                data = {'entry_name':entry_name,'alternative_names':'','nomenclature_trouble':'',
                            'comment':comment,'curated_by':'','source':source,'mol_type':'cds',
                            'blastn_evalue':blastn_evalue, 'dna_snp':'', 'dna_sequence':str(dna_sequence),
                            'dna_accession':dna_accession, 'dna_length':dna_length,
                            'blastp_evalue':blastp_evalue, 'prot_snp':prot_snp, 'prot_sequence':str(prot_sequence),
                            'prot_accession':prot_accession, 'prot_length':prot_length,
                            'reference':'', 'cluster90_grp_name':'', 'operon_grp_name':'',
                            'function_grp_names':'', 'mechanism_names':mechanism_names, 'taxonomy':''}
                if prot_sequence == '':
                        excluded_seqDict[rec_id] = data
                else:
                    cardDict[rec_id] = data
            #except KeyError:
            #    print '\n!!! No sequence in %s record of CARD !!!\n' % d
    return cardDict, excluded_seqDict


def read_resfinder_DBase(resfinderList, resfinderNote):
    excluded_seqDict = {}
    resfinderNoteDict = {}
    f = open(resfinderNote, 'r')
    for line in f:
        if not line.startswith('#'):
            id = line.strip().split(':')[0].encode('ascii','recode')
            comment = line.strip().split(':')[2]
            alternative_names = ''
            if 'Alternate name' in comment:
                alternative_names = comment.replace('Alternate name', '').replace('; ','').replace(', ',',').split(',')
                alternative_names.sort()
                alternative_names = '|'.join(alternative_names)
            if 'Amino acid sequences of ' in comment and 'are identical' in comment:
                alternative_names = comment.replace('Amino acid sequences of ','').replace(' are identical','').replace(', and ',' ').replace(' and ', ' ').split(' ')
                alternative_names.sort()
                alternative_names = '|'.join(alternative_names)
            if 'operon' in comment:
                operon_grp_name = comment
            else:
                operon_grp_name = ''
            func  = line.strip().split(':')[1].replace('(',', ').replace(')','').replace(' resistance','').replace(', and ',', ').replace(' and ', ', ').replace('s, ',', ').split(', ')
            func.sort()
            func = '|'.join(func)
            resfinderNoteDict[id]={'function_grp_names':func,'alternative_names':alternative_names,'operon_grp_name':operon_grp_name}

    resfinderDict = {}
    for filename in resfinderList:
        print '\nFILE: %s' % filename
        for rec in SeqIO.parse(open(filename, 'r'),'fasta'):
            source = 'resfinder_' + rec.id
            r = re.compile("(.+)_([0-9]+)_(.+)")
            match = r.match(rec.id)
            if match:
                entry_name = match.group(1).encode('ascii','recode')
                index = match.group(2)
                dna_accession = match.group(3).encode('ascii','recode')
                print '\n%s ' % rec.id
                print 'entry_name: %s\t index: %s\t dna_accession: %s' % (entry_name, index, dna_accession)
                dna_sequence, prot_sequence = translate(rec.id, filename, rec.seq)
                if entry_name in resfinderNoteDict.keys():
                    alternative_names = resfinderNoteDict[entry_name]['alternative_names']
                    function_grp_names = resfinderNoteDict[entry_name]['function_grp_names']
                    operon_grp_name = resfinderNoteDict[entry_name]['operon_grp_name']
                else:
                    alternative_names = ''
                    function_grp_names = os.path.splitext(os.path.basename(filename))[0]
                    operon_grp_name = ''

                r = re.compile("bla(.+)(-[0-9]+)")
                match = r.match(entry_name)
                if match:
                    entry_name = match.group(1)+match.group(2)

                print 'alternative_names:', alternative_names
                print 'function_grp_names:', function_grp_names
                print 'operon_grp_name:', operon_grp_name
                print 'source: %s' % source
                data = {'entry_name':entry_name,'alternative_names':alternative_names,'nomenclature_trouble':'',
                                     'comment':'','curated_by':'','source':source,'mol_type':'cds',
	                                 'blastn_evalue':'1e-30', 'dna_snp':'', 'dna_sequence':str(dna_sequence),
                                     'dna_accession':dna_accession, 'dna_length':len(dna_sequence),
	                                 'blastp_evalue':'1e-50', 'prot_snp':'', 'prot_sequence':str(prot_sequence),
                                     'prot_accession':'', 'prot_length':len(prot_sequence) ,
                                     'reference':'', 'cluster90_grp_name':'', 'operon_grp_name':operon_grp_name,
                                     'function_grp_names':function_grp_names, 'mechanism_names':'', 'taxonomy':''}
                if prot_sequence != '':
                    print dna_sequence
                    print prot_sequence
                    print ''
                    resfinderDict[rec.id] = data

                else:
                    print '!!! Translation Problem in sequence name %s in file %s; sequence not included !!!' % (rec.id, filename)
                    excluded_seqDict[rec.id] = data
            else:
                print '\n!!! Parsing problem in sequence name %s in file %s !!!\n' % (rec.id, filename)
                exit()
    return resfinderDict, excluded_seqDict


def read_csv_DBase(csvFilename):

    csvDict = {}
    excluded_seqDict = {}
    dicList = []
    f = open(csvFilename, 'r')
    n= 0
    for line in f:
        line = line.strip()
        if line != '' and line.startswith('#') == False:
            if n == 0:
                n = 1
                header = line.split('\t')
            else:
                dic = dict(zip(header, line.split('\t')))
                dicList.append(dic)
    f.close()

    for n, sdic in enumerate(dicList):
        ddic = {'source':'','entry_name':'','alternative_names':'','nomenclature_trouble':'','mol_type':'cds',
                'dna_accession':'','blastn_evalue':'1e-50','dna_snp':'', 'dna_sequence':'',
                'prot_accession':'','blastp_evalue':'1e-50','prot_snp':'', 'prot_sequence':'',
                'function_grp_names':'','mechanism_names':'','operon_grp_name':'','cluster90_grp_name':'',
                'comment':'','taxonomy':'','reference':'','curated_by':'','db_names':''}
        keys = set(ddic.keys()).intersection(sdic.keys())
        if 'entry_name' not in keys or 'prot_sequence' not in keys:
            print '\nParsing csv file message:'
            print 'the column entry_name and prot_sequence are required'
            print 'Abort rec parsing:', sdic
            print ''
            excluded_seqDict[str(n)] = sdic
        else:
            for key in keys:
                ddic[key] = sdic[key]
            print 'entry_name: %s \tprot_sequence: %s' % (ddic['entry_name'], ddic['prot_sequence'][:20]+'...')
            csvDict[str(n)] = ddic
    return csvDict, excluded_seqDict


def update_cnrDict(cnrDict, source_dict):

    dDict = {}
    for n, skey in enumerate(source_dict.keys()):
        #sdna_sequence  = source_dict[skey]['dna_sequence']
        sprot_sequence = source_dict[skey]['prot_sequence']
        if dDict != {}:
            cnrDict = dDict
        dDict = cnrDict

        if dDict != {}:
            cnr_recordList = cnrDict.keys()
            cnr_recordList.sort()
            found = False
            for key in cnr_recordList:
                #dna_sequence = dDict[key]['dna_sequence']
                prot_sequence = dDict[key]['prot_sequence']

                if prot_sequence == sprot_sequence or sprot_sequence in prot_sequence:
                    print '\nSequence: %i' % (n+1)
                    print 'Protein sequence already in the database:'
                    print 'We are just checking the metadata:'
                    #excludedkeyList = ['blastn_evalue','dna_accession','dna_sequence','prot_sequence',
                    #                   'prot_accession','source','nomenclature_trouble','curated_by','taxonomy',
                    #                   'reference','cluster90_grp_name', 'operon_grp_name', 'function_grp_names',
                    #                   'mechanism_names']
                    check_keyList = ['entry_name','alternative_names','mol_type',
                                     'blastp_evalue','prot_snp', 'comment','curated_by']
                    dDict[key] = update_metadata(dDict[key], source_dict[skey], check_keyList)
                    found = True
                    break
                    # 'entry_name','blastp_evalue',
                elif prot_sequence in sprot_sequence:
                    print '\nSequence: %i' % (n+1)
                    print 'Overlapping protein sequence already in the database:'
                    print '\nSource:'
                    print_record(source_dict[skey])
                    print '\nDestination:'
                    print_record(dDict[key])
                    print ''
                    sbjt_seq = prot_sequence
                    qury_seq = sprot_sequence
                    aln_sbjt, aln_qury, idper, sbjt_cov, qury_cov = pairwisealn(sbjt_seq, qury_seq)
                    print 'In:\t%s' % aln_sbjt
                    print 'Nw:\t%s' % aln_qury
                    print 'Identity:', str(idper), '\tIn coverage:', str(sbjt_cov), '\tNw coverage:', str(qury_cov)
                    print ''
                    ipt = raw_input('Import the sequence as a new record Y[n]:')
                    if ipt == '' or ipt == 'n':
                        found = True
                        #excludedkeyList = ['nomenclature_trouble','curated_by','taxonomy',
                        #                   'reference','cluster90_grp_name','operon_grp_name','function_grp_names',
                        #                   'mechanism_names']

                        check_keyList = ['source','entry_name','alternative_names','mol_type',
                                         'dna_accession','blastn_evalue','dna_snp', 'dna_sequence',
                                         'prot_accession','blastp_evalue','prot_snp', 'prot_sequence','comment']

                        updatedRecord = update_metadata(dDict[key], source_dict[skey], check_keyList)
                        if updatedRecord['prot_sequence'] != prot_sequence:
                            del dDict[key]
                            time.sleep(2)
                            dDict[generate_DBkey()] = updatedRecord
                        else:
                            dDict[key] = updatedRecord

                    elif ipt == 'Y':
                        found == False
                        n = n - 1

            if found == False:
                print '\nSequence: %i' % (n+1)
                print 'New protein sequence:'
                print 'The database is updated with the record:\n'
                print_record(source_dict[skey])
                dDict = insert_record(dDict, source_dict, skey)

            #write_json_DBase(cnrDict, 'tmp_cnrDict.json')
            write_csv_DBase(dDict, 'tmp_cnrDict.csv')

        else:
            print '\nFirst sequence of the database:'
            print 'The database is created with the record:\n'
            print_record(source_dict[skey])
            time.sleep(2)
            dDict[generate_DBkey()] = source_dict[skey]
            #write_json_DBase(cnrDict, 'tmp_cnrDict.json')
            write_csv_DBase(dDict, 'tmp_cnrDict.csv')

    return dDict


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


def generate_Date():
    y = time.localtime(time.time()).tm_year
    m = time.localtime(time.time()).tm_mon
    j = time.localtime(time.time()).tm_mday
    date = '%i-%i-%i' % (j,m,y)
    return date


def update_metadata(ddict,sdict, check_keyList=[]):
    for key in ddict.keys():
        if sdict[key] != '' and sdict[key] != ddict[key]:
            if key == 'entry_name' and ddict['curated_by'] == '':
                if ddict[key].upper().replace('-','').replace('(','').replace(')','') == sdict[key].upper().replace('-','').replace('(','').replace(')','') \
                    or sdict[key].upper().replace('-','').replace('(','').replace(')','') in ddict['alternative_names'].upper().replace('-','').replace('(','').replace(')',''):
                    if 'resfinder' in ddict['source'] and 'CARD' in sdict['source']:
                        ddict['curated_by'] == 'CARD_RESFINDER_AGREEMENT'
                    elif ddict['source'] != sdict['source']:
                        ddict['curated_by'] == 'DOUBLE_CHECKED'
                elif ddict[key].upper().replace('-','').replace('(','').replace(')','') != sdict[key].upper().replace('-','').replace('(','').replace(')','') and \
                    sdict[key].upper().replace('-','').replace('(','').replace(')','') not in ddict['alternative_names'].upper().replace('-','').replace('(','').replace(')',''):
                    if ddict['nomenclature_trouble'] == '':
                        ddict['nomenclature_trouble'] = sdict['source'] + '_' + generate_Date() + ':' + sdict['entry_name']
                    else:
                        ddict['nomenclature_trouble'] += '|' + sdict['source'] + '_' + generate_Date() + ':' + sdict['entry_name']

            elif key == 'blastp_evalue':
                if ddict[key] == '':
                    ddict[key] = sdict[key]
                elif float(sdict[key]) < float(ddict[key]):
                    ddict[key] = sdict[key]

            elif key == 'blastn_evalue':
                if ddict[key] == '':
                    ddict[key] = sdict[key]
                elif float(sdict[key]) < float(ddict[key]):
                    ddict[key] = sdict[key]

            elif key == 'comment' and sdict[key] != ddict[key]:
                if ddict[key] == '':
                    ddict[key] = sdict[key]
                else:
                    try:
                        ddict[key] = ddict[key].encode('ascii','ignore') + '. ' +  sdict[key].encode('ascii','ignore')
                    except UnicodeEncodeError:
                        continue

            elif key in check_keyList:
                print '\nSource:'
                print_record(sdict)
                print '\nDestination:'
                print_record(ddict)
                print ''
                print 'Source:', sdict['entry_name'], '-> Destination:', ddict['entry_name']
                print 'The key: %s' % key
                print 'The current value: %s' % ddict[key],'\tThe new value    : %s' % sdict[key]
                print 'r: to replace\tm: to merge\tn: to do nothing'
                transfert = raw_input('Your choise r/m/[n] : ')
                if transfert == 'r':
                    ddict[key] = sdict[key]
                elif transfert == 'm':
                    if key == 'snp_protein' or key == 'snp_dna':
                        ddict[key] +=  ',' + sdict[key]
                    else:
                        ddict[key] += '|' + sdict[key]
                print '\nUpdated record:'
                print_record(ddict)
            else:
                print 'Identical metadata for %s' % key
    return ddict


def insert_record(cnrDict, source_dict, skey):

    time.sleep(2)
    newkey = generate_DBkey()
    cnrDict[newkey] = source_dict[skey]
    print '\nA new record have been insert:\nkey: %s\tentry_name: %s' % (newkey, cnrDict[newkey]['entry_name'])
    write_DB_faa(cnrDict, 'cnrDBase_insert_record.faa', {})
    cmd = '/usr/local/ncbi-blast-2.2.31+/bin/makeblastdb -in %s -dbtype %s -out %s' % ('cnrDBase_insert_record.faa', 'prot', 'cnrDBase_insert_record')
    os.system(cmd)
    print 'Protein blast database done!'

    query_seq = source_dict[skey]['prot_sequence']
    query_length = len(query_seq)
    query_name = '>' + skey + '_' + cnrDict[newkey]['entry_name']
    f = open('query.faa','w')
    f.write(query_name+'\n'+query_seq+'\n')
    f.close()
    print 'Query has been written'

    xmlblastFile = "query.blastpRes.xml"
    cmd = '/usr/local/ncbi-blast-2.2.31+/bin/blastp -out %s -outfmt %i -query %s  -db %s  -num_threads %i' % \
          (xmlblastFile,5,'query.faa','cnrDBase_insert_record', 6)
    os.system(cmd)
    print 'blastp done!\n'

    groupList = []
    result_handle = open(xmlblastFile, 'r')
    from Bio.Blast import NCBIXML
    blast_records = NCBIXML.parse(result_handle)
    for blast_record in blast_records:
        for alignment in blast_record.alignments:
            alignTitle = alignment.title
            spacepos = alignTitle.index(' ')
            hit_infoList = alignTitle[spacepos:].strip().split(' | ')
            hit_key = hit_infoList[0].split('_')[0]
            hit_name = hit_infoList[0].split('_')[1]
            hit_length = len(cnrDict[hit_key]['prot_sequence'])
            for hsp in alignment.hsps:
                blast_matchNb = hsp.identities
                blast_query = hsp.query.encode('ascii', 'replace')
                blast_hit = hsp.sbjct.encode('ascii', 'replace')
                idperc = blast_matchNb / float(len(blast_query.replace('-','')))
                query_cov = len(blast_query.replace('-','')) / float(query_length)
                hit_cov = len(blast_hit.replace('-','')) / float(hit_length)
                if idperc >= 0.9 and query_cov > 0.9 and hit_cov > 0.9:
                    print 'Query key: %s\t Query name: %s' % (newkey, cnrDict[newkey]['entry_name'])
                    print 'Hit key: %s\tHit name: %s' % (hit_key, hit_name)
                    print 'Identity percentage: %.2f' % (idperc)
                    print 'Query coverage: %.2f' % query_cov
                    print 'Hit coverage: %.2f' % hit_cov
                    print ''
                    groupList.append(hit_key)

    #=> collecter et affecter cluster90_grp_name, function_grp_names, mechanism_names, db_names
    cluster90_grp_name = []
    function_grp_names = []
    mechanism_names = []
    db_names = []
    for grp_key in groupList:
        function_grp_names += cnrDict[grp_key]['function_grp_names'].split('|')
        mechanism_names += cnrDict[grp_key]['mechanism_names'].split('|')
        db_names += cnrDict[grp_key]['db_names'].split('|')
        cluster90_grp_name.append(cnrDict[grp_key]['cluster90_grp_name'])

    cluster90_grp_name = list(set(cluster90_grp_name))[-1]
    if cluster90_grp_name == '':
        cluster90_grp_name = cnrDict[groupList[0]]['entry_name']
    function_grp_names = list(set(function_grp_names))
    function_grp_names.sort()
    function_grp_names = '|'.join(function_grp_names)
    mechanism_names = list(set(mechanism_names))
    mechanism_names.sort()
    mechanism_names = '|'.join(mechanism_names)
    db_names = list(set(db_names))
    db_names.sort()
    db_names = '|'.join(db_names)

    for grp_key in groupList:
        cnrDict[grp_key]['cluster90_grp_name'] = cluster90_grp_name
        cnrDict[grp_key]['function_grp_names'] = function_grp_names
        cnrDict[grp_key]['mechanism_names'] = mechanism_names
        cnrDict[grp_key]['db_names'] = db_names

    #os.system('rm query.blastpRes.xml query.faa cnrDBase_insert_record*')
    return cnrDict

# {key:{'entry_name': , 'alternative_names': , 'nomenclature_trouble': , 'comment: ', 'curated_by': ,
# 'source': , 'mol_type':'cds'/'dna',
# 'blastn_evalue': , 'dna_snp': , 'dna_seq': , 'dna_accession': , 'dna_length': ,
# 'blastp_evalue': , 'prot_snp': , 'prot_seq': , 'prot_accession': , 'prot_length': , 'operon_grp_name': ,
# 'reference': , 'cluster90_grp_name: ', 'function_grp_names': , 'mechanism_names': , 'taxonomy': code|species, 'db_names'}


def translate(recid, filename, seq):
    print 'Translation of: %s' % seq
    try:
        prot = seq.translate(table="Bacterial", cds=True, to_stop=True)
    except Bio.Data.CodonTable.TranslationError:
        seq = seq.reverse_complement()
        try:
            prot = seq.translate(table="Bacterial", cds=True, to_stop=True)
        except Bio.Data.CodonTable.TranslationError:
            print '\n!!! Translation problem with sequence %s in file %s !!!\n' % (recid, filename)
            prot = ''

    return seq, prot


def main(args):
    if args.verbose == '1':
        print "[info] logs saved to : " + str(working_directory) + "/app.log"
        logging.basicConfig(filename='app.log', level=logging.DEBUG,
                            format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                            datefmt='%Y-%m-%d %I:%M:%S %p')
        logging.info('main => start updateDBase')

    if os.path.exists('tmp_cnrDict.csv'):
        print 'Load temporary file tmp_cnrDict.csv'
        cnrDict = read_cnr_DBase('tmp_cnrDict.csv')
        for key in cnrDict.keys():
            #print cnrDict[key]
            print_record(cnrDict[key])
        print '\nTemporary file tmp_cnrDict.csv loaded as current cnr database\n'
        out_filename = args.outFile
        if out_filename == '':
            out_filename = 'outDB_new.csv'

    else:
        cnr_db_filename = args.db
        out_filename = args.outFile
        if os.path.exists(cnr_db_filename) == False and out_filename == '':
            cnr_db_filename = ''
            out_filename = 'outDB_new.csv'
        elif os.path.exists(cnr_db_filename) == True and out_filename == '':
            out_dir = os.path.dirname(cnr_db_filename)
            out_filename = os.path.splitext(os.path.basename(cnr_db_filename))[0].split('_')[0]
            version = int(os.path.splitext(os.path.basename(cnr_db_filename))[0].split('_')[1])+1
            out_filename += '_%i.csv' % version
        print '\nCNR database file: %s' % cnr_db_filename
        print 'Output file      : %s' % out_filename
        if cnr_db_filename != '':
            cnrDict = read_cnr_DBase(cnr_db_filename)
        else:
            cnrDict = {}

    if os.path.exists(args.resfinderDB) == True:
        resfinderList = glob.glob(os.path.join(args.resfinderDB, '*.fsa'))
        resfinderNote = os.path.join(args.resfinderDB, 'notes.txt')
        if os.path.exists(resfinderNote) == False:
            resfinderNote = ''
    else:
        resfinderNote = ''
        resfinderList = []
    resfinderList.sort()
    print '\nResfinder file list:', resfinderList
    print 'Resfinder note file:', resfinderNote

    if resfinderList != []:
        resfinderDict, excluded_resfinderDict = read_resfinder_DBase(resfinderList, resfinderNote)
        write_csv_DBase(excluded_resfinderDict, 'resfinder_unparsed_sequence.csv')
        write_csv_DBase(resfinderDict, 'resfinder_parsed_sequence.csv')
    else:
        resfinderDict = {}
        excluded_resfinderDict = {}

    print '\nNumber resfinder parsed record: %i' % len(resfinderDict.keys())
    print 'If any, parsed resfinder sequences have been write in resfinder_parsed_sequence.csv\n'
    print 'Number resfinder unparsed record: %i' % len(excluded_resfinderDict.keys())
    print 'If any, unparsed resfinder sequences have been write in resfinder_unparsed_sequence.csv\n'
    time.sleep(5)

    if resfinderDict != {}:
        print '\nUpdate with resfinder:'
        cnrDict = update_cnrDict(cnrDict, resfinderDict)
        write_csv_DBase(cnrDict, out_filename)
    else:
        print 'No record from resfinder to update the database'

    if os.path.exists(args.cardDB) == True:
        cardList = glob.glob(os.path.join(args.cardDB, 'card.json'))
    else:
        cardList = []
    cardList.sort()
    print '\nCARD file list:', cardList

    if cardList != []:
        cardDict, excluded_cardDict = read_card_DBase(cardList[0])
        write_csv_DBase(excluded_cardDict, 'card_unparsed_sequence.csv')
        write_csv_DBase(cardDict, 'card_parsed_sequence.csv')
    else:
        cardDict = {}
        excluded_cardDict = {}
    print '\nNumber of parsed card record: %i' % len(cardDict.keys())
    print 'If any, parsed card sequences have been write in card_parsed_sequence.csv\n'
    print 'Number of unparsed card record: %i' % len(excluded_cardDict.keys())
    print 'If any, unparsed card sequences have been write in card_unparsed_sequence.csv\n'
    time.sleep(5)

    if cardDict != {}:
        print '\nUpdate with CARD:'
        cnrDict = update_cnrDict(cnrDict, cardDict)
        write_csv_DBase(cnrDict, out_filename)
    else:
        print 'No record from CARD to update the database'

    if os.path.exists(args.csvDB) == True:
        csvFile = args.csvDB
    else:
        csvFile = ''
    print '\nCSV file:', csvFile

    if csvFile != '':
        csvDict, excluded_csvDict = read_csv_DBase(csvFile)
        write_csv_DBase(excluded_csvDict, 'csv_unparsed_sequence.csv')
        write_csv_DBase(csvDict, 'csv_parsed_sequence.csv')
    else:
        csvDict = {}
        excluded_csvDict = {}

    print '\nNumber of parsed csv record: %i' % len(csvDict.keys())
    print 'If any, parsed csv sequences have been write in csv_parsed_sequence.csv\n'
    print 'Number of unparsed csv record: %i' % len(excluded_csvDict.keys())
    print 'If any, unparsed csv sequences have been write in csv_unparsed_sequence.csv\n'
    time.sleep(5)

    if csvDict != {}:
        print '\nUpdate with CSV file:'
        cnrDict = update_cnrDict(cnrDict, csvDict)
        write_csv_DBase(cnrDict, out_filename)
    else:
        print 'No record from CSV file to update the database'


def version():
    return "1.0"


def run():
    parser = argparse.ArgumentParser(description='Resistance Gene Identifier - Version ' + version())
    parser.add_argument('-db', '--dbase', dest="db", default="cnrDB_5.csv", help="The CNR database in csv format (default=cnr_db_5.csv)")
    parser.add_argument('-card', '--cardDB', dest="cardDB", default="card", help='must be the CARD directory (default=card)')
    parser.add_argument('-resf', '--resfinderDB', dest="resfinderDB", default="resfinder", help="must be the resfinder directory (default=resfinder)")
    parser.add_argument('-csv', '--csvFile', dest="csvDB", default="db_arm_5.csv", help="must be a csv file with formated header line (default=db_arm_5.csv)")
    parser.add_argument('-o', '--outFile', dest="outFile", default="", help="must be the ouput file name")
    parser.add_argument('-U', '--update', dest="update", default='0', help="'1' to import and update the database or '0' to import and write the data as csv file (default=0)")
    parser.add_argument('-v', '--verbose', dest="verbose", default="0", help="log process to file. Options are 0 or 1  (default = 0 for no logging)")
    parser.add_argument('-V', '--version', action='version', version='runBlast-' + version(), help="Prints version number")
    args = parser.parse_args()
    main(args)


if __name__ == '__main__':
    run()


