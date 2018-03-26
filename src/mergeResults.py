#!/usr/bin/python
# Write by Richard Bonnet
# Date: 20/01/2016
import fnmatch
import os
import argparse
import logging

def loadFile(directory):
    fileList = []
    for root, dirnames, filenames in os.walk(directory):
        for filename in fnmatch.filter(filenames, '*.blastRes.xml.csv'):
            fileList.append(os.path.join(root, filename))
    return fileList


def loadResult(fileList, idperc, qcov, hcov):
    fresDict = {}
    mresDict = {}
    fhitDict  = {}
    mhitDict  = {}
    db_update_list = []
    seq_list = []
    for filename in fileList:
        print '\nMerge data from %s' % filename
        key = os.path.basename(filename).split('.')[0]
        fresDict[key]={}
        mresDict[key]={}
        f = open(filename, 'r')
        for n, line in enumerate(f):
            line = line.strip()
            if n == 0:
                header = line.split('\t')
            else:
                data = line.split('\t')
                dic = dict(zip(header, data))
                fCluster = dic['Functional cluster'].replace('|','-')
                if fCluster.strip() == '':
                    fCluster = 'Undefined-Cluster'
                mClusters = dic['Mechanism']
                if mClusters.strip() == '' or mClusters.strip() == '|':
                    mClusters = 'Undefined-Cluster'
                iCluster = dic['90 cluster']
                n = 0
                if float(dic['Identity (%)']) >= idperc and float(dic['Hit coverage (%)']) >= hcov: #and float(dic['Query coverage (%)']) > qcov
                    if dic['SNP(s) searched'] == '':
                        print 'Record: %s\tQuery: %s\tHit name: %s\tHit key: %s\tIdentity: %.1f\tQuery Coverage: %.1f\tHit Coverage: %.1f' % \
                              (key, dic['Query name'], dic['Hit name'], dic['Hit key in dtbase'],
                               float(dic['Identity (%)']), float(dic['Query coverage (%)']), float(dic['Hit coverage (%)']))
                    elif (dic['SNP(s) searched'] != '' and (dic['SNP(s) found'] != '' or dic['Other SNP(s)'] != '')):
                        print 'Record: %s\tQuery: %s\tHit name: %s\tHit key: %s\tIdentity: %.1f\tQuery Coverage: %.1f\tHit Coverage: %.1f\tSNP: %s\tSUB: %s' % \
                              (key, dic['Query name'], dic['Hit name'], dic['Hit key in dtbase'],
                               float(dic['Identity (%)']), float(dic['Query coverage (%)']), float(dic['Hit coverage (%)']),
                               dic['SNP(s) found'], dic['Other SNP(s)'])

                    #FUNCTIONAL CLASSIFICATION
                    try:
                        fresDict[key][fCluster].append(dic)
                    except KeyError:
                        fresDict[key][fCluster]=[dic]

                    try:
                        fhitDict[fCluster].append(iCluster+'__'+dic['Hit name']+'__'+dic['Hit key in dtbase'])
                    except KeyError:
                        fhitDict[fCluster] = [iCluster + '__' + dic['Hit name'] + '__' + dic['Hit key in dtbase']]
                    fhitDict[fCluster] = list(set(fhitDict[fCluster]))
                    fhitDict[fCluster].sort()

                    #MECHANISM CLASSIFICATION
                    for mCluster in mClusters.split('|'):
                        if mCluster != '':
                            try:
                                mresDict[key][mCluster].append(dic)
                            except KeyError:
                                mresDict[key][mCluster] = [dic]

                            try:
                                mhitDict[mCluster].append(iCluster + '__' + dic['Hit name'] + '__' + dic['Hit key in dtbase'])
                            except KeyError:
                                mhitDict[mCluster] = [iCluster + '__' + dic['Hit name'] + '__' + dic['Hit key in dtbase']]
                            mhitDict[mCluster]=list(set(mhitDict[mCluster]))
                            mhitDict[mCluster].sort()

                    #SEQUENCE FOR UPDATING DATABASE
                    n += 1
                    upkey = n
                    entry_name = dic['Hit name']
                    comment = 'Imported from CNR data'
                    source = 'CNR_%s' % dic['Query name']
                    mol_type = dic['mol_type']

                    blastn_evalue = dic['blastn_evalue']
                    try:
                        dna_sequence = dic['Query dna sequence']
                    except KeyError:
                        dna_sequence = ''
                    dna_snp = dic['dna_snp']

                    blastp_evalue = dic['blastp_evalue']
                    try:
                        prot_sequence = dic['Query prot sequence']
                    except KeyError:
                        prot_sequence = ''
                    prot_snp = dic['prot_snp']

                    cluster90_grp_name = dic['90 cluster']
                    operon_grp_name = dic['Operon name']
                    function_grp_names = dic['Functional cluster']
                    mechanism_names = dic['Mechanism']
                    taxonomy = ''
                    db_names = dic['db_names']
                    dic ={'key':upkey, 'entry_name':entry_name, 'source': source, 'comment':comment, 'mol_type':mol_type,
                          'blastn_evalue':blastn_evalue,'dna_sequence':dna_sequence,'dna_snp':dna_snp,
                          'blastp_evalue':blastp_evalue,'prot_sequence':prot_sequence,'prot_snp':prot_snp,
                          'cluster90_grp_name':cluster90_grp_name,'operon_grp_name':operon_grp_name,
                          'function_grp_names':function_grp_names,'mechanism_names':mechanism_names,'taxonomy':taxonomy, 'db_names':db_names}
                    #if prot_sequence != '' and prot_sequence not in seq_list:
                    #    db_update_list.append(dic)
                    #    seq_list.append(prot_sequence)
                    db_update_list.append(dic)

    return [(fresDict, fhitDict), (mresDict, mhitDict)], db_update_list


def writeResult(resList, outDir, classType):
    resDict = resList[0]
    hitDict = resList[1]

    f_filename = os.path.join(outDir, 'merged_data_' + classType + '.csv')
    f3_filename = os.path.join(outDir, 'merged_data_' + classType + '_count.csv')

    n = 0
    f = open(f_filename, 'w')
    f3 = open(f3_filename, 'w')

    sample_ids = resDict.keys()
    sample_ids.sort()
    for sample_id in sample_ids:

        header = 'Sample_ID'
        txt = sample_id
        header3 = 'Sample_ID'
        txt3 = sample_id

        res_sample_Dict = resDict[sample_id]
        classList = hitDict.keys()
        classList.sort()
        for cl in classList:

            #f2_filename = os.path.join(outDir, 'merged_data_' + classType + '_' + cl + '.csv')
            #header2 = 'sample_ID'
            #txt2 = sample_id
            #if os.path.exists(f2_filename) == True:
            #    p = 1
            #    f2 = open(f2_filename, 'a')
            #else:
            #    p = 0
            #    f2 = open(f2_filename, 'w')

            try:
                resList = res_sample_Dict[cl]
                cl_nb = len(resList)
            except KeyError:
                resList = [{'Hit key in dtbase':''}]
                cl_nb = 0

            header3 = header3 + '\t' + cl
            txt3 = txt3 + '\t' + str(cl_nb)

            for item in hitDict[cl]:
                iCluster, hit_name, hit_key = item.split('__')
                header = header + '\t' + '%s: %s|%s' % (cl, hit_name, hit_key)
                #header2 = header2 + '\t' + '%s: %s|%s' % (cl, hit_name, hit_key)
                found = False
                for res in resList:
                    if res['Hit key in dtbase'] == hit_key:
                        found = True
                        break
                if found == False:
                    txt = txt + '\t' + 'NAN'
                #    txt2 = txt2 + '\t' + 'NAN'
                else:
                    if res['Match type'] == 'Perfect':
                        txt = txt + '\t' + '%s [Query: %s, id:100, cov: 100]' % (res['Hit name'], res['Query name'])
                #        txt2 = txt2 + '\t' + '%s [Query: %s, id:100, cov: 100]' % (res['Hit name'], res['Query name'])
                    else:
                        idp = round(float(res['Identity (%)']), 1)
                        cov = round((float(res['Query coverage (%)'])+float(res['Hit coverage (%)']))/2, 1)
                        snp = res['SNP(s) found']
                        sub = res['Other SNP(s)']
                        txt = txt + '\t' + '%s [Query: %s, id: %.1f, cov: %.1f, snp: %s, sub: %s]' % (res['Hit name'], res['Query name'], idp, cov, snp, sub)
                #        txt2 = txt2 + '\t' + '%s [Query: %s, id: %.1f, cov: %.1f, snp: %s, sub: %s]' % (res['Hit name'], res['Query name'], idp, cov, snp, sub)
            #if p == 0:
            #    f2.write(header2 + '\n')
            #f2.write(txt2 + '\n')
            #f2.close()

        if n == 0:
            f.write(header + '\n')
            f3.write(header3 + '\n')
            n = 1

        f.write(txt + '\n')
        f3.write(txt3 + '\n')

    f.close()
    f3.close()


def writeUpdate(updateList, outDir):
    f = open(os.path.join(outDir, 'update_db.csv'),'w')
    itemList = ['key','entry_name','source','comment','mol_type','blastp_evalue','prot_snp','prot_sequence',
                'blastn_evalue','dna_snp','dna_sequence','cluster90_grp_name','operon_grp_name','function_grp_names',
                'mechanism_names','taxonomy']
    f.write('\t'.join(itemList)+'\n')
    for rec in updateList:
        txt = ''
        for item in itemList:
            txt = txt + '\t' + str(rec[item])
        f.write(txt[1:]+'\n')
    f.close()


def main(args):
    outDir = args.outDir
    if args.verbose == '1':
        print "[info] logs saved to : " + str(outDir) + "/app.log"
        logging.basicConfig(filename='app.log', level=logging.DEBUG,
                            format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                            datefmt='%Y-%m-%d %I:%M:%S %p')
        logging.info('main => start mergeResults')

    if args.inpDir == None:
        print "[error] missing input(s)"
        print "[info] Try: mergeResults -h"
        logging.error("main => missing input(s)")
        logging.error("main => Try: mergeResults -h")
        exit()

    fileList = loadFile(args.inpDir)
    idperc = float(args.idPerc)
    hcov = float(args.Cov)
    qcov = float(args.Cov)
    print '\nMerge the results using %.2f identity cutoff and %.2f coverage cutoff...' % (idperc, hcov)
    resList, db_update_list = loadResult(fileList, idperc, qcov, hcov)
    writeResult(resList[0], outDir, 'function')
    writeResult(resList[1], outDir, 'mechanism')
    writeUpdate(db_update_list, outDir)


def version():
    return "1.0"


def run():
    inpDir = '/home/bacteriologie/assemblage/workdir_Run9_res'
    parser = argparse.ArgumentParser(description='mergeResults - Version ' + version())
    parser.add_argument('-d', '--inpDir', dest="inpDir", default=inpDir, help='must be the imput directory (default: %s)' % inpDir)
    parser.add_argument('-id','--idPerc', dest="idPerc", default="80", help="Identity percentage cutOff (default: 80)")
    parser.add_argument('-cv','--Cov', dest="Cov", default="80", help="Coverage cutoff (default: 80)")
    parser.add_argument('-o','--outDir', dest="outDir", default=inpDir, help="Coverage cutoff (default: %s)" % inpDir)
    parser.add_argument('-v', '--verbose', dest="verbose", default="0", help="log process to file. Options are 0 or 1  (default = 0 for no logging)")
    parser.add_argument('-V', '--version', action='version', version='mergeResults-' + version(), help="Prints version number")
    args = parser.parse_args()
    main(args)


if __name__ == '__main__':
    run()
