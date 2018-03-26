#!/usr/bin/python
# Write by Richard Bonnet
# Date: 19/01/2016
import fnmatch
import os
import argparse
import re
import time


install_dir = os.path.dirname(os.path.realpath(__file__))


def recursive_fileList(directory, filter):
    fileList = []
    for root, dirnames, filenames in os.walk(directory):
        for filename in fnmatch.filter(filenames, filter):
            fileList.append(os.path.join(root, filename))
    return fileList


def annotGBK(gbkFile, csvFile, idperc, cov, outdir):
    cmd = '%s -gbk %s -csv %s -id %.2f -cv %.2f -o %s' % \
          (os.path.join(install_dir,'annotGBK.py'), gbkFile, csvFile, float(idperc), float(cov), outdir)
    print '\nEdition of gbk file %s' % gbkFile
    os.system(cmd)
    print 'Edition of gbk file %s done!\n' % gbkFile


def read_sample_file(sampleFile):
    sampleDict = {}
    f = open(sampleFile,'r')
    for line in f:
        line = line.strip()
        if line != '' and line.startswith('#') == False:
            sample_id, species = line.split('\t')
            sampleDict[sample_id] = species.lower()
    return sampleDict


def read_mlst_txt(mlst_txt):
    mlstDict = {}
    f = open(mlst_txt,'r')
    for line in f:
        line = line.strip()
        if line != '' and line.startswith('#') == False:
            line = line.split('\t')
            mlstDict[line[0]] = line[1:]
    return mlstDict


def read_mlst_scheme(mlst_csv):
    mlst_name = os.path.splitext(os.path.basename(mlst_csv))[0].lower()
    mlst_gene_Dict = {}
    pattern_st_Dict = {}
    f = open(mlst_csv,'r')
    n = 0
    for line in f:
        line = line.strip().split('\t')
        if n == 0:
            mlst_gene_Dict[mlst_name.lower()]=line[1:]
            for n, item in enumerate(mlst_gene_Dict[mlst_name.lower()]):
                mlst_gene_Dict[mlst_name.lower()][n] = mlst_gene_Dict[mlst_name.lower()][n].lower()
            n = 1
        else:
            pattern = ''
            for item in line[1:]:
                pattern = pattern + '|' + item
            pattern = pattern[1:]
            pattern_st_Dict[pattern] = line[0]
    return mlst_gene_Dict, pattern_st_Dict


def parse_gene_id(recid):
    pattern = re.compile('([A-Za-z_0-9]+)[-]?([0-9]+)')
    match = pattern.match(recid)
    gene_name = match.group(1)
    allele_nb = match.group(2)
    return gene_name, allele_nb


def read_mlst_fasta(mlst_fasta, mlst_gene_list):
    fasta_file_list = []
    gene_length_dict = {}
    from Bio import SeqIO
    for gene_name in mlst_gene_list:
        gene_length_dict[gene_name.lower()]={}
        filename = os.path.join(os.path.dirname(mlst_fasta),gene_name.lower()+'.fasta')
        df = open(filename,'w')
        print 'MLST: %s\tGENE: %s' % (os.path.splitext(os.path.basename(mlst_fasta))[0], filename)
        sf = open(mlst_fasta,'r')
        for rec in SeqIO.parse(sf,'fasta'):
            gene, allele_nb = parse_gene_id(rec.id)
            if gene.lower() == gene_name.lower():
                #print rec
                gene_length_dict[gene_name.lower()][rec.id.lower()]=len(rec.seq)
                seq_name = '>'+rec.id.lower()
                seq      = str(rec.seq)
                #print seq_name
                #print seq
                df.write(seq_name+'\n')
                df.write(seq+'\n')
        fasta_file_list.append(filename)
        df.close()
    print ''
    sf.close()
    return fasta_file_list, gene_length_dict


def makeblastDB(blastDirectory, fasta_file_list, dbtype):
    for fasta_input_name in fasta_file_list:
        cmd = '%s -in %s -dbtype %s -out %s' % (os.path.join(blastDirectory, 'makeblastdb'), fasta_input_name, dbtype, fasta_input_name)
        os.system(cmd)
        print 'Blastn database %s done!\n' % fasta_input_name
    return fasta_file_list


def runBlast(blastDirectory, blastdb_file, query_fasta_file, threads, blast_evalue, out_type, xmlblastResults):
    cmd = '%s -out %s -outfmt %s -query %s -db %s -num_threads %s -evalue %s' % \
          (os.path.join(blastDirectory, 'blastn'), xmlblastResults, out_type,
           query_fasta_file, blastdb_file, threads, blast_evalue)
    #print '\nblastn is started...'
    os.system(cmd)
    #print cmd
    #print 'Blast result file %s written\n' % xmlblastResults


def parsecsvBlast(csvblastResults,gene_length_dict ):
    header = ['query_id', 'subject_id', 'Perc_identity',
              'alignment_length', 'mismatches', 'gap_opens',
              'q_start', 'q_end', 's_start', 's_end',
              'evalue', 'bit score', 'sequence']
    result_dict = {}
    perc_pass = 0.0
    allele_ln = 0
    for line in open(csvblastResults):
        data = dict(zip(header, line.strip().split('\t')))
        subject_id = data['subject_id']
        gene, allele_nb = parse_gene_id(subject_id)
        allele_length = gene_length_dict[gene][subject_id]
        perc_identity = round((int(data['alignment_length'])-int(data['mismatches'])-int(data['gap_opens']))/float(allele_length)*100,2)
        if perc_identity > perc_pass :
            perc_pass = perc_identity
            data['Perc_identity'] = perc_pass
            result_dict = data
            if perc_identity == 100.00:
                data['allele_nb'] = allele_nb
                break
            else:
                data['allele_nb'] = allele_nb + '-like [%s : %s..%s]' % (result_dict['query_id'],result_dict['q_start'],result_dict['q_end'])
                result_dict['subject_id'] = result_dict['subject_id'] + '-like [%s : %s..%s]' % (result_dict['query_id'],result_dict['q_start'],result_dict['q_end'])
                if perc_pass > 100:
                    perc_pass = 99.9
    if result_dict == {}:
        result_dict = {'query_id':'NA', 'subject_id':'NA', 'Perc_identity':'NA',
                       'alignment_length':'NA', 'mismatches':'NA', 'gap_opens':'NA',
                       'q_start':'NA', 'q_end':'NA', 's_start':'NA', 's_end':'NA',
                       'evalue':'NA', 'bit score':'NA', 'sequence':'NA','allele_nb':'NA'}
    print 'Blast result: %s\tQuery: %s\tstart: %s\tend: %s\tPerc. identity: %s' % \
          (result_dict['subject_id'],result_dict['query_id'],result_dict['q_start'],result_dict['q_end'],str(result_dict['Perc_identity']))
    cmd = 'rm %s' % csvblastResults
    os.system(cmd)
    return result_dict


def get_ST(res_mlst_dict, mlst_gene_dict, mlst_name, pattern_st_dict):
    pattern = ''
    for gene in mlst_gene_dict[mlst_name]:
        pattern = pattern + '|' + res_mlst_dict[gene]['allele_nb']
    pattern = pattern[1:]
    try:
        ST = pattern_st_dict[pattern]
    except KeyError:
        ST = '?'
    return ST, pattern


def csvWrite(species, mlst_name, outdir, mlst_gene_dict, sample, ST, pattern):
    outfile = os.path.join(outdir, species + '__' + mlst_name + '__MLST.csv')
    print 'Write MLST results in %s' % outfile
    if os.path.exists(outfile) == False:
        f = open(outfile, 'w')
        f.write('sample_id\tST\t%s\n' % ('\t'.join(mlst_gene_dict[mlst_name.lower()])))
    else:
        f = open(outfile, 'a')
    f.write('%s\t%s\t%s\n' % (sample, ST, pattern.replace('|','\t')))
    f.close()


def main(args):

    sampleID = args.sampleID
    fasFile = args.fasFile

    outdir = args.outDir
    if outdir == '':
        outdir = os.path.dirname(fasFile)

    if sampleID == '':
        sampleID = os.path.splitext(os.path.basename(fasFile))[0].split('_')[0].split('.')[0]
    wkdir = os.path.dirname(fasFile)

    sampleFile = args.sampleFile
    if sampleFile == '':
        sampleFile = os.path.join(wkdir, 'sample.csv')
    if os.path.exists(sampleFile) == True:
        print '\nThe sample file %s have been found.\n' % sampleFile
        sample_species_Dict = read_sample_file(sampleFile)
        sampleList = sample_species_Dict.keys()
        sampleList.sort()
        for key in sampleList:
            print key, sample_species_Dict[key]
        print ''
    else:
        print '\nNo sample file.\nNo MLST can be assignated!\n'
        exit()

    mlstDir = args.mlstDir
    blastDirectory = args.blastDirectory
    dbtype = 'nucl'

    species_mlst_dict = read_mlst_txt(os.path.join(mlstDir,'mlst.txt'))
    mlstList = species_mlst_dict.keys()
    mlstList.sort()
    for key in mlstList:
        print key, species_mlst_dict[key]
    print ''

    species = sample_species_Dict[sampleID].lower()
    for mlst_name in species_mlst_dict[species]:
        res_mlst_dict = {}
        mlst_gene_dict, pattern_st_dict   = read_mlst_scheme(os.path.join(mlstDir, mlst_name+'.csv'))
        fasta_file_list,gene_length_dict = read_mlst_fasta(os.path.join(mlstDir, mlst_name+'.fasta'), mlst_gene_dict[mlst_name.lower()])
        blastdb_file_list = makeblastDB(blastDirectory, fasta_file_list, dbtype)
        for blastdb_file in blastdb_file_list:
            gene_name = os.path.splitext(os.path.basename(blastdb_file))[0]
            csvblastResults = fasFile + '.blastResults.csv'
            runBlast(blastDirectory, blastdb_file, fasFile, '8', '1e-1', '6', csvblastResults)
            result = parsecsvBlast(csvblastResults,gene_length_dict )
            res_mlst_dict[gene_name] = result
        for fasta_file in fasta_file_list:
            cmd = 'rm %s*' % fasta_file
            os.system(cmd)
        ST, pattern = get_ST(res_mlst_dict, mlst_gene_dict, mlst_name.lower(), pattern_st_dict)
        print '\nMLST name: %s\tSequence type: %s\tPattern: %s\n' % (mlst_name.lower(), ST, pattern)
        time.sleep(3)
        csvWrite(species.replace(' ','_'), mlst_name, outdir, mlst_gene_dict, sampleID, ST, pattern)


def version():
    return "1.0"


def run():
    parser = argparse.ArgumentParser(description='seqMLST - Version ' + version())
    parser.add_argument('-f','--fasFile', dest="fasFile", help='The fasta file')
    parser.add_argument('-o','--outDir', dest="outDir", default='', help='The fasta file')
    parser.add_argument('-id','--sampleID', dest="sampleID", default='', help='The fasta file')
    parser.add_argument('-s','--sampleFile', dest="sampleFile", default='', help="The tab sample file with corresponding bacterial species")
    parser.add_argument('-mlst','--mlstDir', dest="mlstDir", default='/usr/local/seqblast/db_mlst', help="A directory for mlst schemes and data (default: /usr/local/seqblast/mlst)")
    parser.add_argument('-bd','--blastDirectory', dest="blastDirectory", default="/usr/local/ncbi-blast-2.6.0+/bin", help="Path to blast directory (default=/usr/local/ncbi-blast-2.6.0+/bin)")
    parser.add_argument('-e','--evalue',  dest="evalue", default='1e-10', help="blast evalue cutoff (default: 1e-30)")
    parser.add_argument('-t','--threads', dest="threads", default='8', help="Thread number for blast (default: 8)")
    parser.add_argument('-V','--version', action='version', version='seqMLST-' + version(), help="Prints version number")
    args = parser.parse_args()
    main(args)


if __name__ == '__main__':
    run()
