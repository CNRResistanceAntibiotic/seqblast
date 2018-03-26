#!/usr/bin/python
# Write by Richard Bonnet
# Date: 20/12/2015

import os, argparse, logging
from Bio import SeqIO


def loadData(dataFile):
    resDict = {}
    if os.path.exists(dataFile) == True:
        n = 0
        f = open(dataFile, 'r')
        for line in f:
            line = line.strip()
            if n == 0:
                keys = line.split('\t')
                n = 1
            else:
                values = line.split('\t')
                dic = dict(zip(keys, values))
                try:
                    resDict[dic['Query name']].append(dic)
                except KeyError:
                    resDict[dic['Query name']] = [dic]
        f.close()
    else:
        print '\nNo csv blast result file %s for annotation\n' % dataFile
    return resDict


def read_sequence_file(filename):
    if os.path.splitext(filename)[1] in ['.gbk', '.gb', '.gbf']:
        format = 'genbank'
    elif os.path.splitext(filename)[1] == '.faa':
        format = 'fasta'
    else:
        print '\nUnknown file format: %s' % os.path.splitext(filename)[1]
        print 'Please use .gbk, .gbf, .gb or .faa extention\n'
        exit()

    rec_lst = []
    #print filename, format
    for record in SeqIO.parse(open(filename, 'r'), format):
        rec_lst.append(record)
        #print record.id, record.description

    print "%i records from %s\n" % (len(rec_lst), filename)
    return rec_lst, format


def write_sequence_file(record_lst, outdir, prefix, format):
    """
    Write sequence record in different format
    :param lst: record list
    :param outfile: output file name
    :param format: fasta or genbank
    :return: output file name
    """
    if format == 'genbank':
        ext = '.gbk'
    else:
        ext = '.%s' % format

    outfile = os.path.join(outdir, prefix + ext)
    SeqIO.write(record_lst, outfile, format)
    return outfile


def editGbk(gbkFile, resDict, idperc_cutoff, hcov_cutoff, qcov_cutoff, outdir):
    sample = os.path.splitext(os.path.basename(gbkFile))[0]
    fasprefix = os.path.join(os.path.dirname(gbkFile), sample)
    faa = open(fasprefix+'_seqD_faa.fasta','w')
    fna = open(fasprefix+'_seqD_fna.fasta','w')
    indir = os.path.dirname(gbkFile)
    outprefix = sample + '_seqD'
    if os.path.exists(gbkFile) == True and resDict != {}:
        record_lst, format = read_sequence_file(gbkFile)
        for rec in record_lst:
            contig_name = rec.id
            for feature in rec.features:
                if feature.type == 'CDS':
                    for resname in resDict.keys():
                        gbklocus_tag = feature.qualifiers['locus_tag'][0]
                        if gbklocus_tag == resname:
                            #print 'feature', feature.qualifiers, resname
                            for dic in resDict[resname]:
                                hitname = dic['Hit name']
                                dbkey = dic['Hit key in dtbase']
                                function = dic['Functional cluster']
                                mechanism = dic['Mechanism']
                                idperc = round(float(dic['Identity (%)']),2)
                                qcov = float(dic['Query coverage (%)'])
                                hcov = float(dic['Hit coverage (%)'])
                                cov = round((qcov + hcov) / 2.0, 2)

                                snp = dic['SNP(s) found']
                                if snp == '': snp = 'None'
                                sub = dic['Other SNP(s)']
                                if sub == '': sub = 'None'

                                if dic['SNP(s) searched'] != '':
                                    txt = '%s [id: %.2f, cov: %.2f, Antibiotic resistance-associated SNP(s): %s]' % (hitname, idperc, cov, snp)
                                else:
                                    txt = '%s [id: %.2f, cov: %.2f]' % (hitname, idperc, cov)

                                if idperc >= idperc_cutoff and hcov >= hcov_cutoff:
                                    feature.qualifiers['inference'].append('ab initio CNR prediction by similarity')
                                    feature.qualifiers['product'].append(hitname)
                                    try:
                                        feature.qualifiers['note'].append(txt)
                                    except KeyError:
                                        feature.qualifiers['note'] = [txt]
                                    print 'Inference: ab initio prediction by similarity'
                                    print 'Feature: %s' % resname
                                    print 'Product:', feature.qualifiers['product']
                                    print 'Note:', feature.qualifiers['note']
                                    print ''

                                    locus_start = feature.location.nofuzzy_start
                                    locus_end = feature.location.nofuzzy_end
                                    if feature.strand == -1:
                                        locus_dna_seq = rec.seq[locus_start:locus_end].reverse_complement()
                                    else:
                                        locus_dna_seq = rec.seq[locus_start:locus_end]
                                    locus_prot_seq = locus_dna_seq.translate(table='Bacterial', cds=True)
                                    if dic['SNP(s) searched'] == '':
                                        locus_fas_id  = '>%s__%s__%s__%s__%s func:%s,mechanism:%s,id:%.2f,cov:%.2f' % \
                                                    (sample,contig_name,gbklocus_tag,hitname.replace(' ','_'),dbkey,
                                                     function,mechanism,idperc, cov)
                                    else:
                                         locus_fas_id  = '>%s__%s__%s__%s__%s func:%s,mechanism:%s,id:%.2f,cov:%.2f,snp:%s,sub:%s' % \
                                                    (sample,contig_name,gbklocus_tag,hitname.replace(' ','_'),dbkey,
                                                     function,mechanism,idperc, cov, snp.replace(',','|'), sub.replace(',','|'))

                                    faa_txt = str(locus_fas_id+'\n'+locus_prot_seq+'\n')
                                    #print faa_txt
                                    fna_txt = str(locus_fas_id+'\n'+locus_dna_seq+'\n')
                                    #print fna_txt
                                    faa.write(faa_txt)
                                    fna.write(fna_txt)
        faa.close()
        fna.close()
        outfile = write_sequence_file(record_lst, os.path.dirname(gbkFile), outprefix, 'genbank')
        print 'Genbank annotated file %s written' % outfile


def main(args):
    outdir = args.outDir
    if outdir == '':
        outdir = os.path.dirname(args.gbkFile)
    if args.verbose == '1':
        print "[info] logs saved to : " + str(outdir) + "/app.log"
        logging.basicConfig(filename='app.log', level=logging.DEBUG,
                            format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                            datefmt='%Y-%m-%d %I:%M:%S %p')
        logging.info('main => start annotGBK')

    if args.csvFile == None:
        print "[error] missing input(s)"
        print "[info] Try: runBlast -h"
        logging.error("main => missing input(s)")
        logging.error("main => Try: annotGBK -h")
        exit()

    if args.gbkFile == None:
        print "[error] missing input(s)"
        print "[info] Try: runBlast -h"
        logging.error("main => missing input(s)")
        logging.error("main => Try: annotGBK -h")
        exit()

    idperc = float(args.idPerc)
    hcov = float(args.Cov)
    qcov = float(args.Cov)

    dataFile = args.csvFile
    gbkFile = args.gbkFile
    resDict = loadData(dataFile)
    editGbk(gbkFile, resDict, idperc, hcov, qcov, outdir)


def version():
    return "1.0"


def run():
    parser = argparse.ArgumentParser(description='annotGBK - Version ' + version())
    parser.add_argument('-csv', '--csvFile', dest="csvFile", help='must be the csv blast result file')
    parser.add_argument('-gbk', '--gbkFile', dest="gbkFile", help='must be the corresponding gbk file')
    parser.add_argument('-id','--idPerc', dest="idPerc", default="80", help="Identity percentage cutOff (default: 0)")
    parser.add_argument('-cv','--Cov', dest="Cov", default="80", help="Coverage cutoff (default: 0)")
    parser.add_argument('-o','--outDir', dest="outDir", default='', help="Blast type. Options are blastp or blastn")
    parser.add_argument('-v', '--verbose', dest="verbose", default="0", help="log process to file. Options are 0 or 1  (default = 0 for no logging)")
    parser.add_argument('-V', '--version', action='version', version='annotGBK-' + version(), help="Prints version number")
    args = parser.parse_args()
    main(args)


if __name__ == '__main__':
    run()
