#!/usr/bin/python
# Write by Richard Bonnet
# Date: 20/12/2015

import os
import argparse
import filepaths
import logging

path = filepaths.determine_path()

def getSequence(afile):
	sequences_dict = {}
	if os.stat(afile).st_size != 0:
		from Bio import SeqIO
		for record in SeqIO.parse(afile, 'fasta'):
			sequences_dict[record.id] = str(record.seq)

	return sequences_dict


def makeBlastDB(blastDirectory, fasta_input_name, out_type, blast_output_name):
    #if os.path.isfile(os.path.join(path,fasta_input_name)) == True and \
    #                os.path.exists(os.path.join(path, blast_output_name+".phr")) == True and \
    #                os.path.exists(os.path.join(path, blast_output_name+".pin")) == True and \
    #                os.path.exists(os.path.join(path, blast_output_name+".psq")) == True :
    #    print "DB exists"
    #else:
    cmd = '%s -in %s -dbtype %s -out %s' % (os.path.join(blastDirectory, 'makeblastdb'), os.path.join(path, fasta_input_name), out_type, os.path.join(path, blast_output_name))
    os.system(cmd)

    return os.path.join(path, blast_output_name)


def runBlast(blastDirectory, outdir, blastType, querySeq, db_name, threads, evalue):

    outfile = os.path.join(outdir, os.path.basename(querySeq)  + '_' + os.path.splitext(os.path.basename(db_name))[0] + ".blastRes.xml")
    cmd = '%s -out %s -outfmt 5 -query %s -db %s -num_threads %s -evalue %s' % \
              (os.path.join(blastDirectory, blastType), outfile,
               querySeq, os.path.join(path, db_name), threads, evalue)
    if blastType == 'blastp':
        submitted_proteins_dict = getSequence(querySeq)
        logging.info("runBlast => start blast for inType: " + blastType)
        os.system(cmd)

    elif blastType == 'blastn':
        submitted_proteins_dict = getSequence(querySeq)
        logging.info("runBlast => start blast for inType: " + blastType)
        os.system(cmd)

    #elif blastType == 'blastx':
    #    logging.info("runBlast => fqToFsa => start")
    #    fqToFsa.main(inputSeq)
    #    logging.info("runBlast => fqToFsa => done")
    #    clean_files.append(working_directory + "/" + file_name + ".read.fsa")
    #    logging.info("runBlast => start blastX for inType: " + inType)
    #    from Bio.Blast.Applications import NcbiblastxCommandline
    #    blastCLine = NcbiblastxCommandline(query=os.path.join(working_directory, file_name + ".read.fsa"), db=os.path.join(path, "protein.db"),
    #                                       outfmt=5, out=working_directory + "/" + file_name + ".blastRes.xml",
    #                                       num_threads=threads)
    #    stdt, stdr = blastCLine()
    return outfile

def main(args):
    querySeq = args.querySeq

    outdir = args.outDir
    if outdir == '':
        outdir = os.path.dirname(querySeq)

    if args.verbose == '1':
        print "[info] logs saved to : " + str(outdir) + "/app.log"
        logging.basicConfig(filename='app.log', level=logging.DEBUG,
                            format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                            datefmt='%Y-%m-%d %I:%M:%S %p')
        logging.info('main => start runBlast')

    if args.querySeq == None:
        print "[error] missing input(s)"
        print "[info] Try: runBlast -h"
        logging.error("main => missing input(s)")
        logging.error("main => Try: runBlast -h")
        exit()

    blastType = args.blastType
    if blastType == 'blastn':
        output_type = 'nucl'
        dbBlast_output_name =  os.path.splitext(args.dbFasta)[0] + '.fna'
    elif blastType == 'blastp':
        output_type = 'prot'
        dbBlast_output_name = os.path.splitext(args.dbFasta)[0] + '.faa'
    blastDirectory = args.blastDirectory
    print '\nCreate %s blast database...' % output_type
    makeBlastDB(blastDirectory, args.dbFasta, output_type, dbBlast_output_name)
    print '%s blast database done!\n' % args.dbFasta


    threads = args.threads
    evalue = args.evalue
    print '\nStart %s-type blast with %s against %s...' % (blastType, dbBlast_output_name, querySeq)
    outfile = runBlast(blastDirectory, outdir, blastType, querySeq, dbBlast_output_name, threads, evalue)
    print 'blast done!\n'
    return outfile


def version():
    return "1.0"


def run():
    parser = argparse.ArgumentParser(description='runBlast - Version ' + version())
    parser.add_argument('-q', '--querySeq', dest="querySeq", default="query.faa", help='must be the queries in fasta format (default: query.faa)')
    parser.add_argument('-db','--dbaseFasta', dest="dbFasta", default='/usr/local/seqblast/cnr/armDB_1_GN.faa', help="fasta database name (default=/usr/local/seqblast/cnr/armDB_1_GN.faa)")
    parser.add_argument('-bt','--blastType', dest="blastType", default="blastp", help="Blast type. Options are blastp or blastn (default=blastp)")
    parser.add_argument('-t', '--threads', dest="threads", default="8", help="Number of threads used by blast. (default=8)")
    parser.add_argument('-e', '--evalue', dest="evalue", default="1e-30", help="Number of threads used by blast. (default=1e-30)")
    parser.add_argument('-o', '--outDir', dest="outDir", default="", help="Number of threads used by blast. (default=./)")
    parser.add_argument('-bd','--blastDirectory', dest="blastDirectory", default="/usr/local/ncbi-blast-2.6.0+/bin", help="Path to blast directory (default=/usr/local/ncbi-blast-2.6.0+/bin)")
    parser.add_argument('-v', '--verbose', dest="verbose", default="0", help="log process to file. Options are 0 or 1  (default = 0 for no logging)")
    parser.add_argument('-V', '--version', action='version', version='runBlast-' + version(), help="Prints version number")
    args = parser.parse_args()
    main(args)


if __name__ == '__main__':
    run()

