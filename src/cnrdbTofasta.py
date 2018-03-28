#!/usr/bin/python
# Write by Richard Bonnet
# Date: 20/12/2015

import os
import argparse
import re
import logging
from dbtools import read_cnr_DBase

path = os.path.dirname(os.path.realpath(__file__))
working_directory = os.getcwd()


def checkKeyExisted(key, my_dict):
	try:
		nonNone = my_dict[key] is not None
	except KeyError:
		nonNone = False
	return nonNone


def writeFASTAfromCSV(cnr_db_name, fasta_output_name, output_type, catDB, filter_type, filter):
    noSeqList = []
    if filter_type == '' and filter == '':
        filter_type == 'entry_name'
        filter == '.*'
    cnr_data = read_cnr_DBase(os.path.join(path,cnr_db_name))
    with open(os.path.join(path, fasta_output_name), 'w') as wp:
        for item in cnr_data:
            # filter the data
            if catDB in cnr_data[item]['db_names'].split('|'):
                if re.match(filter, cnr_data[item][filter_type]):
                    if output_type == 'dna':
                        pass_eval = 1e-30
                        # search_type: blastn only (pass_evalue)
                        if cnr_data[item]["dna_snp"] == "" :
                            if checkKeyExisted("blastn_evalue", cnr_data[item]):
                                pass_eval = cnr_data[item]["blastn_evalue"]
                            if checkKeyExisted("dna_sequence", cnr_data[item]):
                                #print '>' + str(item) + "_" + cnr_data[item]["entry_name"] + \
                                #      " | search_type: blast" + " | pass_evalue: " + str(pass_eval)
                                print>> wp, ('>' + str(item) + "_" + cnr_data[item]["entry_name"] +
                                      " | search_type: blast" + " | pass_evalue: " + str(pass_eval))
                                #print cnr_data[item]["dna_sequence"]
                                print>> wp, (cnr_data[item]["dna_sequence"])
                            else:
                                noSeqList.append(item)

                        # model_type: blastn + SNP (pass_evalue + snp)
                        else:
                            snpList = cnr_data[item]['dna_snp']
                            if checkKeyExisted("blastn_evalue", cnr_data[item]):
                                pass_eval = cnr_data[item]["blastn_evalue"]
                            if checkKeyExisted("dna_sequence", cnr_data[item]):
                                #print snpList
                                print>> wp, ('>' + str(item) +  "_" + cnr_data[item]["entry_name"] + " | search_type: snp" + " | pass_evalue: " + str(pass_eval) + " | SNP: " + snpList)
                                print>> wp, (cnr_data[item]["dna_sequence"])
                            else:
                                noSeqList.append(item)

                    elif output_type == 'prot':
                        pass_eval = 1e-30
                        # search_type: blastp only (pass_evalue)
                        if cnr_data[item]["prot_snp"] == "" :
                            if checkKeyExisted("blastp_evalue", cnr_data[item]):
                                pass_eval = cnr_data[item]["blastp_evalue"]
                            if checkKeyExisted("prot_sequence", cnr_data[item]):
                                print>>wp, ('>' + str(item) +  "_" + cnr_data[item]["entry_name"] +
                                            " | search_type: blast" + " | pass_evalue: " + str(pass_eval))
                                print>>wp, (cnr_data[item]["prot_sequence"])
                            else:
                                 noSeqList.append(item)

                        # model_type: blastP + SNP (pass_evalue + snp)
                        else:
                            snpList = cnr_data[item]['prot_snp']
                            if checkKeyExisted("blastp_evalue", cnr_data[item]):
                                pass_eval = cnr_data[item]["blastp_evalue"]
                            if checkKeyExisted("prot_sequence", cnr_data[item]):
                                print>> wp, ('>' + str(item) +  "_" + cnr_data[item]["entry_name"] +
                                             " | search_type: snp" + " | pass_evalue: " +
                                             str(pass_eval) + " | SNP: " + snpList)
                                print>> wp, (cnr_data[item]["prot_sequence"])
                            else:
                                noSeqList.append(item)

	return os.path.join(path, fasta_output_name)
	# get a list of models who are incomplete noSeqList


def main(args):
	if args.verbose == '1':
		print "[info] logs saved to : " + str(working_directory) + "/app.log"
		logging.basicConfig(filename='app.log', level=logging.DEBUG, format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s', datefmt='%Y-%m-%d %I:%M:%S %p')
		logging.info('main => start cnrdbTofasta')

	if args.cnrDB == None:
		print "[error] missing input(s)"
		print "[info] Try: cnrdbTofasta -h"
		logging.error("main => missing input(s)")
		logging.error("main => Try: cnrdbTofasta -h")
		exit()

	cnr_db_name = args.cnrDB
	output_type = args.type
	catDB = args.catDB
	if output_type == 'dna':
		fasta_output_name = args.output + '_' + catDB + '.fna'
	else:
		fasta_output_name = args.output + '_' + catDB + '.faa'

	filter_type, filter = args.filter_type, args.filter
	print '\nExtract %s fasta sequence from %s database as %s file...' % (output_type, cnr_db_name, fasta_output_name)
	fasta_file_name = writeFASTAfromCSV(cnr_db_name, fasta_output_name, output_type, catDB, filter_type, filter)
	print '%s file done!\n' % fasta_output_name
	return fasta_file_name


def version():
	return "1.0"


def run():
    cnr_database = "cnr/armDB_1.csv"
    out_database = os.path.splitext(cnr_database)[0]
    parser = argparse.ArgumentParser(description='cnrdbTofasta - Version ' + version())
    parser.add_argument('-cnr', '--cnrDB', dest="cnrDB", default=cnr_database, help='must be the CNR database in csv format (default: %s)' % cnr_database)
    parser.add_argument('-o', '--outFile', dest="output", default=out_database, help="Prefix of output fasta file (default= %s)" % out_database)
    parser.add_argument('-t', '--type', dest="type", default='dna', help="This option is used to write DNA or protein sequences. Options are prot or dna (default=prot)")
    parser.add_argument('-cat', '--catDB', dest="catDB", default='all', help="tag for a category in database (default= complete)")
    parser.add_argument('-ft', '--filterType', dest="filter_type", default="entry_name", help="Indicate the field of database used for apply a filter. (default='entry_name' for none)")
    parser.add_argument('-f', '--filter', dest="filter", default=".*", help="Specify the python regular expression to apply as filter (default='.*' for none)")
    parser.add_argument('-v', '--verbose', dest="verbose", default="0", help="log process to file. Options are 0 or 1  (default = 0 for no logging)")
    parser.add_argument('-V', '--version', action='version', version='cnrdbTofasta-' + version(), help="Prints version number")
    args = parser.parse_args()
    main(args)


if __name__ == '__main__':
	run()
