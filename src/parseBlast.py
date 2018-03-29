#!/usr/bin/python
# Write by Richard Bonnet
# Date: 20/12/2015

import os
#import json
import argparse
import logging
import re
from dbtools import pairwisealn, read_cnr_DBase

path = os.path.dirname(os.path.realpath(__file__))


def isfloat(value):
  try:
	float(value)
	return True
  except ValueError:
	return False


def findNumDash(subject, index):
	numDash = 0
	toFind = "-"
	stringList = list(subject)
	output = []

	for i in range(0, len(stringList)):
		if (stringList[i] == toFind):
			numDash += 1
		else:
			output.append(stringList[i])
		if (len(output) == index):
			break
	return numDash


def getSequence(afile):
	submitted_dict = {}
	if os.stat(afile).st_size != 0:
		from Bio import SeqIO
		for record in SeqIO.parse(afile, 'fasta'):
			submitted_dict[record.id] = str(record.seq)

	return submitted_dict


def parse_query_info(query_info):
    infoList = query_info.strip().split(' ')
    infoDict = {}
    infoDict['locustag'] = infoList[0]
    for item in infoList:
        if ':' in item:
            try:
                key, value = item.split(':')
                if key in ['scaffold','contig','start','end','strand','description']:
                    infoDict[key.strip()] = value.strip().encode('ascii', 'replace')
            except ValueError:
                continue
	return infoDict


def found_snp(aln_query, aln_subject, snpDict, start=1):
	start = start - 1
	subdicList = []
	snpdicList = []
	subj, indexj, subi, indexi = '', '', '', ''
	for i in [i for i in xrange(len(aln_subject)) if aln_subject[i] != aln_query[i]]:
		indexSubject = i - findNumDash(aln_subject, i)
		# indexQuery   = i - findNumDash(aln_query, i)
		if aln_subject[i] == '-':
			subi = 'in'
			indexi = i
			ori = 'i'
		else:
			subi = ''
			ori = aln_subject[i]

		if aln_query[i] == '-':
			chan = 'd'
		else:
			chan = aln_query[i]

		if subi == 'in' and subi == subj and indexi == indexj + 1:
			if indexSubject + 1 not in snpDict.keys():
				try:
					subdicList[-1]["change"] = subdicList[-1]["change"] + chan
				except IndexError:
					subdicList.append({"original": ori, "change": chan, "position": indexSubject + 1 + start})
			else:
				if chan in snpDict[indexSubject + 1]["change"]:
					try:
						snpdicList[-1]['change'] = snpdicList[-1]['change'] + chan
					except IndexError:
						snpdicList.append({"original": ori, "change": chan, "position": indexSubject + 1 + start})
				else:
					try:
						snpdicList[-1]['change'] = snpdicList[-1]['change'] + '[' + chan + ']'
					except IndexError:
						snpdicList.append({"original": ori, "change": '[' + chan + ']', "position": indexSubject + 1 + start})
		else:
			if indexSubject + 1 not in snpDict.keys():
				subdicList.append({"original": ori, "change": chan, "position": indexSubject + 1 + start})
			else:
				if chan in snpDict[indexSubject + 1]["change"]:
					snpdicList.append({"original": ori, "change": chan, "position": indexSubject + 1 + start})
				else:
					snpdicList.append({"original": ori, "change": '[' + chan + ']', "position": indexSubject + 1 + start})
		subj = subi
		indexj = indexi
	return snpdicList, subdicList


def readBlast(querySeq, sequenceType, xmlblastFile, cnrDB, criteria):

	submitted_dict = getSequence(querySeq)

	result_handle  = open(xmlblastFile, 'r')

	logging.info("parseBlast => start readBlast")
	from Bio.Blast import NCBIXML
	blast_records = NCBIXML.parse(result_handle)
	blastResults = {}
	pjson = {}

	json_data = read_cnr_DBase(cnrDB)

	for blast_record in blast_records:
		#print 'Blast record:', blast_record
		perfect = {}
		strict = {}
		loose = {}
		if sequenceType != '2' and blast_record.alignments:
			blast_record.alignments.sort(key= lambda align: -max(hsp.score for hsp in align.hsps))
			score = -99999.99
			best_score = False
			
			n = 0
			for alignment in blast_record.alignments:
				for hsp in alignment.hsps:
					if float(hsp.score) >= score: 
						score = float(hsp.score)
					else:
						blast_record.alignments = blast_record.alignments[:n]
						blast_record.alignments.sort(key= lambda align: -max(hsp.identities for hsp in align.hsps))						
				n +=1

		n = 0
		for alignment in blast_record.alignments:
			if n == 1 and sequenceType != '2':
				break
			n += 1
			query_info = blast_record.query.encode('ascii', 'replace')
			print '\nQUERY: '
			print ' info: %s' % query_info
			query_infoDict = parse_query_info(query_info)
			for item in ['scaffold','contig','locustag','start','end','strand','description']:
				if item in query_infoDict.keys():
					print ' '+item+':', query_infoDict[item]

			alignTitle = alignment.title
			spacepos = alignTitle.index(' ')
			hitid = alignTitle[0:spacepos].encode('ascii', 'replace')
			hit_infoList = alignTitle[spacepos:].strip().split(' | ')
			keyDB = hit_infoList[0].split('_')[0]
			entry_name = hit_infoList[0].split('_')[1].encode('ascii', 'replace')
			search_type = hit_infoList[1].split('search_type: ')[1]
			if search_type == '':
				search_type = 'blast'
			passevalue = float(hit_infoList[2].split('pass_evalue: ')[1].strip())
			if passevalue == '':
				passevalue = 1e-30
			if isfloat(passevalue):
				truePassEvalue = passevalue
			else:
				truePassEvalue = float(passevalue[0:passevalue.find(' ')])
			print 'HIT:'
			print ' title: %s' % alignTitle
			print ' hitid: %s' % hitid
			print ' keyDB: %s' % keyDB
			print ' Entry_name: %s' % entry_name
			print ' search type: %s' % search_type
			print ' passe value:', truePassEvalue

			snpDict = {}
			if search_type == 'snp':
				snpL = hit_infoList[3].strip().split(',')
				snpL = alignTitle.split(' | SNP: ')[-1].split(',')
				for eachsnp in snpL:
					r = re.compile("([a-zA-Z]+)([0-9]+)([a-zA-Z]+)")
					m = r.match(eachsnp)
					if m:
						ori = m.group(1)
						chan = m.group(3)
						if ori == 'i':
							ori = '-'
						if chan == 'd':
							chan = '-'
						pos = int(m.group(2)) #eachsnp[1:-1])
						if pos not in snpDict.keys():
							snpDict[pos] = {"orginal":ori, "change":[chan]}
						else:
							snpDict[pos]["change"] = snpDict[pos]["change"] + [chan]
						print " original:", eachsnp[0], "\tchange:", eachsnp[-1], "\tposition:", eachsnp[1:-1]
			logging.info("parseBlast => [info] | modelTypeID = " + str(alignTitle))

			init = 0
			for hsp in alignment.hsps:
				querySeq = hsp.query.replace('-', '')
				realQueryLength = len(querySeq)

				sbjctSeq = hsp.sbjct.replace('-', '')
				realSbjctLength = len(sbjctSeq)

				if   sequenceType == '1':
					db_sequence = str(json_data[keyDB]["dna_sequence"]).upper()
				elif sequenceType == '0':
					db_sequence = str(json_data[keyDB]["prot_sequence"]).upper()
				elif sequenceType == '2':
					db_sequence = str(json_data[keyDB]["dna_sequence"]).upper()
				query_sequence = str(submitted_dict[query_infoDict['locustag']]).upper()

				if db_sequence == query_sequence:
					pinsidedict = {}
					pinsidedict['sequence_type'] = sequenceType
					pinsidedict["type_match"] = "Perfect"
					pinsidedict["keyDB"] = keyDB
					pinsidedict["entry_name"] = json_data[keyDB]["entry_name"]
					pinsidedict["search_type"] = search_type
					pinsidedict["pass_evalue"] = passevalue
					pinsidedict["protDB_sequence"] = json_data[keyDB]["prot_sequence"]
					pinsidedict["dnaDB_sequence"] = json_data[keyDB]["dna_sequence"]
					if sequenceType == '0':
						pinsidedict["SNP_search"] = json_data[keyDB]["prot_snp"]
					else:
						pinsidedict["SNP_search"] = json_data[keyDB]["dna_snp"]
					pinsidedict["mol_type"] = json_data[keyDB]["mol_type"]
					pinsidedict["dna_snp"] = json_data[keyDB]["dna_snp"]
					pinsidedict["prot_snp"] = json_data[keyDB]["dna_snp"]
					pinsidedict["blastp_evalue"] = json_data[keyDB]["blastp_evalue"]
					pinsidedict["blastn_evalue"] = json_data[keyDB]["blastn_evalue"]
					pinsidedict["cluster90_grp_name"] = json_data[keyDB]["cluster90_grp_name"]
					pinsidedict["operon_grp_name"] = json_data[keyDB]["operon_grp_name"]
					pinsidedict["function_grp_names"] = json_data[keyDB]["function_grp_names"]
					pinsidedict["mechanism_names"] = json_data[keyDB]["mechanism_names"]
					pinsidedict["db_names"] = json_data[keyDB]["db_names"]

					pinsidedict["query_locustag"] = query_infoDict['locustag']

					if sequenceType == '0' or sequenceType == '1':
						for item in [("query_strand",'strand'), ("query_start",'start'), ("query_end",'end'),
									 ("query_from",'contig'),("query_from",'scaffold')]:
							key, value = item[0], item[1]
							if value in query_infoDict.keys():
								pinsidedict[key] = query_infoDict[value]
							else:
								pinsidedict[key] = ''

					elif sequenceType == '2':
						pinsidedict['query_from']   = blast_record.query.encode('ascii', 'replace')
						pinsidedict['query_start']  = hsp.query_start
						pinsidedict['query_end']    = hsp.query_start + realQueryLength - 1
						if hsp.sbjct_start < hsp.sbjct_end and hsp.query_start < hsp.query_end:
							query_strand = '1'
						elif hsp.sbjct_start > hsp.sbjct_end and hsp.query_start > hsp.query_end:
							query_strand = '1'
						else:
							query_strand = '-1'
						pinsidedict['query_strand'] = query_strand
						print ' query_from', blast_record.query.encode('ascii', 'replace')
						print ' query_start', hsp.query_start
						print ' query_end', hsp.query_start + realQueryLength
						print ' query_strand', query_strand

					pinsidedict["blast_evalue"] = hsp.expect
					print "Match: perfect",
					print ' blast evalue:', hsp.expect
					pinsidedict["blast_bit-score"] = hsp.bits
					pinsidedict["blast_max-identities"] = hsp.identities
					pinsidedict["blast_query"] = hsp.query.encode('ascii', 'replace')
					pinsidedict["blast_match"] = hsp.match.encode('ascii', 'replace')
					pinsidedict["blast_sbjct"] = hsp.sbjct.encode('ascii', 'replace')
					pinsidedict["blast_query_start"] = hsp.query_start
					pinsidedict["blast_query_end"] = hsp.query_start + realQueryLength
					pinsidedict["blast_query_From"] = blast_record.query.encode('ascii', 'replace')
					pinsidedict["query_sequence"] = query_sequence

					pinsidedict["id_per"] = 100.00
					pinsidedict["query_cov"] = 100.00
					pinsidedict["subject_cov"] = 100.00
					print " id_perc:", 100.00
					print " query_cov:", 100.00
					print " subject_cov:", 100.00

					pinsidedict["SNP_found"] = ''
					pinsidedict["substitutions_found"] = ''

					init += 1
					perfect[hitid + "|hsp_num:" + str(init)] = pinsidedict
					print 'Store ->', hitid + "|hsp_num:" + str(init) + ' as perfect match'

				elif hsp.expect <= float(truePassEvalue):
					sinsidedict = {}
					sinsidedict['sequence_type'] = sequenceType
					sinsidedict["type_match"] = "Strict"
					sinsidedict["keyDB"] = keyDB
					sinsidedict["entry_name"] = json_data[keyDB]["entry_name"]
					sinsidedict["search_type"] = search_type
					sinsidedict["pass_evalue"] = passevalue
					sinsidedict["protDB_sequence"] = json_data[keyDB]["prot_sequence"]
					sinsidedict["dnaDB_sequence"] = json_data[keyDB]["dna_sequence"]
					if sequenceType == '0':
						sinsidedict["SNP_search"] = json_data[keyDB]["prot_snp"]
					else:
						sinsidedict["SNP_search"] = json_data[keyDB]["dna_snp"]
					sinsidedict["mol_type"] = json_data[keyDB]["mol_type"]
					sinsidedict["dna_snp"] = json_data[keyDB]["dna_snp"]
					sinsidedict["prot_snp"] = json_data[keyDB]["dna_snp"]
					sinsidedict["blastp_evalue"] = json_data[keyDB]["blastp_evalue"]
					sinsidedict["blastn_evalue"] = json_data[keyDB]["blastn_evalue"]
					sinsidedict["cluster90_grp_name"] = json_data[keyDB]["cluster90_grp_name"]
					sinsidedict["operon_grp_name"] = json_data[keyDB]["operon_grp_name"]
					sinsidedict["function_grp_names"] = json_data[keyDB]["function_grp_names"]
					sinsidedict["mechanism_names"] = json_data[keyDB]["mechanism_names"]
					sinsidedict["db_names"] = json_data[keyDB]["db_names"]

					sinsidedict["query_locustag"] = query_infoDict['locustag']
					if sequenceType == '0' or sequenceType == '1':
						for item in [("query_strand",'strand'), ("query_start",'start'), ("query_end",'end'),
									 ("query_from",'contig'),("query_from",'scaffold')]:
							key, value = item[0], item[1]
							if value in query_infoDict.keys():
								sinsidedict[key] = query_infoDict[value]
							else:
								sinsidedict[key] = ''
					elif sequenceType == '2':
						sinsidedict['query_from']   = blast_record.query.encode('ascii', 'replace')
						sinsidedict['query_start']  = hsp.query_start
						sinsidedict['query_end']    = hsp.query_start + realQueryLength - 1
						if hsp.sbjct_start < hsp.sbjct_end and hsp.query_start < hsp.query_end:
							query_strand = '1'
						elif hsp.sbjct_start > hsp.sbjct_end and hsp.query_start > hsp.query_end:
							query_strand = '1'
						else:
							query_strand = '-1'
						sinsidedict['query_strand'] = query_strand
						print ' query_from', blast_record.query.encode('ascii', 'replace')
						print ' query_start', hsp.query_start
						print ' query_end', hsp.query_start + realQueryLength
						print ' query_strand', query_strand

					sinsidedict["blast_evalue"] = hsp.expect
					print "Match: Strict"
					print ' blast evalue:', hsp.expect
					sinsidedict["blast_bit-score"] = hsp.bits
					sinsidedict["blast_max-identities"] = hsp.identities
					sinsidedict["blast_query"] = hsp.query.encode('ascii', 'replace')
					sinsidedict["blast_match"] = hsp.match.encode('ascii', 'replace')
					sinsidedict["blast_sbjct"] = hsp.sbjct.encode('ascii', 'replace')
					sinsidedict["blast_query_start"] = hsp.query_start
					sinsidedict["blast_query_end"] = hsp.query_start + realQueryLength
					sinsidedict["blast_query_From"] = blast_record.query.encode('ascii', 'replace')
					sinsidedict["query_sequence"] = query_sequence

					if sequenceType == '0' or sequenceType == '1':
						aln_subject, aln_query, idper, subject_cov, query_cov = pairwisealn(db_sequence.upper(), query_sequence.upper())
						sinsidedict["id_per"] = idper
						sinsidedict["query_cov"] = query_cov
						sinsidedict["subject_cov"] = subject_cov
						snpdicList, subdicList = found_snp(aln_query, aln_subject, snpDict)
					elif sequenceType == '2':
						sinsidedict["id_per"] = round((100 * hsp.identities)/(float(len(hsp.query.encode('ascii', 'replace')))), 2)
						sinsidedict["query_cov"] = round((len(hsp.query.encode('ascii', 'replace').replace('-', '')) * 100) / float(len(query_sequence)),2)
						sinsidedict["subject_cov"] = round((len(hsp.sbjct.encode('ascii', 'replace').replace('-', '')) * 100) / float(len(db_sequence)),2)
						snpdicList, subdicList = found_snp(hsp.query.encode('ascii', 'replace'), hsp.sbjct.encode('ascii', 'replace'), snpDict, hsp.query_start)

					print " id_perc:", sinsidedict["id_per"]
					print " query_cov:", sinsidedict["query_cov"]
					print " subject_cov:", sinsidedict["subject_cov"]

					snp = ''
					for snpDict in snpdicList:
						snp = snp + snpDict["original"] + str(snpDict["position"]) + snpDict["change"] + ','
					sinsidedict["SNP_found"] = snp[:-1]
					print ' SNP:', snp[:-1]

					sub = ''
					for subDict in subdicList:
						sub = sub + subDict["original"] + str(subDict["position"]) + subDict["change"] + ','
					sinsidedict["substitutions_found"] = sub[:-1]
					print ' SUB:', sub[:-1]

					if snp == '' and sub == '' and \
					   sinsidedict["id_per"] == 100.0 and \
					   sinsidedict["query_cov"] == 100.00 and \
					   sinsidedict["subject_cov"] == 100.00:
						sinsidedict["type_match"] = "Perfect"
						print 'Correction -> Match: Perfect'

					init += 1
					if sinsidedict["type_match"] == "Perfect":
						perfect[hitid + "|hsp_num:" + str(init)] = sinsidedict
						print 'Store ->', hitid + "|hsp_num:" + str(init) + ' as perfect match'
					else:
						strict[hitid + "|hsp_num:" + str(init)] = sinsidedict
						print 'Store ->', hitid + "|hsp_num:" + str(init) + ' as strict match'

				else:
					linsidedict = {}
					linsidedict['sequence_type'] = sequenceType
					linsidedict["type_match"] = "Bad"
					linsidedict["keyDB"] = keyDB
					linsidedict["entry_name"] = json_data[keyDB]["entry_name"]
					linsidedict["search_type"] = search_type
					linsidedict["pass_evalue"] = passevalue
					linsidedict["protDB_sequence"] = json_data[keyDB]["prot_sequence"]
					linsidedict["dnaDB_sequence"] = json_data[keyDB]["dna_sequence"]
					if sequenceType == '0':
						linsidedict["SNP_search"] = json_data[keyDB]["prot_snp"]
					else:
						linsidedict["SNP_search"] = json_data[keyDB]["dna_snp"]
					linsidedict["mol_type"] = json_data[keyDB]["mol_type"]
					linsidedict["dna_snp"] = json_data[keyDB]["dna_snp"]
					linsidedict["prot_snp"] = json_data[keyDB]["dna_snp"]
					linsidedict["blastp_evalue"] = json_data[keyDB]["blastp_evalue"]
					linsidedict["blastn_evalue"] = json_data[keyDB]["blastn_evalue"]
					linsidedict["cluster90_grp_name"] = json_data[keyDB]["cluster90_grp_name"]
					linsidedict["operon_grp_name"] = json_data[keyDB]["operon_grp_name"]
					linsidedict["function_grp_names"] = json_data[keyDB]["function_grp_names"]
					linsidedict["mechanism_names"] = json_data[keyDB]["mechanism_names"]
					linsidedict["db_names"] = json_data[keyDB]["db_names"]

					linsidedict["query_locustag"] = query_infoDict['locustag']
					if sequenceType == '0' or sequenceType == '1':
						for item in [("query_strand", 'strand'), ("query_start", 'start'), ("query_end", 'end'),
									 ("query_from", 'contig'), ("query_from", 'scaffold')]:
							key, value = item[0], item[1]
							if value in query_infoDict.keys():
								linsidedict[key] = query_infoDict[value]
							else:
								linsidedict[key] = ''
					elif sequenceType == '2':
						linsidedict['query_from']   = blast_record.query.encode('ascii', 'replace')
						linsidedict['query_start']  = hsp.query_start
						linsidedict['query_end']    = hsp.query_start + realQueryLength - 1
						if hsp.sbjct_start < hsp.sbjct_end and hsp.query_start < hsp.query_end:
							query_strand = '1'
						elif hsp.sbjct_start > hsp.sbjct_end and hsp.query_start > hsp.query_end:
							query_strand = '1'
						else:
							query_strand = '-1'
						linsidedict['query_strand'] = query_strand
						print ' query_from', blast_record.query.encode('ascii', 'replace')
						print ' query_start', hsp.query_start
						print ' query_end', hsp.query_start + realQueryLength
						print ' query_strand', query_strand

					linsidedict["blast_evalue"] = hsp.expect
					print "Match: Bad"
					print ' blast evalue:', hsp.expect
					linsidedict["blast_bit-score"] = hsp.bits
					linsidedict["blast_max-identities"] = hsp.identities
					linsidedict["blast_query"] = hsp.query.encode('ascii', 'replace')
					linsidedict["blast_match"] = hsp.match.encode('ascii', 'replace')
					linsidedict["blast_sbjct"] = hsp.sbjct.encode('ascii', 'replace')
					linsidedict["blast_query_start"] = hsp.query_start
					linsidedict["blast_query_end"] = hsp.query_start + realQueryLength
					linsidedict["blast_query_From"] = blast_record.query.encode('ascii', 'replace')
					linsidedict["protQR_sequence"] = query_sequence

					if sequenceType == '0' or sequenceType == '1':
						aln_subject, aln_query, idper, subject_cov, query_cov = pairwisealn(db_sequence.upper(), query_sequence.upper())
						linsidedict["id_per"] = idper
						linsidedict["query_cov"] = query_cov
						linsidedict["subject_cov"] = subject_cov
					elif sequenceType == '2':
						linsidedict["id_per"] = round((100 * hsp.identities)/(float(len(hsp.query.encode('ascii', 'replace')))), 2)
						linsidedict["query_cov"] = round((len(hsp.query.encode('ascii', 'replace').replace('-',''))*100) / float(len(query_sequence)),2)
						linsidedict["subject_cov"] = round((len(hsp.sbjct.encode('ascii', 'replace').replace('-',''))*100) / float(len(db_sequence)),2)

					print " id_perc:", linsidedict["id_per"]
					print " query_cov:", linsidedict["query_cov"]
					print " subject_cov:", linsidedict["subject_cov"]

					linsidedict["SNP_found"] = 'NA'
					linsidedict["Substitutions_found"] = 'NA'

					init += 1
					loose[hitid + "|hsp_num:" + str(init)] = linsidedict
					print 'Store ->', hitid + "|hsp_num:" + str(init) + ' as bad match'

		if len(perfect) == 0 and len(strict) == 0:
			if criteria == "0":
				blastResults[blast_record.query.encode('ascii', 'replace')] = loose
				logging.info("parseBlast => hit = " + str(blast_record.query.encode('ascii', 'replace')) + " => Loose")

		elif len(perfect) == 0:
			blastResults[blast_record.query.encode('ascii', 'replace')] = strict
			logging.info("parseBlast => hit = " + str(blast_record.query.encode('ascii', 'replace')) + " => Strict")

		else:
			blastResults[blast_record.query.encode('ascii', 'replace')] = perfect
			logging.info("parseBlast => hit = " + str(blast_record.query.encode('ascii', 'replace')) + " => Perfect")
	return blastResults


def write_excel(blastResults, outfile):
	pjson_tagList = ['query_locustag','query_from','query_start','query_end','query_strand',
				 'entry_name','type_match','blast_evalue','id_per','query_cov', 'subject_cov',
				 'SNP_found','SNP_search','substitutions_found',
				 'keyDB','pass_evalue','function_grp_names','cluster90_grp_name','mechanism_names','operon_grp_name',
				 'blast_query_From','blast_query_start','blast_query_end',
				 'blast_bit-score', 'blast_max-identities','blast_query','blast_match','blast_sbjct',
				 'query_sequence','protDB_sequence','dnaDB_sequence', 'blastn_evalue', 'blastp_evalue', 'dna_snp','prot_snp',
                     'mol_type','db_names']

	excel_colname = ['Query name','Query location','Query start', 'Query end', 'Query strand',
					 'Hit name', 'Match type', 'blast evalue', 'Identity (%)', 'Query coverage (%)', 'Hit coverage (%)',
					 'SNP(s) found', 'SNP(s) searched', 'Other SNP(s)',
					 'Hit key in dtbase','Pass evalue','Functional cluster','90 cluster','Mechanism', 'Operon name',
					 'blast query from','blast query start', 'blast query end',
					 'blast bit-score', 'blast max-identities', 'blast query sequence', 'blast match sequence', 'blast hit sequence',
					 'Query sequence','Hit dtbase protein sequence','Hit dtbase DNA sequence', 'blastn_evalue','blastp_evalue',
                     'dna_snp','prot_snp','mol_type','db_names']

	f = open(outfile, 'w')
	f.write('\t'.join(excel_colname)+'\n')

	queryList = blastResults.keys()
	queryList.sort()
	for query in queryList:
		print '\nQuery name: ' + query
		hitDict = blastResults[query]
		hit_nameList = hitDict.keys()
		hit_nameList.sort()
		for hit_name in hit_nameList:
			print ' Result:', hit_name
			txtList = []
			for tag in pjson_tagList:
				txtList.append(str(hitDict[hit_name][tag]))
			f.write('\t'.join(txtList)+'\n')
	f.close()


#def write_json(blastResults, outdir):
#    pjson = json.dumps(blastResults)
#    f = open(os.path.join(outdir, "blastResults.json"), 'w')
#    f.write(pjson)
#    f.close()

def main(args):
    outdir = args.outDir
    if args.verbose == '1':
        print "[info] logs saved to : " + str(outdir) + "/app.log"
        logging.basicConfig(filename='app.log', level=logging.DEBUG,format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',datefmt='%Y-%m-%d %I:%M:%S %p')
        logging.info('main => start parseBlast')
    cnrDB = args.cnrDB
    querySeq = args.querySeq
    xmlblastFile = args.xmlblastFile
    csvblastFile = os.path.join(outdir, os.path.basename(xmlblastFile) + '.csv')

    sequenceType = args.sequenceType
    criteria = args.criteria
    print '\nParsing xml blast result of %s file...' % xmlblastFile
    blastResults = readBlast(querySeq, sequenceType, xmlblastFile, cnrDB, criteria)
    write_excel(blastResults, csvblastFile)
    print '\nParsed blast result are stored in %s file!\n' % csvblastFile
    logging.info("parseBlast => done")


def version():
	return "1.0"


def run():
    #querySeq = 'query.faa'
    #xmlblastFile = querySeq + ".blastRes.xml"
    #outdir = os.path.dirname(querySeq)
    parser = argparse.ArgumentParser(description='Resistance Gene Identifier - Version ' + version())
    parser.add_argument('-q', '--querySeq', dest="querySeq", help='must be the protein queries in fasta format')
    parser.add_argument('-xml', '--xmlblastFile', dest="xmlblastFile", help="Xml blast result file")
    parser.add_argument('-t', '--sequenceType', dest="sequenceType", default="0", help="0 for protein (blastp), 1 for CDS (blastn), 2 for contigs (blastn) (default=0")
    parser.add_argument('-db', '--cnrDB', dest="cnrDB", default="/usr/local/seqblast/cnr/armDB_1.csv", help="Blast database name (default=/usr/local/seqblast/cnr/armDB_1.csv)")
    parser.add_argument('-e', '--exclude_loose',  dest="criteria", default="1", help="This option is used to include or exclude the loose hits. Options are 0 or 1 (default=1 for exclude)")
    parser.add_argument('-o', '--outDir',  dest="outDir", default="", help="The ouput directory. (default= ./")
    parser.add_argument('-v', '--verbose', dest="verbose", default="0", help="log process to file. Options are 0 or 1  (default = 0 for no logging)")
    parser.add_argument('-V', '--version', action='version', version='fastaToblastdb-' + version(), help="Prints version number")
    args = parser.parse_args()
    main(args)


if __name__ == '__main__':
	run()
