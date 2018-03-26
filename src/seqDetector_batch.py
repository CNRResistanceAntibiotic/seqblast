#!/usr/bin/python
# Write by Richard Bonnet
# Date: 19/01/2016
import fnmatch
import os
import argparse


def recursive_fileList(directory, filter):
    fileList = []
    for root, dirnames, filenames in os.walk(directory):
        for filename in fnmatch.filter(filenames, filter):
            fileList.append(os.path.join(root, filename))
    return fileList


def setupblastDB(csvdbName, catDB, dbType, force):
    outdbFasta = os.path.splitext(csvdbName)[0]
    if dbType == 'dna':
        outdbFasta = outdbFasta + '.fna'
    else:
        outdbFasta = outdbFasta + '.faa'
    if os.path.exists(outdbFasta) == False or force == True:
        cmd = './cnrdbTofasta.py -cnr %s -cat %s -o %s -t %s ' % (csvdbName, catDB, os.path.splitext(csvdbName)[0], dbType)
        #print cmd
        os.system(cmd)
    else:
        print '\nFasta DB file %s already exists!\n' % outdbFasta


def runBlast(fasFile, fasdbName, outdir, blastType, threads, blastevalue, force):
    blastDone = False
    xmlFile = os.path.join(outdir, os.path.basename(fasFile) + '_' + os.path.splitext(os.path.basename(fasdbName))[0] + '.blastRes.xml')
    csvFile = xmlFile + '.csv'
    if os.path.exists(csvFile) == False or force == True:
        cmd = './runBlast.py -q %s -db %s -o %s -bt %s -t %s -e %s' % (fasFile, fasdbName, outdir, blastType, threads, blastevalue)
        #print cmd
        os.system(cmd)
        blastDone = True
    else:
        print '\nBlast parsing file %s detected\nBlast already done!\n' % csvFile
    return xmlFile, csvFile, blastDone


def parseBlast(fasFile, csvdbName, blastresFile, outdir, sequenceType):
    cmd = './parseBlast.py -q %s -db %s -xml %s -o %s -t %s ' % (fasFile, csvdbName, blastresFile, outdir, sequenceType)
    #print cmd
    os.system(cmd)
    cmd = 'rm %s' % blastresFile
    os.system(cmd)


def mergeResults(workdir, idperc, cov, outdir):
    cmd = './mergeResults.py -d %s -id %.2f -cv %.2f -o %s' % \
          (workdir, float(idperc), float(cov), outdir)
    #print cmd
    os.system(cmd)


def annotGBK(gbkFile, csvFile, idperc, cov, outdir):
    cmd = './annotGBK.py -gbk %s -csv %s -id %.2f -cv %.2f -o %s' % \
          (gbkFile, csvFile, float(idperc), float(cov), outdir)
    print '\nEdition of gbk file %s' % gbkFile
    os.system(cmd)
    print 'Edition of gbk file %s done!\n' % gbkFile


def main(args):

    force = args.Force
    wkdir = args.wkDir

    csvdbName = args.cnrDB
    catDB = args.catDB
    sequenceType = args.seqType
    if sequenceType == '1' or sequenceType == '2':
        dbType = 'dna'
        fasdbPrefix = os.path.splitext(csvdbName)[0] + '_' + catDB
        fasdbName = fasdbPrefix + '.fna'
        blastType = 'blastn'
        if sequenceType == '1':
            fasFiles = recursive_fileList(wkdir, '*.fna')
        if sequenceType == '2':
            fasFiles = recursive_fileList(wkdir, '*.contigs.fasta')
    else:
        dbType = 'prot'
        fasdbPrefix = os.path.splitext(csvdbName)[0]  + '_' + catDB
        fasdbName = fasdbPrefix + '.faa'
        blastType = 'blastp'
        fasFiles = recursive_fileList(wkdir, '*.faa')

    setupblastDB(csvdbName, catDB, dbType, force)

    threads = args.threads
    blastevalue = args.evalue
    idperc = args.idPerc
    cov = args.Cov

    gbkFiles = recursive_fileList(wkdir, '*.gbk')
    fileList = []
    fasFiles.sort()
    for fasFile in fasFiles:
        gbkfound = False
        for gbkFile in gbkFiles:
            if os.path.splitext(fasFile)[0] == os.path.splitext(gbkFile)[0]:
                fileList.append((fasFile,gbkFile))
                gbkfound = True
                break
        if gbkfound == False:
            fileList.append((fasFile,''))

    for fasFile, gbkFile in fileList:
        outdir = os.path.dirname(fasFile)
        xmlblastresFile, csvblastresFile, blastDone = runBlast(fasFile, fasdbName, outdir, blastType, threads, blastevalue, force)
        if blastDone == True:
            parseBlast(fasFile, csvdbName, xmlblastresFile, outdir, sequenceType)
        if gbkFile != '':
            annotGBK(gbkFile, csvblastresFile, idperc, cov, outdir)
    mergeResults(wkdir, idperc, cov, wkdir)


def version():
    return "1.0"


def run():
    parser = argparse.ArgumentParser(description='seqDetector - Version ' + version())
    parser.add_argument('-d', '--wkDir', dest="wkDir", default='/home/bacteriologie/assemblage/workdir_run10',help='The work directory containing the data')
    parser.add_argument('-ext','--contigExt', dest="contigExt", default="fasta", help="The extension of contig files (default: fasta)")
    parser.add_argument('-s','--sampleFile', dest="sampleFile", default='', help="The tab sample file with corresponding bacterial species")
    parser.add_argument('-db','--cnrDB', dest="cnrDB", default='cnr/armDB_dna_1.csv', help="A CNR database as csv file (default: cnr/armDB_1.csv)")
    parser.add_argument('-cat','--catDB', dest="catDB", default='GN', help="Category of record in CNR database (default: all)")
    parser.add_argument('-st','--seqType', dest="seqType", default='2', help="0 for protein (blastp), 1 for CDS (blastn), 2 for contigs (blastn) (default: 0)")
    parser.add_argument('-e','--evalue', dest="evalue", default='1e-30', help="blast evalue cutoff (default: 1e-30)")
    parser.add_argument('-t','--threads', dest="threads", default='8', help="Thread number for blast (default: 8)")
    parser.add_argument('-id','--idPerc', dest="idPerc", default='80', help="Identity percentage cutoff for annotation (default: 80)")
    parser.add_argument('-cv','--Cov', dest="Cov", default='80', help="Coverage percentage cutoff for annotation (default: 80)")
    parser.add_argument('-F','--Force', dest="Force", action='store_true', default=True, help="Force file overwrite (default: False)")
    parser.add_argument('-V', '--version', action='version', version='seqDetector-' + version(), help="Prints version number")
    args = parser.parse_args()
    main(args)


if __name__ == '__main__':
    run()
