#!/usr/bin/python
# Write by Richard Bonnet
# Date: 20/01/2016

import os
import argparse
import logging
import fnmatch

def loadFile(directory, contigExt):
    fileList = []
    for root, dirnames, filenames in os.walk(directory):
        for filename in fnmatch.filter(filenames, '*.'+contigExt):
            fileList.append(os.path.join(root, filename))
    return fileList


def loadSample(sampleFile):
    sampleDict = {}
    if os.path.exists(sampleFile) == True:
        for line in open(sampleFile):
            line = line.strip()
            if line != '' and line.startswith('#') == False:
                sampleDict[line.split('\t')[0]] = line.split('\t')[1].lower()
                print line.split('\t')[0], line.split('\t')[1].lower()
    else:
        print '\nNo sample file\n'
    return sampleDict


def runAnnotation(directory, contigExt, sampleDict, exe, hmms):

    # SEARCH CONTIG FILE
    file_list = loadFile(directory, contigExt)

    if file_list != []:
        file_list.sort()
        for filename in file_list:
            wkdir = os.path.dirname(filename)
            sample_id = os.path.splitext(os.path.basename(filename))[0].split('_')[0].split('.')[0]
            if os.path.exists(sample_id + '.gbk') == False and os.path.exists(sample_id + '.faa') == False and os.path.exists(sample_id + '.fna') == False:
                try:
                    genus = sampleDict[sample_id].split(' ')[0]
                    species = sampleDict[sample_id].split(' ')[1]

                    print '\nAnnotation of %s in %s as %s %s species' % (sample_id, wkdir, genus, species)
                    cmd = '%s --outdir %s --prefix %s --locustag %s %s --evalue 1e-30 --metagenome \
                        --usegenus --genus %s --species %s --strain %s %s --force ' %\
                        (exe, wkdir, sample_id, sample_id, hmms, genus, species, sample_id, filename)

                    print cmd
                    print 'with file %s in process...' % (filename)
                    os.system('PATH=${PATH}:/usr/local/prokka-1.11/bin; ' + cmd)
                    #print >> sys.stderr
                    print 'Annotation of %s done!\n' % sample_id
                except KeyError:
                    print '\n%s and/or genus/species not found in sample file' % sample_id
                    print 'No annotation for %s\n' % sample_id
                    continue
            else:
                print '\nAnnotation already done for %s\n' % sample_id
    else:
        print '\nNo contig file found in %s\n' % directory
        exit()


def main(args):
    inpDir = args.inpDir

    if args.inpDir == None:
        print "[error] missing input(s)"
        print "[info] Try: mergeResults -h"
        exit()

    exe = args.exe
    hmms = '--hmms ' + args.hmms
    sampleDict = loadSample(args.sampleFile)
    contigExt = args.contigExt

    print '\nStart the edition of %s file...' % args.sampleFile
    runAnnotation(inpDir, contigExt, sampleDict, exe, hmms)
    print '\nEdition of %s file done!\n' % args.sampleFile


def version():
    return "1.0"


def run():
    inpDir = './'
    parser = argparse.ArgumentParser(description='runProkka - Version ' + version())
    parser.add_argument('-d', '--inpDir', dest="inpDir", default=inpDir, help='The input directory (default: %s)' % inpDir)
    parser.add_argument('-ext','--contigExt', dest="contigExt", default="fasta", help="The extension of contig files (default: fasta)")
    parser.add_argument('-s','--sampleFile', dest="sampleFile", default=os.path.join(inpDir, 'sample.csv'), help="The tab sample file with corresponding bacterial species (default: %s)" % os.path.join(inpDir, 'sample.csv'))
    parser.add_argument('-exe','--exe', dest="exe", default="/usr/local/prokka-1.11/bin/prokka", help="Prokka location (default: /usr/local/prokka-1.11/bin/prokka)")
    parser.add_argument('-hmms','--hmms', dest="hmms", default="/usr/local/prokka-1.11/db/hmm/Resfams.hmm", help="hmms profil location (default: /usr/local/prokka-1.11/db/hmm/Resfams.hmm)")
    parser.add_argument('-V', '--version', action='version', version='runProkka-' + version(), help="Prints version number")
    args = parser.parse_args()
    main(args)


if __name__ == '__main__':
    run()
