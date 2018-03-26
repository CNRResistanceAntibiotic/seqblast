#!/usr/bin/python
# Write by Richard Bonnet
# Date: 20/01/2016

import os
import argparse
import logging
import fnmatch

def runAnnotation(exe, fastaFile, sampleID, evalue, genus, species, hmms):

    wkdir = os.path.dirname(fastaFile)

    if hmms != '':
        hmms = '--hmms %s' % hmms
    if genus != '':
        genus = '--usegenus --genus %s' % genus
    if species != '':
        species = '--species %s' % species
    if evalue != '':
        evalue = '--evalue %s' % evalue

    print '\nAnnotation of %s in %s as %s %s species' % (sampleID, wkdir, genus, species)
    cmd = '%s --outdir %s --prefix %s --locustag %s --metagenome %s %s %s %s --strain %s %s --force ' %\
        (exe, wkdir, sampleID, sampleID, hmms, genus, species, evalue, sampleID, fastaFile)
    print cmd
    print 'with file %s in process...' % (fastaFile)
    os.system('PATH=${PATH}:/usr/local/prokka-1.11/bin; ' + cmd)
    #print >> sys.stderr
    cmd = 'cp %s %s' % (os.path.join(wkdir, sampleID+'.gbf'), os.path.join(wkdir, sampleID+'.gbk'))
    os.system(cmd)
    print 'Annotation of %s done!\n' % fastaFile



def main(args):

    fastaFile = args.fastaFile
    sampleID = args.sampleID
    genus  = args.genus
    species = args.genus
    evalue = args.evalue
    exe = args.exe
    hmms = args.hmms

    print '\nStart the annotation from %s file...' % fastaFile
    runAnnotation(exe, fastaFile, sampleID, evalue, genus, species, hmms)
    print '\nAnnotation from %s file is done!\n' % fastaFile


def version():
    return "1.0"


def run():
    inpDir = '/home/bacteriologie/assemblage/workdir_Run9_res'
    parser = argparse.ArgumentParser(description='runProkka - Version ' + version())
    parser.add_argument('-f', '--fastaFile', dest="fastaFile", help='The input fasta file')
    parser.add_argument('-id', '--sampleID', dest="sampleID", help='The input fasta file')
    parser.add_argument('-g','--genus', dest="genus", default="", help="Genus of the strain")
    parser.add_argument('-s','--species', dest="species", default="", help="Species of the strain")
    parser.add_argument('-e','--evalue', dest="evalue", default="1e-30", help="E value (default: 1e-30)")
    parser.add_argument('-exe','--exe', dest="exe", default="/usr/local/prokka-1.11/bin/prokka", help="Prokka location (default: /usr/local/prokka-1.11/bin/prokka)")
    parser.add_argument('-hmms','--hmms', dest="hmms", default="/usr/local/prokka-1.11/db/hmm/Resfams.hmm", help="hmms profil location (default: /usr/local/prokka-1.11/db/hmm/Resfams.hmm)")
    parser.add_argument('-V', '--version', action='version', version='runProkka-' + version(), help="Prints version number")
    args = parser.parse_args()
    main(args)


if __name__ == '__main__':
    run()
