#!/usr/bin/python
import os
import argparse
import subprocess


def launch(exe, quast, plasmid, cv, pefile1, pefile2, sfiles, pacbio, sanger, trcontig, uncontig, outdir):

    if plasmid == True:
        ass_dir = os.path.join(outdir, 'plasmidspades')
        i = 0
        while os.path.exists(ass_dir):
            i += 1
            ass_dir = os.path.join(outdir, 'plasmidspades_%i' % i)

        cmd = '%s --plasmid --careful' % exe
    else:
        ass_dir = os.path.join(outdir, 'spades')
        i = 0
        while os.path.exists(ass_dir):
            i += 1
            ass_dir = os.path.join(outdir, 'spades_%i' % i)
        cmd = '%s --careful' % exe

    if pefile1 != '' and pefile2 != '':
        cmd = cmd + ' -1 %s -2 %s' % (pefile1, pefile2)

    for index, f in enumerate(sfiles):
        index += 1
        if f != '':
            cmd = cmd + ' --s%i %s' % (index, f)

    if pacbio != '':
        cmd = cmd + ' --pacbio %s' % pacbio

    if sanger != '':
        cmd = cmd + ' --sanger %s' % sanger

    if trcontig != '':
        cmd = cmd + ' --trusted-contigs %s' % trcontig

    if uncontig != '':
        cmd = cmd + ' --untrusted-contigs %s' % uncontig

    if cv != 'off':
        cmd = cmd + ' --cov-cutoff %s' % cv

    cmd = cmd + ' -o %s' % ass_dir
    print '\nSpades:\n%s\n' % cmd
    print 'Spades in process...'
    out_str = subprocess.check_output(cmd, shell=True)
    print out_str

    cmd = '%s %s -o %s' % (quast, os.path.join(ass_dir, 'scaffolds.fasta'), os.path.join(ass_dir, 'quast'))
    print '\nQuast in process...\n'
    out_str = subprocess.check_output(cmd, shell=True)
    print out_str

    cmd = 'reapr perfectmap assembly.fa short_1.fq short_2.fq 300 perfect; ' \
          'reapr smaltmap assembly.fa long_1.fq long_2.fq long_mapped.bam;' \
          'reapr pipeline assembly.fa long_mapped.bam outdir perfect'

    cmd = 'rm -Rf'
    for item in ['corrected', 'K*', 'tmp', 'misc', 'mismatch*']:
        cmd = cmd + ' %s' % os.path.join(ass_dir, item)
    os.system(cmd)

def main(args):
    pefile1 = args.pefile1
    pefile2 = args.pefile2
    sfiles  = args.upfiles.split(',')
    pacbio  = args.pacbio
    sanger  = args.sanger
    trcontig = args.trcontig
    uncontig = args.uncontig
    outdir = args.outdir
    plasmid = args.plasmid
    cv = args.cv
    exe    = args.exe
    quast  = args.quast
    print '\nOutput dir: %s' % outdir
    print 'Input file names: %s %s %s %s %s %s %s' % (pefile1, pefile2, ' '.join(sfiles), pacbio, sanger, trcontig, uncontig)
    #backup_assembly(outdir, sample)
    launch(exe, quast, plasmid, cv, pefile1, pefile2, sfiles, pacbio, sanger, trcontig, uncontig, outdir)


def version():
    return "0.0.1"


def run():
    parser = argparse.ArgumentParser(description='launch A5 assembler - Version '+version())
    parser.add_argument('-1','--forwardPE', dest='pefile1', action='store', default='sk_s1_pe.fastq', help='Forward fastq or fastqz file')
    parser.add_argument('-2','--reversePE', dest='pefile2', action='store', default='sk_s2_pe.fastq', help='Reverse fastq or fastqz file')
    parser.add_argument('-s','--unpaired', dest='upfiles', action='store', default='sk_s3_up.fastq', help='Unpaired fastq or fastqz file')
    parser.add_argument('-pb','--pacbio', dest='pacbio', action='store', default='', help='PacBio file')
    parser.add_argument('-sg','--sanger', dest='sanger', action='store', default='', help='Sanger file')
    parser.add_argument('-tc','--trustedcontig', dest='trcontig', action='store', default='', help='Trusted contig file')
    parser.add_argument('-uc','--untrustedContig', dest='uncontig', action='store', default='', help='UnTrusted contig file')
    parser.add_argument('-pl', '--plasmid', dest='plasmid', action='store_true', help='Plamide assembly')
    parser.add_argument('-cv','--cv', dest='cv', action='store', default='off', help='Coverage cutoff value (a positive float number or \'auto\' or \'off\') [default: \'off\']')
    parser.add_argument('-o', '--outputDir',  dest='outdir', action='store', help="Name of output directory")
    parser.add_argument('-exe', '--executable',  dest='exe', action='store', default='/usr/local/SPAdes-3.10.0/bin/spades.py', help="Python script location")
    parser.add_argument('-q', '--quast_executable',  dest='quast', action='store', default='/usr/local/quast-4.4/quast.py', help="Python script location")
    parser.add_argument('-V', '--version', action='version', version='rgi-'+version(), help = "Prints version number")
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = run()
    main(args)
