#!/usr/bin/python
import glob
import os
import sys
import argparse


def backup_assembly(outdir, sample):
    # BACKUP CONTIG, SCAFFOLD AND STAT FILES
    for ext in ['.assembly_stats.csv','.contigs.fasta','.contigs.fastq',\
                '.final.scaffolds.fasta', '.final.scaffolds.fastq', 'pe.sort.bam', '.pe.sort.bam.bai']:
        for filename in glob.glob(os.path.join(outdir, sample + ext)):
            if filename:
                print 'Previous assembly file %s detected' % filename
                print 'Backup of file %s as file %s_previous' % (filename, filename)
                os.rename(filename, filename + '_previous')


def launch(sample, exe, file1, file2, outdir):
    print '\nAssembly of %s in %s with %s and %s' % (sample, outdir, file1, file2)
    outdir = os.path.abspath(outdir)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    #cmd = 'cd %s; %s %s %s %s'  % (outdir, exe, file1, file2, sample)
    cmd = 'cd %s; %s --end=5 %s %s %s; rm -Rf %s; rm -Rf %s.s*; rm %s.agp; rm %s.ec*; rm %s.library*; rm %s.preproc*; rm %s.raw*; rm %s.tmp*' \
    % (outdir, exe, file1, file2, sample, sample, sample, sample, sample, sample, sample, sample, sample)
    #print cmd
    print 'File1: %s' % file1
    print 'File2: %s' % file2
    if file1.split('_R1') == file2.split('_R2'):
        print 'In process...'
        print cmd
        os.system(cmd)
        print >> sys.stderr
        print 'Assembly of %s done!' % sample


        inpstatFile = os.path.join(outdir, sample + '.assembly_stats.csv')
        if os.path.exists(inpstatFile):
            outstatFile = os.path.join(os.path.dirname(os.path.dirname(outdir)), 'a5_assembly_stats.csv')
            if os.path.exists(outstatFile):
                f = open(outstatFile, 'a')
            else :
                f = open(outstatFile, 'w')
                header = 'ID\tNombre de contigs\tNombre de scaffolds\tTaille du genome\tScaffold le plus long\tN50\t' \
                     'Nombre de reads\tNombre de reads conserves\t% de reads conservees\tNombre de nucleotides\t' \
                     'Nombre de nucleotides conserves\t% de nucleotides conserves\t' \
                     'Profondeur moyenne\tProfondeur moyenne filtree\tProfondeur mediane\tProfondeur au 10eme percentile\t' \
                     'Nombre de bases >= Q40\tGC %\n'
                f.write(header)
            for n, line in enumerate(open(inpstatFile, 'r')):
                if n > 0:
                    f.write('\t'.join(line.split('\t')[:-1])+'\n')
            f.close()


def main(args):
    file1 = args.file1
    file2 = args.file2
    sample = args.sample
    outdir = args.outdir
    exe    = args.exe
    print '\nSample: %s' % sample
    print 'Input file names: %s %s' % (file1, file2)
    print 'Output dir: %s\n' % outdir
    backup_assembly(outdir, sample)
    launch(sample, exe, file1, file2, outdir)


def version():
    return "0.0.1"


def run():
    parser = argparse.ArgumentParser(description='launch A5 assembler - Version '+version())
    parser.add_argument('-file1','--readFile_input1', dest='file1', default='file1.fastq', help='Forward fastq or fastqz file')
    parser.add_argument('-file2','--readFile_input2', dest='file2', default='file2.fastq', help='Reverse fastq or fastqz file')
    parser.add_argument('-id', '--sampleName',  dest='sample', help="ID of sample")
    parser.add_argument('-out', '--outputDir',  dest='outdir', help="Name of output directory")
    parser.add_argument('-exe', '--executable',  dest='exe', default='/usr/local/a5_miseq_linux_20140604/bin/a5_pipeline.pl', help="Perl script location")
    parser.add_argument('-V', '--version', action='version', version='rgi-'+version(), help = "Prints version number")
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = run()
    main(args)
