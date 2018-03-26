#!/usr/bin/python
import argparse
import os
import subprocess


def merging(inF, inR, output_dir, m, r, f, s):
    # http://scottmyourstone.blogspot.fr/2013/10/merging-overlapping-paired-end-reads.html
    # flash -m 25 -r 150 -f 300 -s 100 <forward reads file path> <reverse reads file path>
    outPrefix = os.path.join(output_dir, os.path.splitext(os.path.basename(inF))[0] + '_merged' )
    cmd = 'flash -m %s -r %s -f %s -s %s %s %s -d %s -o %s ' % (m, r, f, s, inF, inR, output_dir, outPrefix)
    print '\nMerging %s %s with flash:\n%s\n' % (inF, inR, cmd)
    out_str = subprocess.check_output(cmd, shell=True)
    print out_str
    return outPrefix + '.fastq'


def subset(inF, inR, output_dir, subset_size):
    # https://github.com/lh3/seqtk/
    # Subsample 10000 read pairs from two large paired FASTQ files (remember to use the same random seed to keep pairing):
    # seqtk sample -s100 read1.fq 10000 > in_F.fq
    # seqtk sample -s100 read2.fq 10000 > in_R.fq
    subset_list = []
    for index, infile in enumerate([inF, inR]):
        if infile:
            outfile = os.path.join(output_dir, 'subset_'+ str(index + 1) +'.fastq.gz')
            cmd = 'seqtk sample -s100 %s 10000 > %s' % (infile, outfile)
            print '\nSubseting %s in %s:\n%s\n' % (infile, outfile, cmd)
            out_str = subprocess.check_output(cmd, shell=True)
            print out_str
            subset_list.append(outfile)
    return subset_list


def sickle(inF, inR, output_dir, read_type, t, q, l, g, x, n):
    # https://github.com/najoshi/sickle
    # sickle se -t sanger -q 30 -l 40 -x -n -g -f in.fastq -o tr.fastq.gz
    # sickle pe -t illumina -l 80 -q 20 -g -f in_F.fastq -r in_R.fastq -o tr_F.fastq.gz -p tr_R.fastq.gz -s tr_S.fastq.gz
    if read_type == 'se':
        trim_list = [os.path.join(output_dir, 'sk_s1_se.fastq.gz')]
        cmd = 'sickle se -t %s -q %s -l %s' % (t, q, l)
        option = (' -g',' -x',' -n')
        for index, item in enumerate([g, x, n]):
            if item == 'True':
                cmd = cmd + option[index]
        cmd = cmd + ' -f %s -o %s ' % (inF, trim_list[0])
        print '\nTrimming unpaired file %s with sickle:\n%s\n' % (inF, cmd)
        out_str = subprocess.check_output(cmd, shell=True)
        print out_str
    else:
        trim_list = [os.path.join(output_dir, 'sk_s1_pe.fastq.gz'), os.path.join(output_dir, 'sk_s2_pe.fastq.gz'),
                     os.path.join(output_dir, 'sk_s3_up.fastq.gz')]
        cmd = 'sickle pe -t %s -q %s -l %s' % (t, q, l)
        option = (' -g',' -x',' -n')
        for index, item in enumerate([g, x, n]):
            if item == 'True':
                cmd = cmd + option[index]
        cmd = cmd + ' -f %s -r %s -o %s -p %s -s %s' % (inF, inR, trim_list[0], trim_list[1], trim_list[2])
        print '\nTrimming paired-end files %s and %s with sickle:\n%s\n' % (inF, inR, cmd)
        out_str = subprocess.check_output(cmd, shell=True)
        print out_str
    return trim_list


def trimmomatic(inF, inR, output_dir, read_type, trimc, triml, trimt, trimw, trimq, trimm):
    # java -jar $TRIM/trimmomatic.jar PE inF inR s1_pe s1_se s2_pe s2_se
    # ILLUMINACLIP:$TRIM/adapters/TruSeq3-PE.fa:2:30:10
    # LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:40

    # java -jar $TRIM/trimmomatic-0.30.jar SE inF s1_se.fq.gz
    # ILLUMINACLIP:$TRIM/adapters/TruSeq3-SE:2:30:10
    # LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:40

    trimmomatic_dir  = '/usr/local/Trimmomatic-0.36'
    adapters_dir = os.path.join(trimmomatic_dir, 'adapters')
    if read_type == 'se':
        trim_list = [os.path.join(output_dir, 'tr_s1_se.fastq.gz')]
        cmd = 'java -jar %s SE -threads 4 -phred33 %s %s ILLUMINACLIP:%s:%s LEADING:%s TRAILING:%s SLIDINGWINDOW:%s:%s MINLEN:%s' %\
              (os.path.join(trimmomatic_dir, 'trimmomatic.jar'), inF, trim_list[0],
               os.path.join(adapters_dir, 'adapter.fasta'), trimc, triml, trimt, trimw, trimq, trimm)
        print 'Trimming unpaired file %s with trimmomatic:\n%s\n' % (inF, cmd)
        out_str = subprocess.check_output(cmd, shell=True)
        print out_str

    else:
        trim_list = [os.path.join(output_dir, 'tr_s1_pe.fastq.gz'), os.path.join(output_dir, 'tr_s1_up.fastq.gz'),
                     os.path.join(output_dir, 'tr_s2_pe.fastq.gz'), os.path.join(output_dir, 'tr_s2_up.fastq.gz'), 
                     os.path.join(output_dir, 'tr_s3_up.fastq.gz')]
        cmd = 'java -jar %s PE -threads 4 -phred33 %s %s %s %s %s %s ILLUMINACLIP:%s:%s LEADING:%s TRAILING:%s SLIDINGWINDOW:%s:%s MINLEN:%s' % \
              (os.path.join(trimmomatic_dir, 'trimmomatic.jar'), inF, inR, trim_list[0], trim_list[1], trim_list[2], trim_list[3],
               os.path.join(adapters_dir, 'adapter.fasta'), trimc, triml, trimt, trimw, trimq, trimm)
        print 'Trimming paired-end files %s and %s with trimmomatic:\n%s\n' % (inF, inR, cmd)
        out_str = subprocess.check_output(cmd, shell=True)
        print out_str
        cmd = 'cat %s %s > %s' % (trim_list[1],trim_list[3],trim_list[4])
        out_str = subprocess.check_output(cmd, shell=True)
        cmd = 'rm %s %s' % (trim_list[1],trim_list[3])
        out_str = subprocess.check_output(cmd, shell=True)
        trim_list = [os.path.join(output_dir, 'tr_s1_pe.fastq.gz'), 
                     os.path.join(output_dir, 'tr_s2_pe.fastq.gz'), 
                     os.path.join(output_dir, 'tr_s3_up.fastq.gz')]
    return trim_list


def main(args):
    #PARSE INPUT FILE NAME
    inF = args.in_f
    inR = args.in_r
    if inR == '':
        read_type = 'se'
        args.merge = 'False'
        args.subset_size = 'all'
    else:
        read_type = 'pe'

    #PARSE OUTPUT DIRECTORY
    output_dir = args.output_dir
    if output_dir == '':
        output_dir = os.path.abspath(os.path.splitext(inF)[0])
    if not os.path.exists(output_dir):
        cmd = 'mkdir %s' % output_dir
        print 'Make output directory:\n%s\n' % cmd
        out_str = subprocess.check_output(cmd, shell=True)
        print out_str

    #PARSE MERGING ARGUMENTS FOR FLASH
    merged_file = ''
    if args.flash == True:
        flashm = args.flashm
        flashr = args.flashr
        flashf = args.flashf
        flashs = args.flashs
        merged_file = merging(inF, inR, output_dir, flashm, flashr, flashf, flashs)

    #PARSE SUBSET SIZE FOR SEQTK
    if args.subset_size == 'all':
        input_list = [inF, inR]
    else:
        input_list = subset(inF, inR, output_dir, args.subset_size)

    #PARSE TRIMMOMATIC ARGUMENTS FOR READ TRIMMING
    trim_list = []
    if args.trimmomatic == True:
        trimc = args.trimc
        triml = args.triml
        trimt = args.trimt
        trimw = args.trimw
        trimq = args.trimq
        trimm = args.trimm
        trim_list = trim_list + trimmomatic(input_list[0], input_list[1], output_dir, read_type, trimc, triml, trimt, trimw, trimq, trimm)

    # PARSE SICKLE ARGUMENTS FOR READ TRIMMING
    if args.sickle == True:
        sicklet = args.sicklet
        sickleq = args.sickleq
        sicklel = args.sicklel
        sickleg = args.sickleg
        sicklen = args.sicklen
        sicklex = args.sicklex
        trim_list = trim_list + sickle(input_list[0], input_list[1], output_dir, read_type, sicklet, sickleq, sicklel, sickleg, sicklex, sicklen)
    trim_list.append(merged_file)
    return trim_list


def version():
    return "0.0.1"


def run():
    parser = argparse.ArgumentParser(description='Timming')
    #READS INPUT FILES
    parser.add_argument('-f', '--in_forward', dest='in_f', action="store", default='input1.fastq.gz',
                        help='Input forward read file (fastq or fastq.gz) [input1.fastq.gz]')
    parser.add_argument('-r', '--in_reverse', dest='in_r', action="store", default='input2.fastq.gz',
                        help='Input reverse read file (fastq or fastq.gz) [input2.fastq.gz]')

    # OUTPUT DIR
    parser.add_argument('-o', '--output_dir', dest='output_dir', action="store", default='',
                        help='Output directory')

    #MERGE WITH FLASH
    parser.add_argument('-m', '--merge', dest='flash', action="store_true",
                        help='Merge reads with flash')
    parser.add_argument('--flashm', dest='flashm', action="store", default='25',
                        help='flash option m [25]')
    parser.add_argument('--flashr', dest='flashr', action="store", default='150',
                        help='flash option r [150]')
    parser.add_argument('--flashf', dest='flashf', action="store", default='300',
                        help='flash option f [300]')
    parser.add_argument('--flashs', dest='flashs', action="store", default='100',
                        help='flash option s [100]')

    #SUBSET SIZE WITH SEQTK
    parser.add_argument('-s', '--subset', dest='subset_size', action="store", default='all',
                        help='Size of the read subset. default:[all]')

    #TRIMMING WITH TRIMMOMATIC
    parser.add_argument('-tr', '--trimmomatic', dest='trimmomatic', action="store_true",
                        help='Trimming with trimmomatic')
    parser.add_argument('--trimc', dest='trimc', action="store", default='2:30:10',
                    help='sickle option trailing [2:30:10]')
    parser.add_argument('--triml', dest='triml', action="store", default='3',
                        help='sickle option leading [3]')
    parser.add_argument('--trimt', dest='trimt', action="store", default='3',
                        help='sickle option trailing [3]')
    parser.add_argument('--trimw', dest='trimw', action="store", default='4',
                        help='sickle option windows_size [4]')
    parser.add_argument('--trimq', dest='trimq', action="store", default='15',
                        help='sickle option quality [15]')
    parser.add_argument('--trimm', dest='trimm', action="store", default= '36',
                        help='sickle option minlen [36]')

    #TRIMMING WITH SICKLE
    parser.add_argument('-sk', '--sickle', dest='sickle', action="store_true",
                        help='Trimming with sickle')
    parser.add_argument('--sicklet', dest='sicklet', action="store", default='sanger',
                        help='sickle option type [sanger]')
    parser.add_argument('--sickleq', dest='sickleq', action="store", default='20',
                        help='sickle option quality [20]')
    parser.add_argument('--sicklel', dest='sicklel', action="store", default='40',
                        help='sickle option lenght [40]')
    parser.add_argument('--sickleg', dest='sickleg', action="store_true",
                        help='sickle option g')
    parser.add_argument('--sicklen', dest='sicklen', action="store_true",
                        help='sickle option n')
    parser.add_argument('--sicklex', dest='sicklex', action="store_true",
                        help='sickle option x')
    #VERSION
    parser.add_argument('-V', '--version', action='version', version='rgi-'+version(), help = "Prints version number")
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = run()
    trim_list = main(args)
    print 'Trimmed files:'
    print trim_list







