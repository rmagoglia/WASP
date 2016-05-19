#!/usr/bin/python

# python BAMsplitter_chr.py -b <BAMfile> -n <number_chomosomes>

import sys, argparse, logging, os, operator, subprocess, time

def splitter(inbam, chrom):    
    start_time = time.time()
    os.system("module load samtools")
    os.system("samtools index " +inbam)
    os.system('for inp in {1..'+chrom+'}; do echo chr"$inp" ; samtools view -b '+inbam+' chr"$inp" > chr"$inp".'+inbam+' ; done')
    print("minutes to complete: ", round(float(time.time()-start_time)/60),2)
    return()


def main(args, loglevel):
    logging.basicConfig(format="%(levelname)s: %(message)s", level=loglevel)
    splitter(args.bam, args.numchrom)
    
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = "Splits a (-n sorted) BAM file into -l lines per file. If reads are paired-end, make sure to split by an even number of reads so that pairs are kept together.  ### python BAMsplitter.py -b <BAMfile> -l <lines per file> ### ")
    # Specify your parameters here.
    parser.add_argument("-v", "--verbose", help="increase output verbosity", action="store_true")
    parser.add_argument('-b', '--BAM', required=True, dest='bam', help='BAM file to split, samtools sort -n ')
    parser.add_argument('-n', '--chrom', required=True, dest='numchrom', help='Number of chromosomes to split')
    
    args = parser.parse_args()
  
    # Setup logging
    if args.verbose:
        loglevel = logging.DEBUG
    else:
        loglevel = logging.INFO
  
    
    main(args, loglevel) 
