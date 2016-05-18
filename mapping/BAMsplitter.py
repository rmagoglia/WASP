#!/usr/bin/python

# python BAMsplitter.py -b <BAMfile> -l <lines per file(even for PE reads)>

# splitter("asw1.map1ns.bam", 10000000)
import sys, argparse, logging, os, operator, subprocess, time

def splitter(inbam, numberlines):    
    start_time = time.time()
    ### remove .bam extension
    prefix = os.path.splitext(inbam)[0] 
    #only output mapped reads
    os.system("samtools view -F4 " + inbam + " > " + prefix + ".sam")
    os.system("split -d --additional-suffix=.splitsam -l "+ str(numberlines) +" " + prefix + ".sam" + " " + prefix + "_")
    ### save as BAM, reheader, rename
    os.system("samtools view -H " + inbam + " > header.sam")
    os.system('for inp in *.splitsam; do echo "$inp" ; cat header.sam "$inp" > "$inp"H ; samtools view -bS "$inp"H > "$inp"bam ; done')
    os.system("rename .splitsambam .bam " + prefix + "*.splitsambam")
    os.system("rm " + prefix + "*.splitsam")
    os.system("rm " + prefix + "*.splitsamH")
    os.system("rm " + prefix + "*.sam")
    print("minutes to complete: ", round(float(time.time()-start_time)/60),2)
    return()


def main(args, loglevel):
    logging.basicConfig(format="%(levelname)s: %(message)s", level=loglevel)
    splitter(args.bam, args.numlines)
    
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = "Splits a (-n sorted) BAM file into -l lines per file. If reads are paired-end, make sure to split by an even number of reads so that pairs are kept together.  ### python BAMsplitter.py -b <BAMfile> -l <lines per file> ### ")
    # Specify your parameters here.
    parser.add_argument("-v", "--verbose", help="increase output verbosity", action="store_true")
    parser.add_argument('-b', '--BAM', required=True, dest='bam', help='BAM file to split, samtools sort -n ')
    parser.add_argument('-l', '--lines', required=True, dest='numlines', help='Number of lines per file')
    
    args = parser.parse_args()
  
    # Setup logging
    if args.verbose:
        loglevel = logging.DEBUG
    else:
        loglevel = logging.INFO
  
    
    main(args, loglevel) 
 
