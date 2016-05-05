from __future__ import print_function
import argparse
import gzip
import pysam
from collections import defaultdict, Counter
from glob import glob
from os import path
import pandas as pd
import sys
import itertools as it

try:
    from functools import reduce
    from operator import mul
except ImportError as exc:
    # We better hope we're in Python 2.
    print(exc)
    pass

MAX_SEQS_PER_READ = 32

def product(iterable):
    return reduce(mul, iterable, 1)

def get_snps(snpdir):
    snp_dict = defaultdict(lambda : defaultdict(set))
    for fname in glob(path.join(snpdir, '*.txt.gz')):
        chrom = path.basename(fname).split('.')[0]
        i = -1
        for i, line in enumerate(gzip.open(fname, 'rt', encoding='ascii')):
            pos, ref, alt = line.split()
            snp_dict[chrom][int(pos)-1].update([ref, alt])
    return snp_dict

def get_indels(snpdict):
    indel_dict = defaultdict(lambda : defaultdict(bool))
    for chrom in snp_dict:
        for pos, alleles in snp_dict[chrom].items():
            indel_dict[chrom][pos] =  ('-' in alleles) or (max(len(i) for i in alleles) > 1)
    return indel_dict

def get_read_seqs(read, snp_dict, indel_dict, dispositions):
    num_snps = 0
    seqs = [read.seq]

    chrom = read.reference_name
    for (read_pos, ref_pos, read_base) in (read.get_aligned_pairs(matches_only=True, with_seq=True)):
        if ref_pos == None:
            continue
        if indel_dict[chrom][ref_pos]:
            dispositions['toss_indel'] += 1
            return []
        if len(seqs) > MAX_SEQS_PER_READ:
            dispositions['toss_manysnps'] += 1
            return []

        if ref_pos in snp_dict[chrom]:
            if read_base.upper() in snp_dict[chrom][ref_pos]:
                dispositions['ref_match'] += 1
                num_snps += 1
                for new_allele in snp_dict[chrom][ref_pos].difference({read_base}):
                    for seq in list(seqs):
                        # Note that we make a copy up-front to avoid modifying
                        # the list we're iterating over
                        new_seq = seq[:read_pos] + new_allele + seq[read_pos+1:]
                        seqs.append(new_seq)
            else:
                dispositions['no_match'] += 1
        else:
            # No SNP
            pass
    if len(seqs) == 1:
        dispositions['no_snps'] += 1
    else:
        dispositions['has_snps'] += 1
    return seqs

def assign_reads(insam, snp_dict, indel_dict, is_paired=True):
    fname = insam.filename
    if isinstance(fname, bytes):
        fname = fname.decode('ascii')
    basename = fname.rsplit('.', 1)[0]
    keep = pysam.Samfile('.'.join([basename, 'keep.bam']), 
            'wb', 
            template=insam)
    remap_bam = pysam.Samfile('.'.join([basename, 'to.remap.bam']),
            'wb',
            template=insam)
    if is_paired:
        fastqs = [
                gzip.open('.'.join([basename, 'remap.fq1.gz']), 'wt'),
                gzip.open('.'.join([basename, 'remap.fq2.gz']), 'wt'),
                ]
    else:
        fastqs = [ gzip.open('.'.join([basename, 'remap.fq.gz']), 'wt'),]
    unpaired_reads = [{}, {}]
    read_results = Counter()
    for i, read in enumerate(insam):
        read_seqs = get_read_seqs(read, snp_dict, indel_dict, read_results)
        if not is_paired:
            write_read_seqs([(read, read_seqs)], keep, remap_bam, fastqs)
        elif read.is_proper_pair:
            if read.qname in unpaired_reads[read.is_read1]:
                both_read_seqs = [None, None]
                both_read_seqs[read.is_read2] = read, read_seqs
                both_read_seqs[read.is_read1] = unpaired_reads[read.is_read1].pop(read.qname)
                write_read_seqs(both_read_seqs, keep, remap_bam, fastqs)
            else:
                unpaired_reads[read.is_read2][read.qname] = read, read_seqs
        else:
            read_results['not_proper_pair'] += 1
            # Most tools assume reads are paired and do not check IDs. Drop it out.
            continue
    print()
    print(len(unpaired_reads[0]), len(unpaired_reads[1]))
    print(read_results)

            
            


def write_read_seqs(both_read_seqs, keep, remap_bam, fastqs, dropped=None):
    reads, seqs = zip(*both_read_seqs)
    assert len(reads) == len(fastqs)

    num_seqs = product(len(r[1]) for r in both_read_seqs)
    if num_seqs == 0 or num_seqs > MAX_SEQS_PER_READ:
        if dropped is not None:
            for read in reads:
                dropped.write(read)
            return
        else:
            pass
    elif num_seqs == 1:
        for read, seqs in both_read_seqs:
            keep.write(read)
    else:
        for read in reads:
            remap_bam.write(read)
        left_pos = min(r.pos for r in reads)
        right_pos = max(r.pos for r in reads)

        first = True
        # Some python fanciness to deal with single or paired end reads (or
        # n-ended reads, if such technology ever happens.
        for remap_num, read_seqs in enumerate(it.product(*seqs)):
            if first:
                first = False
                continue
            for seq, read, fastq in zip(read_seqs, reads, fastqs):
                loc_line = '{}:{}:{}:{}:{}'.format(
                        remap_num,
                        read.reference_name,
                        left_pos,
                        right_pos,
                        num_seqs-1,
                        )
                fastq.write(
                        "@{loc_line}\n{seq}\n+{loc_line}\n{qual}\n"
                        .format(
                            loc_line=loc_line,
                            seq=seq,
                            qual=read.qual)
                        )



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    
    parser.add_argument("-p", "--paired_end",
                        action='store_true',
                        dest='is_paired_end', default=False,
                        help=('Indicates that reads are '
                              'paired-end (default is single).'))
    
    parser.add_argument("-s", "--sorted",
                        action='store_true', dest='is_sorted',  default=False,
                        help=('Indicates that the input bam file'
                              ' is coordinate sorted (default is False)'))


    parser.add_argument("infile", type=pysam.Samfile, help=("Coordinate sorted bam "
                        "file."))
    snp_dir_help = ('Directory containing the SNPs segregating within the '
                    'sample in question (which need to be checked for '
                    'mappability issues).  This directory should contain '
                    'sorted files of SNPs separated by chromosome and named: '
                    'chr<#>.snps.txt.gz. These files should contain 3 columns: '
                    'position RefAllele AltAllele')
    
    parser.add_argument("snp_dir", action='store', help=snp_dir_help)

    options = parser.parse_args()
    
    snp_dict = get_snps(options.snp_dir)
    indel_dict = get_indels(snp_dict)

    assign_reads(options.infile, snp_dict, indel_dict, options.is_paired_end)
