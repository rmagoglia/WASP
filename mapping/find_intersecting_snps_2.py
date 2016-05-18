""" find_intersecting_snps

This is a rewrite of the official WASP code.  It has a more straightforward
design, that is faster and has no maximum window size, but requires loading all
SNPs into memory at once. For reference, 70 million SNPs requires about 10GB of
RAM.

"""
from __future__ import print_function
import argparse
import gzip
import itertools as it
import bz2
import sqlite3
import os
import sys
from collections import defaultdict, Counter
from pysam import AlignmentFile as Samfile

try:
    from functools import reduce
    from operator import mul
except ImportError as exc:
    # We better hope we're in Python 2.
    print(exc)

MAX_SEQS_PER_READ = 1024


class SNPDB(object):

    """A wrapper for an sqlite SNP database."""

    def __init__(self, snp_file_dir, db_file=None, overwrite=False):
        """Initialize a SNPS object from a snp_dir or bed file.

        :snp_file_dir: Either a dir with one file per chromosome in the format:
                            <pos> <ref> <alt>
                       or a simple vcf file.
        :db_file:      Optional file path for the db.
        :overwrite:    Delete db and create a new one.

        """
        self.chromosomes = {}

        if db_file:
            self.db = os.path.abspath(db_file)

        if os.path.isdir(snp_file_dir):
            if not self.db:
                self.db = os.path.abspath(os.path.join(snp_file_dir, 'snps.db'))
            if overwrite and os.path.isfile(self.db):
                os.remove(self.db)
            if os.path.isfile(self.db):
                self._initdb()
                self.length = len(self)  # Make sure chromosome dict is set
                return
            files = [os.path.abspath(os.path.join(snp_file_dir, i)) for i in
                     os.listdir(snp_file_dir)
                     if i.split('.')[0].startswith('chr')
                     or i.split('.')[0].isdigit()]
            if not files:
                raise self.SNP_Error('{} contains no recognizable SNP files'
                                     .format(snp_file_dir))
            sys.stderr.write('Creating SNP database {}\n'.format(self.db))
            # Create database here
            self._initdb()
            for fl in files:
                chrom = os.path.basename(fl).split('.')[0]
                self._add_table_if_none(chrom)
                self.chromosomes[chrom] = 0
                print("Loading SNPs from", fl)
                with open_zipped(fl) as fin:
                    tln = fin.readline().rstrip().split('\t')
                    if len(tln) != 3 or not tln[0].isdigit():
                        sys.stderr.write("{}\n".format(tln))
                        raise self.SNP_Error('File {} is not a snp file.'
                                             .format(fl))
                    fin.seek(0)
                    self._c.executemany("INSERT INTO '{}' VALUES (?, ?, ?)"
                                        .format(chrom),
                                        (line.strip().split('\t') for line in fin)
                                       )
                    self.chromosomes[chrom] += 1
                self._conn.commit()

        elif os.path.isfile(snp_file_dir):
            sys.stderr.write('Not implemented yet.\n')
            return
        else:
            raise OSError('File not found {}'.format(snp_file_dir))

        # Create indicies
        self._create_indices()

    def find(self, chromosome, location=None):
        """ Return ref, alt if found, None if not. """
        if chromosome not in self.chromosomes:
            return None

        tables = 'ref,alt' if location else 'pos,ref,alt'
        expr = "SELECT {0} FROM '{1}'".format(tables, chromosome)
        if location:
            expr += (" INDEXED BY '{0}_pos' WHERE pos = {1}"
                     .format(chromosome, location))

        try:
            self._c.execute(expr)
        except sqlite3.OperationalError as e:
            if str(e).startswith('no such table'):
                sys.stderr.write("WARNING --> Chromosome '{}' is not in "
                                 "the lookup table, lookup failed.\n"
                                 .format(chromosome))
                return None
            else:
                sys.stderr.write(expr + '\n')
                raise(e)

        answer = self._c.fetchall()
        if answer:
            return answer[0] if location else answer
        else:
            return None

    #######################
    #  Private Functions  #
    #######################

    def _add_table_if_none(self, table):
        """ Add a table if it does not already exist. """
        expr = ("SELECT * FROM sqlite_master WHERE name = '{}' " +
                "and type='table';").format(table)
        self._c.execute(expr)
        if not self._c.fetchall():
            exp = ("CREATE TABLE '{}' (pos int, ref text, alt text);"
                   .format(table))
            self._c.execute(exp)
            self._conn.commit()

    def _initdb(self):
        """ Ininitialize the database connection. """
        self._conn = sqlite3.connect(self.db)
        self._c    = self._conn.cursor()

    def _create_indices(self):
        """ Add the indicies needed for fast lookups. """
        self._c.execute("SELECT name FROM sqlite_master WHERE type='table';")
        for i in self._c.fetchall():
            exp = ("CREATE INDEX '{0}_pos' ON '{0}' " +
                    "(pos)").format(i[0])
            self._c.execute(exp)
            self._conn.commit()

    #################
    #  Error Class  #
    #################

    class SNP_Error(Exception):

        """ An exception class for SNP parsing issues. """

    ###############
    #  Internals  #
    ###############

    def __getitem__(self, item):
        """Alias for find, syntax: SNPDB[chr,pos]."""
        if isinstance(item, (list, tuple)):
            if len(item) is 2:
                return self.find(item[0], item[1])
            elif len(item) is 1:
                return self.find(item[0])
            else:
                raise self.SNP_Error('list must be chr,pos or chr only')
        elif isinstance(item, str):
            return self.find(item)
        else:
            raise TypeError('{} should be chr or [chr, pos], is {}'
                            .format(item, type(item)))

    def __getattr__(self, attr):
        """Prevent lookup of attributes if already set."""
        if attr == "length":
            return self.length if hasattr(self, "length") else len(self)

    def __len__(self):
        """The total number of SNPs."""
        length = 0
        self._c.execute("SELECT name FROM sqlite_master WHERE type='table';")
        chromosomes = self._c.fetchall()
        for chrom in chromosomes:
            chrom = chrom[0]
            self._c.execute("SELECT Count(*) FROM '{}';".format(chrom))
            l = self._c.fetchone()[0]
            self.chromosomes[chrom] = l
            length += l
        self.length = length
        return length

    def __iter__(self):
        """Iterate over every SNP."""
        for chrom in self.chromosomes:
            print("Iterating over chromosome {}".format(chrom))
            self._c.execute("SELECT * FROM {};".format(chrom))
            for row in self._c:
                pos, ref, alt = row
                yield chrom, pos, ref, alt

    def __repr__(self):
        """Basic information about this class."""
        return "SNPDB<len: {}>".format(self.length)

    def __str__(self):
        """Print chromosome lengths."""
        output = ('SNPDB Object. Total SNPs: {}. Chromosomes::\n'
                  .format(self.length))
        for chrom, len in self.chromosomes.items():
            output += '\tchromosome {}: {} SNPs\n'.format(chrom, len)
        return output


def open_zipped(infile, mode='r'):
    """ Return file handle of file regardless of zipped or not
        Text mode enforced for compatibility with python2 """
    mode   = mode[0] + 't'
    p2mode = mode
    if hasattr(infile, 'write'):
        return infile
    if isinstance(infile, str):
        if infile.endswith('.gz'):
            return gzip.open(infile, mode)
        if infile.endswith('.bz2'):
            if hasattr(bz2, 'open'):
                return bz2.open(infile, mode)
            else:
                return bz2.BZ2File(infile, p2mode)
        return open(infile, p2mode)


def product(iterable):
    "Returns the product of all items in the iterable"
    return reduce(mul, iterable, 1)

def get_indels(snps):
    """Returns a dict-of-dicts with positions of indels."""
    indel_dict = defaultdict(dict)
    for chrom, pos, ref, alt in snps:
        if ('-' in [ref, alt]) or (max(len(i) for i in [ref, alt]) > 1):
            indel_dict[chrom][pos] = True
    return indel_dict

RC_TABLE = {
    ord('A'):ord('T'),
    ord('T'):ord('A'),
    ord('C'):ord('G'),
    ord('G'):ord('C'),
}

def reverse_complement(seq):
    "Reverse complements the string input"
    return seq.translate(RC_TABLE)[::-1]


def get_dual_read_seqs(read1, read2, snps_db, indel_dict, dispositions):
    """ For each pair of reads, get all concordant SNP substitutions

    Note that if the reads overlap, the matching positions in read1 and read2
    will get the same subsitution as each other.
    """
    if read1.is_unmapped or read2.is_unmapped:
        dispositions['unmapped read'] += 1
        return [[], []]
    seq1 = read1.seq
    seq2 = read2.seq
    seqs1, seqs2 = [read1.seq], [read2.seq]

    chrom = read1.reference_name
    snps = {}
    read_posns = defaultdict(lambda: [None, None])

    for (read_pos1, ref_pos) in read1.get_aligned_pairs(matches_only=True):
        if indel_dict[chrom].get(ref_pos, False):
            dispositions['toss_indel'] += 1
            return [[], []]
        allele_info = snps_db[chrom, ref_pos]
        if allele_info:
            snps[ref_pos] = allele_info
            read_posns[ref_pos][0] = read_pos1

    for (read_pos2, ref_pos) in read2.get_aligned_pairs(matches_only=True):
        if indel_dict[chrom].get(ref_pos, False):
            dispositions['toss_indel'] += 1
            return [[], []]
        allele_info = snps_db[chrom, ref_pos]
        if allele_info:
            snps[ref_pos] = allele_info
            read_posns[ref_pos][1] = read_pos2

    if product(len(i) for i in snps.values()) > MAX_SEQS_PER_READ:
        dispositions['toss_manysnps'] += 1
        return [[], []]

    for ref_pos in snps:
        alleles = snps[ref_pos]
        pos1, pos2 = read_posns[ref_pos]
        new_seqs1 = []
        new_seqs2 = []
        if pos1 is None:
            for allele in alleles:
                if allele == seq2[pos2]:
                    continue
                for seq1, seq2 in zip(seqs1, seqs2):
                    new_seqs1.append(seq1)
                    new_seqs2.append(''.join([seq2[:pos2], allele, seq2[pos2+1:]]))

        elif pos2 is None:
            for allele in alleles:
                if allele == seq1[pos1]:
                    continue
                for seq1, seq2 in zip(seqs1, seqs2):
                    new_seqs1.append(''.join([seq1[:pos1], allele, seq1[pos1+1:]]))
                    new_seqs2.append(seq2)
        else:
            if seq1[pos1] != seq2[pos2]:
                dispositions['toss_anomalous_phase'] += 1
                return [[], []]
            for allele in alleles:
                if allele == seq2[pos2]:
                    continue
                for seq1, seq2 in zip(seqs1, seqs2):
                    new_seqs1.append(''.join([seq1[:pos1], allele, seq1[pos1+1:]]))
                    new_seqs2.append(''.join([seq2[:pos2], allele, seq2[pos2+1:]]))
        seqs1.extend(new_seqs1)
        seqs2.extend(new_seqs2)

    if len(seqs1) == 1:
        dispositions['no_snps'] += 1
    else:
        dispositions['has_snps'] += 1
    return seqs1, seqs2

def get_read_seqs(read, snps, indel_dict, dispositions):
    """ For each read, get all possible SNP substitutions

    for N biallelic snps in the read, will return 2^N reads
    """
    num_snps = 0
    seqs = [read.seq]

    chrom = read.reference_name
    for (read_pos, ref_pos) in read.get_aligned_pairs(matches_only=True):
        if ref_pos is None:
            continue
        if indel_dict[chrom].get(ref_pos, False):
            dispositions['toss_indel'] += 1
            return []
        if len(seqs) > MAX_SEQS_PER_READ:
            dispositions['toss_manysnps'] += 1
            return []

        alleles = snps[chrom, ref_pos]
        if alleles:
            read_base = read.seq[read_pos]
            if read_base in alleles:
                dispositions['ref_match'] += 1
                num_snps += 1
                for new_allele in alleles:
                    if new_allele == read_base:
                        continue
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

def assign_reads(insam, snps, indel_dict, is_paired=True):
    """ Loop through all the reads in insam and output them to the appropriate file


    """
    fname = insam.filename
    if isinstance(fname, bytes):
        fname = fname.decode('ascii')
    basename = fname.rsplit('.', 1)[0]
    keep = Samfile('.'.join([basename, 'keep.bam']),
                   'wb',
                   template=insam)
    remap_bam = Samfile('.'.join([basename, 'to.remap.bam']),
                        'wb',
                        template=insam)
    dropped_bam = Samfile('.'.join([basename, 'dropped.bam']),
                          'wb',
                          template=insam)
    if is_paired:
        fastqs = [
            gzip.open('.'.join([basename, 'remap.fq1.gz']), 'wt'),
            gzip.open('.'.join([basename, 'remap.fq2.gz']), 'wt'),
        ]
    else:
        fastqs = [gzip.open('.'.join([basename, 'remap.fq.gz']), 'wt'),]
    unpaired_reads = [{}, {}]
    read_results = Counter()
    remap_num = 1
    for i, read in enumerate(insam):
        if i % 10000 == 0:
            pass
        if not is_paired:
            read_seqs = get_read_seqs(read, snps, indel_dict, read_results)
            write_read_seqs([(read, read_seqs)], keep, remap_bam, fastqs)
        elif read.is_proper_pair:
            slot_self = read.is_read2 # 0 if is_read1, 1 if read2
            slot_other = read.is_read1
            if read.qname in unpaired_reads[slot_other]:
                both_reads = [None, None]
                both_reads[slot_self] = read
                both_reads[slot_other] = unpaired_reads[slot_other].pop(read.qname)
                both_seqs = get_dual_read_seqs(both_reads[0], both_reads[1],
                                               snps, indel_dict, read_results)
                both_read_seqs = list(zip(both_reads, both_seqs))
                remap_num += write_read_seqs(both_read_seqs, keep, remap_bam,
                                             fastqs, dropped_bam, remap_num)
            else:
                unpaired_reads[slot_self][read.qname] = read
        else:
            read_results['not_proper_pair'] += 1
            # Most tools assume reads are paired and do not check IDs. Drop it out.
            continue
    print()
    print(len(unpaired_reads[0]), len(unpaired_reads[1]))
    print(read_results)


def write_read_seqs(both_read_seqs, keep, remap_bam, fastqs, dropped=None, remap_num=0):
    """Write the given reads out to the appropriate file


    If there are no SNPs in the read, write it to the BAM file `keep`

    If there are too many SNPs, and `dropped` is provided, write the original
    read out to `dropped`

    Otherwise, output all possible substitutions to the fastqs for remapping,
    as well as a bam file containing the original read.
    """
    reads, seqs = zip(*both_read_seqs)
    assert len(reads) == len(fastqs)

    num_seqs = product(len(r[1]) for r in both_read_seqs)
    if num_seqs == 0 or num_seqs > MAX_SEQS_PER_READ:
        if dropped is not None:
            for read in reads:
                dropped.write(read)
            return 0
        else:
            return 0
    elif num_seqs == 1:
        for read, seqs in both_read_seqs:
            keep.write(read)
    else:
        assert len(reads) > 0
        for read in reads:
            remap_bam.write(read)
        left_pos = min(r.pos for r in reads)
        right_pos = max(r.pos for r in reads)
        loc_line = '{}:{}:{}:{}:{}'.format(
            remap_num,
            read.reference_name,
            left_pos,
            right_pos,
            num_seqs-1,
        )

        if left_pos == 16053407:
            print(seqs)
        first = True
        # Some python fanciness to deal with single or paired end reads (or
        # n-ended reads, if such technology ever happens.
        for read_seqs in it.product(*seqs):
            if first:
                first = False
                continue
            for seq, read, fastq in zip(read_seqs, reads, fastqs):
                fastq.write(
                    "@{loc_line}\n{seq}\n+{loc_line}\n{qual}\n"
                    .format(
                        loc_line=loc_line,
                        seq=reverse_complement(seq) if read.is_reverse else seq,
                        qual=read.qual)
                    )
        return 1
    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--paired_end",
                        action='store_true',
                        dest='is_paired_end', default=False,
                        help=('Indicates that reads are '
                              'paired-end (default is single).'))

    parser.add_argument("-s", "--sorted",
                        action='store_true', dest='is_sorted', default=False,
                        help=('Indicates that the input bam file'
                              ' is coordinate sorted (default is False)'))


    parser.add_argument("infile", type=Samfile, help=("Coordinate sorted bam "
                                                      "file."))
    snp_dir_help = ('Directory containing the SNPs segregating within the '
                    'sample in question (which need to be checked for '
                    'mappability issues).  This directory should contain '
                    'sorted files of SNPs separated by chromosome and named: '
                    'chr<#>.snps.txt.gz. These files should contain 3 columns: '
                    'position RefAllele AltAllele')

    parser.add_argument("snp_dir", action='store', help=snp_dir_help)

    options = parser.parse_args()

    global SNPS
    snps       = SNPDB(options.snp_dir)
    print("Finding indels")
    indel_dict = get_indels(snps)

    print("Done with SNPs")

    assign_reads(options.infile, snps, indel_dict, options.is_paired_end)
