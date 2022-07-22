#!/usr/bin/env python


__author__ = 'duceppemo'
__version__ = 'v0.1'


import os
from argparse import ArgumentParser
from multiprocessing import cpu_count
from time import time
from seqcounter_methods import SeqCounterMethods


class SeqCounter(object):
    def __init__(self, args):
        # I/O
        self.input_folder = args.input
        self.ref = args.ref
        self.r1 = args.r1
        self.r2 = args.r2
        self.output = args.output
        self.kmer_size = args.kmer_size
        self.mismatch = args.mismatch

        # Performance
        self.cpu = args.threads
        self.parallel = args.parallel

        # Run
        self.run()

    def run(self):
        start_time = time()

        # Check input
        SeqCounterMethods.check_input(self.input_folder, self.r1, self.r2)

        # Get sample(s)
        sample_dict = SeqCounterMethods.get_samples(self.input_folder, self.r1, self.r2)

        # Count spacers
        refseq_count_dict = SeqCounterMethods.count_refseqs_seal_parallel(self.ref, sample_dict, self.output,
                                                                          self.kmer_size, self.mismatch,
                                                                          self.cpu, self.parallel)

        # Print report on screen and to file
        SeqCounterMethods.print_report(refseq_count_dict, self.output)

        # Print elapsed time
        print('Elapsed time: {}'.format(SeqCounterMethods.elapsed_time(time() - start_time)))


if __name__ == "__main__":
    max_cpu = cpu_count()

    parser = ArgumentParser(description='Count occurrence of reference sequence(s) in fastq or fasta.')
    parser.add_argument('-i', '--input', metavar='/path/to/folder/with/fastq/or/fasta/files',
                        required=False, type=str,
                        help='Folder containing fastq or fasta file to search. '
                             'Gzipped or not. Mandatory.')
    parser.add_argument('-r1', metavar='/path/to/sample_R1.[fastq|fasta]',
                        required=False, type=str,
                        help='R1 fastq file from paired-end or single end sequencing data or fasta file. '
                             'Gzipped or not. Mandatory.')
    parser.add_argument('-r2', metavar='/path/to/R2/fastq',
                        required=False, type=str,
                        help='R2 fastq file from paired-end. Optional.')
    parser.add_argument('-ref', metavar='/path/to/ref.fasta',
                        required=True, type=str,
                        help='Fasta file containing reference sequence(s) to find and count occurrence. Optional.')
    parser.add_argument('-k', '--kmer-size', metavar='25',
                        required=False, default=25, type=int,
                        help='Kmer length used for finding reference sequences.  reference sequences shorter than k '
                             'will not be found. Default 25. Optional.')
    parser.add_argument('-m', '--mismatch', metavar='1',
                        required=False, default=1, type=int,
                        help='Maximum number of mismatches between the reference and the input sequences '
                             '(Maximum Hamming distance for ref kmers). Default 1. Optional.')
    parser.add_argument('-o', '--output', metavar='/path/to/output/folder/',
                        required=True, type=str,
                        help='Folder to hold the result files. Mandatory.')
    parser.add_argument('-t', '--threads', metavar=str(max_cpu),
                        required=False,
                        type=int, default=max_cpu,
                        help='Number of threads. Default is maximum available({}). Optional'.format(max_cpu))
    parser.add_argument('-p', '--parallel', metavar='2',
                        required=False,
                        type=int, default=2,
                        help='Number of samples to process in parallel. Default is 2. Optional.')
    parser.add_argument('-v', '--version', action='version',
                        version=f'{os.path.abspath(__file__)}: version {__version__}')

    # Get the arguments into an object
    arguments = parser.parse_args()

    SeqCounter(arguments)
