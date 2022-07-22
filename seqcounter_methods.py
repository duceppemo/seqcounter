import os
import sys
from collections import OrderedDict
import subprocess
import pathlib
from multiprocessing import cpu_count
from psutil import virtual_memory
from concurrent import futures
import pandas as pd


class SeqCounterMethods(object):
    ext = ['.fna', '.fna.gz',
           '.fa', '.fa.gz',
           '.fasta', '.fasta.gz',
           '.fastq', '.fastq.gz',
           '.fq', '.fq.gz']

    @staticmethod
    def check_input(input_folder, r1, r2):
        if not input_folder and not r1:
            raise Exception('Please use "--input" OR "-r1" (and "-r2") as input.')

        if input_folder and (r1 or r2):
            raise Exception('Please use "--input" OR "-r1" (and "-r2") as input.')

        if r2 and not r1:
            raise Exception('Please specify "-r1" before specifying "-r2"')

        if input_folder:
            if not os.path.isdir(input_folder):
                raise Exception('Please choose a folder for "--input".')
            if not os.path.exists(input_folder):
                raise Exception('Please choose an existing folder for "--input".')

    @staticmethod
    def check_cpus(requested_cpu, n_proc):
        total_cpu = cpu_count()

        # Number of threads
        if requested_cpu > total_cpu:
            requested_cpu = total_cpu
            sys.stderr.write("Requested threads was set higher than available system threads ({})".format(total_cpu))
            sys.stderr.write("Number of threads was set to {}".format(requested_cpu))
        elif requested_cpu < 1:
            requested_cpu = 1
            sys.stderr.write("Number of threads was set to {}".format(1))

        # Number of parallel processes
        if n_proc > requested_cpu:
            n_proc = total_cpu
            sys.stderr.write("Requested parallel processes was set higher than requested threads ({})".format(
                requested_cpu))
            sys.stderr.write("Number of parallel processes was set to {}".format(total_cpu))
        elif n_proc < 1:
            n_proc = 1
            sys.stderr.write("Number of threads was parallel processes to {}".format(1))

        return requested_cpu, n_proc

    @staticmethod
    def check_mem(requested_mem):
        max_mem = int(virtual_memory().total * 0.85 / 1000000000)  # in GB
        if requested_mem:
            if requested_mem > max_mem:
                requested_mem = max_mem
                sys.stderr.write("Requested memory was set higher than available system memory ({})".format(max_mem))
                sys.stderr.write("Memory was set to {}".format(requested_mem))
        else:
            requested_mem = max_mem

        return requested_mem

    @staticmethod
    def get_sample_name(input_file):
        sample = '.'.join(os.path.basename(input_file).split('.')[:-1])  # basename minus what is after last "."
        if input_file.endswith('gz'):  # Need to drop what after the last 2 "."
            sample = '.'.join(sample.split('.')[:-1])
        if any([x in sample for x in ['_R1', '_1']]) or any([x in sample for x in ['_R2', '_2']]):
            sample = '_'.join(sample.split("_")[:-1])
        return sample

    @staticmethod
    def get_samples(input_folder, r1, r2):
        sample_dict = dict()

        if input_folder:
            # Walk the input directory recursively and look for fastq files
            for root, directories, filenames in os.walk(input_folder):
                for filename in filenames:
                    absolute_path = os.path.join(root, filename)
                    if os.path.isfile(absolute_path) and filename.endswith(tuple(SeqCounterMethods.ext)):
                        # Get sample name
                        sample = SeqCounterMethods.get_sample_name(absolute_path)

                        # Add to sample dictionary
                        if sample not in sample_dict:
                            sample_dict[sample] = []
                        if any([x in filename for x in ['_R1', '_1']]):
                            sample_dict[sample].insert(0, absolute_path)
                        elif any([x in filename for x in ['_R2', '_2']]):
                            sample_dict[sample].insert(1, absolute_path)
                        else:
                            sample_dict[sample].insert(0, absolute_path)

        elif r1:
            sample = SeqCounterMethods.get_sample_name(r1)
            sample_dict[sample] = [os.path.abspath(r1)]
            if r2:
                sample_dict[sample].append(os.path.abspath(r1))
        else:
            raise Exception('Please provide valid input. See help.')

        return sample_dict

    @staticmethod
    def make_folder(folder):
        # Will create parent directories if don't exist and will not return error if already exists
        pathlib.Path(folder).mkdir(parents=True, exist_ok=True)

    @staticmethod
    def count_refseqs_seal(ref, output_folder, fastq_list, sample, kmer_size, mismatch, cpu):
        print('Processing {}'.format(sample))

        # Stats file
        stats_file = output_folder + '/' + sample + '_stats.tsv'

        if len(fastq_list) == 1:
            cmd = ['seal.sh',
                   'in={}'.format(fastq_list[0]),
                   'ref={}'.format(ref),
                   'k={}'.format(kmer_size),
                   'rcomp=t',
                   'hdist={}'.format(mismatch),
                   'clearzone=999999',
                   'ambiguous=all',
                   'maskmiddle=f',  # Do not treat the middle base of a kmer as a wildcard
                   'ow=t',
                   'stats={}'.format(stats_file),
                   'nzo=f',
                   'threads={}'.format(cpu)]
        else:
            cmd = ['seal.sh',
                   'in={}'.format(fastq_list[0]),
                   'in2={}'.format(fastq_list[1]),
                   'ref={}'.format(ref),
                   'k={}'.format(kmer_size),
                   'rcomp=t',
                   'hdist={}'.format(mismatch),
                   'clearzone=999999',
                   'ambiguous=all',
                   'maskmiddle=f',  # Do not treat the middle base of a kmer as a wildcard
                   'ow=t',
                   'stats={}'.format(stats_file),
                   'nzo=f',
                   'threads={}'.format(cpu)]

        subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

        # Parse stats file
        """
        #File /home/bioinfo/analyses/vsnp3_test_spoligo/SRR16058435_R1.fastq.gz
        #Total  822714
        #Matched    799 0.09712%
        #Name   Reads   ReadsPct
        spacer25    62  0.00754%
        spacer02    48  0.00583%
        """

        refseq_count_dict = dict()

        with open(stats_file, 'r') as f:
            for line in f:
                line = line.rstrip()
                if not line:
                    continue

                if line.startswith('#'):
                    continue
                else:
                    field_list = line.split('\t')
                    refseq_count_dict[field_list[0]] = int(field_list[1])

        # Remove temporary stat file
        os.remove(stats_file)

        return sample, OrderedDict(sorted(refseq_count_dict.items()))  # Order dictionary by keys

    @staticmethod
    def count_refseqs_seal_parallel(ref, sample_dict, output_folder, kmer_size, mismatch, cpu, parallel):
        SeqCounterMethods.make_folder(output_folder)
        sample_count_dict = dict()

        with futures.ThreadPoolExecutor(max_workers=int(parallel)) as executor:
            # ref, output_folder, fastq_list, sample, cpu
            args = ((ref, output_folder, fastq_list, sample, kmer_size, mismatch, int(cpu / parallel))
                    for sample, fastq_list in sample_dict.items())
            for sample, count_dict in executor.map(lambda x: SeqCounterMethods.count_refseqs_seal(*x), args):
                sample_count_dict[sample] = count_dict

        return sample_count_dict

    @staticmethod
    def print_report(refseq_count_dict, output_folder):
        report_file = output_folder + '/' + 'refseqcounts.tsv'

        sample_list = list()
        df_list = list()
        for sample, count_dict in refseq_count_dict.items():
            sample_list.append(sample)
            df_list.append(pd.DataFrame.from_dict(count_dict, orient='index'))

        df = pd.concat(df_list, axis='columns', names=sample_list)

        # Rename columns
        df.columns = sample_list
        df.index.name = 'Reference'

        df.to_csv(report_file, sep='\t', header=True, index=True)

        # Add command line used for QA
        command_line = ' '.join(sys.argv)

        with open(report_file, 'r+') as f:
            # Slurp file content
            file_content = f.read()

            # Go to start of file
            f.seek(0)

            # Prepend file content with command line used
            f.write('#{}\n{}'.format(command_line, file_content))

    @staticmethod
    def elapsed_time(seconds):
        """
        Transform a time value into a string
        :param seconds: Elapsed time in seconds
        :return: Formatted time string
        """
        minutes, seconds = divmod(seconds, 60)
        hours, minutes = divmod(minutes, 60)
        days, hours = divmod(hours, 24)
        periods = [('d', days), ('h', hours), ('m', minutes), ('s', seconds)]
        time_string = ''.join('{}{}'.format(int(round(value)), name) for name, value in periods if value)
        return time_string
