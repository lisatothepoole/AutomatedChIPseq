#!/usr/bin/python

"""
ChIP-seq data analysis for Cortez Lab
Authors: Lisa Poole and James Pino

Requirements for this experiment - sequencing data files must be in a folder
entitled "fastq" within a folder corresponding to the experiment name;
modify general information in the first part of the script
(starting with sample base).

Takes zipped fasta/q files and results in sorted bam files that are indexed to be viewed via a genome browser and an Excel file of peaks with a BED file of peak summits.
"""

import os
import subprocess

# General information - modify appropriately for each experiment

number_of_samples = '4'  # Number of sequencing samples
species = 'human'
read_length = 75
sample_suffix = 'fastq'
compression_suffix = 'gz'
read_type = 'SE'  # Valid options are SE (single-end) or PE (paired-end)
n_cpus = 8
pc = 'lisa'
experiment_name = 'E76b_RPA_ChIP-seq'
q_value_cutoff = '0.01'  # Default is 0.01, use 0.05 for broad peaks
duplication_choice = 'all'  # Valid options are 'auto' to have macs2 determine the max number of reads possible at one site, 'all' to keep all tags at the same location, or an integer value selected by the user
peak_type = 'normal'  # Valid choices are 'broad' or 'normal'(narrow)

if pc == 'cortez_mac':
    # identifying path of programs
    samtools = '/usr/local/bin/samtools'
    samstat = '/usr/local/bin/samstat'
    bwa = '/usr/local/bin/bwa-0.7.15/bwa'
    picard = '/Users/temporary/Sources/picard.jar'
    homer = '/Users/temporary/Sources/homer/bin/'
    fastqc = '/Users/temporary/Sources/FastQC.app/Contents/MacOS/fastqc'
    flexbar = '/usr/local/bin/flexbar_v2.5_macosx'
    macs2 = '/Users/lisapoole/Sources/MACS2-2.1.1.20160309/macs2'
    adaptors = '/Users/temporary/genomes/adapters/illumina_truseq.fasta'

    # experiment specific information
    output_directory = "/Users/temporary/projects/{}".format(experiment_name)
    fasta_directory = "/Users/temporary/projects/{}/fastq".format(experiment_name)

    # Reference files
    if species == 'mouse':
        # This protocol is adapted for the mm10 version of the mouse genome
        bwa_index = '/Users/temporary/genomes/Mus_musculus/UCSC/mm10/Sequence/BWAIndex/version0.7.15/genome_indel'
        reference_genome = '/Users/temporary/genomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome_indel.fa'

    elif species == 'human':
        # This pipeline uses the hg38 version of the human genome
        bwa_index = '/Users/temporary/genomes/Homo_sapiens_hg38/Homo_sapiens/UCSC/hg38/Sequence/BWAIndex/genome'
        reference_genome = '/Users/temporary/genomes/Homo_sapiens_hg38/Homo_sapiens/UCSC/hg38/Sequence/BWAIndex/genome.fa'

    else:
        print(
            "Error - invalid species. Valid arguments are 'human' or 'mouse'")
        quit()

elif pc == 'lisa':
    output_directory = '/Users/lisapoole/Desktop/E76b_asisi_rpa_chip-seq'
    # Setting up programs
    macs2 = '/Users/lisapoole/Sources/MACS2-2.1.1.20160309/bin/macs2'

if species == 'human':
    genome_size = '2.7e9'
elif species == 'mouse':
    genome_size = '1.87e9'

def fastqc_analysis(sample_base):
    # this step analyzes the fastq file obtained from sequencing to get some preliminary quality analyses
    print("Starting fastq analysis of {}".format(sample_base))
    if not os.path.exists('{}/quality_control'.format(output_directory)):
        os.mkdir('{}/quality_control'.format(output_directory))

    path_to_executable = fastqc
    input_files = '{}/{}.{}.{}'.format(fasta_directory, sample_base, sample_suffix, compression_suffix)
    output = '-o {}/quality_control'.format(output_directory)
    important_options = '-f {} -t {}'.format(sample_suffix,
                                             n_cpus)  # -t indicates number of threads/cpus -f specifies the format of the file
    command = [path_to_executable, output, important_options, input_files]
    call_code = ' '.join(command)
    print(call_code)
    process = subprocess.Popen([call_code], shell=True,
                               stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    while True:
        output = process.stdout.readline()
        if output == '' and process.poll() is not None:
            break
            # if output:
            #         # global_output += output
            # print output.strip()
        rc = process.poll()
    print("Done analyzing fastq for {}".format(sample_base))


# Adaptor Trimming via Flexbar
def flexbar_trim(sample_base):
    # Flexbar options
    # '-a' = adaptor sequences (fasta format)
    # '-n' = number of threads
    # '-u' = max uncalled bases for each read to pass filtering
    # '-m' = min read length to remain after filtering/trimming
    # '-t' = prefix for output file names
    # '-ao' = adapter min overlap
    # '-ae' = adapter trim end
    print("Start trimming {}".format(sample_base))
    path_to_executable = flexbar
    suffix_for_output = '-t {}/{}-trimmed'.format(fasta_directory, sample_base)
    adaptor_trim_end = '-ae ANY'
    adaptor_overlap = '-ao 5'
    path_to_adaptors = "-a {}".format(adaptors)
    number_max_uncalled_bases_pass = '-u {}'.format(read_length)
    main_read_length_to_remain = '-m 18'
    threads = '-n {}'.format(n_cpus)

    if read_type not in ('SE', 'PE'):
        print("Error. Invalid read type. Valid options are SE or PE. Aborting.")
        quit()

    elif read_type == 'SE':
        reads = '-r {}/{}.{}.{}'.format(fasta_directory, sample_base, sample_suffix, compression_suffix)

    elif read_type == 'PE':
        reads = ' -r {}/{}_1.{}.{} -p {}/{}_2.{}.{}'.format(fasta_directory, sample_base, sample_suffix,
                                                            compression_suffix, fasta_directory, sample_base,
                                                            sample_suffix, compression_suffix)

    command = [path_to_executable, reads, path_to_adaptors, threads, suffix_for_output, adaptor_overlap,
               adaptor_trim_end, number_max_uncalled_bases_pass, main_read_length_to_remain]
    call_code = ' '.join(command)
    print(call_code)
    process = subprocess.Popen([call_code], shell=True,
                               stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    while True:
        output = process.stdout.readline()
        if output == '' and process.poll() is not None:
            break
        if output:
            print output.strip()
        rc = process.poll()
    print("Done trimming {}".format(sample_base))


def bwa_alignment(sample_base):
    print("Starting alignment for {}".format(sample_base))
    if not os.path.exists('{}/BWA_BAM_files'.format(output_directory)):
        os.mkdir('{}/BWA_BAM_files'.format(output_directory))

    path_to_executable = '{} mem'.format(bwa)

    if read_type == 'PE':
        reads = "{0}/{1}-trimmed_1.{2} {0}/{1}-trimmed_2.{2}".format(fasta_directory, sample_base, sample_suffix)
    elif read_type == 'SE':
        reads = '{}/{}-trimmed.{}'.format(fasta_directory, sample_base, sample_suffix)
    else:
        print("Error - invalid read type. Valid arguments are 'SE' or 'PE'")
        quit()

    important_options = '-t {}'.format(n_cpus)  # -t indicates number of threads/cpus
    path_to_reference = reference_genome
    export_to_file = '> {}/BWA_BAM_files/{}.sam'.format(output_directory, sample_base)
    command = [path_to_executable, important_options, path_to_reference, reads, export_to_file]
    call_code = ' '.join(command)
    print(call_code)
    process = subprocess.Popen([call_code], shell=True,
                               stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    while True:
        output = process.stdout.readline()
        if output == '' and process.poll() is not None:
            break
            # if output:
            #         # global_output += output
            # print output.strip()
        rc = process.poll()
    print("Done aligning {}".format(sample_base))


def sam_read_group_addition(sample_base):
    print("Starting read group addition for {}".format(sample_base))
    path_to_executable = 'java -jar {}'.format(picard)
    picard_program = "AddOrReplaceReadGroups"
    input_files = 'I={}/BWA_BAM_files/{}.sam'.format(output_directory, sample_base)
    output_files = 'O={}/BWA_BAM_files/{}.rg.sam'.format(output_directory, sample_base)
    necessary_parameters = "RGID={0} RGLB={0} RGPL=ILLUMINA RGPU=ILLUMINA RGSM={0}".format(sample_base)
    command = [path_to_executable, picard_program, input_files, output_files, necessary_parameters]
    call_code = ' '.join(command)
    print(call_code)
    process = subprocess.Popen([call_code], shell=True,
                               stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    while True:
        output = process.stdout.readline()
        if output == '' and process.poll() is not None:
            break
        if output:
            print output.strip()
        rc = process.poll()
    print("Done adding read group {}".format(sample_base))


def sam_to_bam(sample_base):
    print("Start sam to bam conversion {}".format(sample_base))
    path_to_executable = '{} view'.format(samtools)
    path_to_samples = '-S -b {}/BWA_BAM_files/{}.sam'.format(output_directory,
                                                             sample_base)  # -S indicates input is SAM file -b indicates output will be BAM
    output_filename = '-o {}/BWA_BAM_files/{}.bam'.format(output_directory, sample_base)
    threads = '--threads {}'.format(n_cpus)
    command = [path_to_executable, path_to_samples, threads, output_filename]
    call_code = ' '.join(command)
    print(call_code)
    process = subprocess.Popen([call_code], shell=True,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT)
    while True:
        output = process.stdout.readline()
        if output == '' and process.poll() is not None:
            break
        if output:
            print output.strip()
        rc = process.poll()
    print("Done sam to bam conversion {}".format(sample_base))


def bam_sort(sample_base):
    print("Start sorting {}".format(sample_base))
    path_to_executable = '{} sort'.format(samtools)
    path_to_samples = '{}/BWA_BAM_files/{}.bam'.format(output_directory, sample_base)
    output_filename = '-o {}/BWA_BAM_files/{}.sorted.bam'.format(output_directory, sample_base)
    threads = '--threads {}'.format(n_cpus)
    command = [path_to_executable, threads, output_filename, path_to_samples]
    call_code = ' '.join(command)
    print(call_code)
    process = subprocess.Popen([call_code], shell=True,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT)
    while True:
        output = process.stdout.readline()
        if output == '' and process.poll() is not None:
            break
        if output:
            print output.strip()
        rc = process.poll()
    print("Done sorting {}".format(sample_base))


def bam_index(sample_base):
    print("Start indexing {}".format(sample_base))
    path_to_executable = '{} index'.format(samtools)
    path_to_samples = '{}/BWA_BAM_files/{}.sorted.bam'.format(output_directory, sample_base)
    command = [path_to_executable, path_to_samples]
    call_code = ' '.join(command)
    print(call_code)
    process = subprocess.Popen([call_code], shell=True,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT)
    while True:
        output = process.stdout.readline()
        if output == '' and process.poll() is not None:
            break
        if output:
            print output.strip()
        rc = process.poll()
    print("Done indexing {}".format(sample_base))


def unaligned_reads_removal(sample_base):
    print("Removing reads that didn't align or aligned to multiple places in the genome")
    path_to_executable = '{} view'.format(samtools)
    options = '-bF 4'
    input_file = '{}/BWA_BAM_files/{}.sorted.bam'.format(output_directory, sample_base)
    output_file = '> {}/BWA_BAM_files/{}.filtered.bam'.format(output_directory, sample_base)
    command = [path_to_executable, options, input_file, output_file]
    call_code = ' '.join(command)
    print(call_code)
    process = subprocess.Popen([call_code], shell=True,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT)
    while True:
        output = process.stdout.readline()
        if output == '' and process.poll() is not None:
            break
        if output:
            print output.strip()
        rc = process.poll()
    print("Finished removing reads that didn't align or aligned to multiple places in the genome")


def samstat_analysis(sample_base):
    print("Start SAMSTAT check {}".format(sample_base))

    path_to_executable = samstat
    path_to_samples = '{}/BWA_BAM_files/{}.sorted.bam'.format(output_directory, sample_base)
    command = [path_to_executable, path_to_samples]
    call_code = ' '.join(command)
    print(call_code)
    process = subprocess.Popen([call_code], shell=True,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT)
    while True:
        output = process.stdout.readline()
        if output == '' and process.poll() is not None:
            break
        if output:
            print output.strip()
        rc = process.poll()
    print("Done SAMSTAT check {}".format(sample_base))
    os.rename('{}/BWA_BAM_files/{}.sorted.bam.samstat.html'.format(output_directory, sample_base),
              '{}/quality_control/{}.sorted.bam.samstat.html'.format(output_directory, sample_base))


def macs2_peak_calling(sample_base):
    print("Begin calling peaks for {}".format(sample_base))
    if not os.path.exists('{}/peaks'.format(output_directory)):
        os.mkdir('{}/peaks'.format(output_directory))
    path_to_executable = '{} callpeak'.format(macs2)
    treatment_file = '-t {}/{}_rpa_ip.sorted.bam'.format(output_directory, sample_base)
    input_file = '-c {}/{}_input.sorted.bam'.format(output_directory, sample_base)
    input_format = '-f BAM'
    duplicates = '--keep-dup {}'.format(duplication_choice)
    genome = '-g {}'.format(genome_size)
    naming = '-n {}'.format(sample_base)
    output_dir = '--outdir {}/peaks'.format(output_directory)
    enrichment = '--qvalue {}'.format(q_value_cutoff)
    if peak_type == 'broad':
        peak_indicator = '--broad'
    elif peak_type == 'normal':
        peak_indicator = ''
    else:
        print('Invalid peak_type choice.')
    command = [path_to_executable, treatment_file, input_file, input_format, duplicates, enrichment, genome, naming,
               peak_indicator, output_dir]
    call_code = ' '.join(command)
    print(call_code)
    process = subprocess.Popen([call_code], shell=True,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT)
    while True:
        output = process.stdout.readline()
        if output == '' and process.poll() is not None:
            break
        if output:
            print output.strip()
        rc = process.poll()
    print("Done calling peaks for {}".format(sample_base))


def excess_file_clean_up(sample_base):
    os.remove('{}/BWA_BAM_files/{}.sam'.format(output_directory, sample_base))
    os.remove('{}/BAM_files/{}.rg.sam'.format(output_directory, sample_base))
    os.remove('./BAM_files/{}.bam'.format(sample_base))


def automated_chip_seq_analysis(sample_base):
    fastqc_analysis(sample_base)
    flexbar_trim(sample_base)
    bwa_alignment(sample_base)
    sam_read_group_addition(sample_base)
    sam_to_bam(sample_base)
    bam_sort(sample_base)
    bam_index(sample_base)
    samstat_analysis(sample_base)
    excess_file_clean_up(sample_base)

macs2_peak_calling('asisi_etoh')