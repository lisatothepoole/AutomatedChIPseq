#!/usr/bin/python2
"""
This pipeline is for analyzing ChIP-seq peaks that were obtained using the automated_chip-seq_analysis.py

"""

import pandas as pd


#  Experiment specific
species = 'human'  # Valid options are 'human' or 'mouse'
read_length = 75
read_type = 'SE'  # Valid options are 'SE' (single-end) or 'PE' (paired-end)
pc = 'cortez_mac'
experiment_name = 'E76_RPA_ChIP-seq'
output_directory = '/Users/lisapoole/Desktop/{}'.format(experiment_name)
pc = 'lisa'
cutoff_of_reads = '90'

# Step 1: import blacklist files
if pc == 'lisa':
    if species == 'human':
        blacklist_file = '/Users/lisapoole/Sources/hg38.blacklist.bed'
        repeat_file = '/Users/lisapoole/Sources/hg38_repeat_telomere_centromere.txt'
    elif species == 'mouse':
        blacklist_file = '/Users/lisapoole/Sources/mm10.blacklist.bed'
    else:
        print('Species value invalid.')

elif pc == 'cortez_mac':
    if species == 'human':
        blacklist_file = '/Users/temporary/genomes/Homo_sapiens_hg38/Homo_sapiens/UCSC/hg38/Sequence/AbundantSequences/hg38.blacklist.bed'
        repeat_file = '/Users/temporary/genomes/Homo_sapiens_hg38/Homo_sapiens/UCSC/hg38/Sequence/AbundantSequences/hg38_repeat_telomere_centromere.txt'
    elif species == 'mouse':
        blacklist_file = '/Users/temporary/genomes/Mus_musculus/UCSC/mm10/Sequence/AbundantSequences/mm10.blacklist.bed'
    else:
        print('Species value invalid.')

else:
    print('Invalid pc choice')

if species == 'human':
    genome_size = 'hs'  # indicates genome size will be 2.7e9
elif species == 'mouse':
    genome_size = 'mm'  # indicates genome size will be 1.87e9
else:
    print('Species value invalid.')

# Indexing blacklist and hg38 repeat files
blacklist = pd.read_table(blacklist_file, usecols=[0, 1, 2], names=['chr', 'start', 'end'])
blacklist['index_name'] = blacklist['chr'].map(str) + '_' + blacklist['start'].map(str) + '_' + blacklist['end'].map(
    str)

repeats = pd.read_table(repeat_file, usecols=[0, 1, 2, 3], names=['chr', 'start', 'end', 'repeat_type'])
keep_list = ['telomere']
repeats = repeats[repeats['repeat_type'].isin(keep_list)]

def peak_blacklist_filter(sample_base):
    #  This function removes the peaks that fall within the ENCODE blacklist of regions of the genome and peaks that are
    #  in the alternative contigs of the hg38 genome complication
    data_file = '{}/peaks/{}_peaks.xls'.format(output_directory, sample_base)

    peak_file = pd.read_table(data_file, comment='#')
    print('total number of peaks', peak_file.shape)
    quit()
    keep_list = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12',
                 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX',
                 'chrY']

    peak_file = peak_file[peak_file['chr'].isin(keep_list)]
    print('refined number of peaks', peak_file.shape)

    peak_file['index_name'] = peak_file['chr'].map(str) + '_' + peak_file['start'].map(str) + '_' + peak_file[
        'end'].map(str)
    overlaps = []
    for i in peak_file['chr'].unique():
        bl_index = blacklist[blacklist['chr'] == i]
        if len(bl_index) == 0:
            continue
        bl_ranges = list()
        for ind, row in bl_index.iterrows():
            rng = set(range(row['start'], row['end']))
            bl_ranges.append(rng)
        peak_index = peak_file[peak_file['chr'] == i]
        for ind, row in peak_index.iterrows():
            rng = set(range(row['start'], row['end']))
            for bl_r in bl_ranges:
                if len(bl_r.intersection(rng)) > 0:
                    # print('Intersection', row['index_name'])
                    overlaps.append(row['index_name'])
    new_data = peak_file[~peak_file['index_name'].isin(overlaps)]

    new_data.to_csv('{}/peaks/{}_blacklist_filtered_peaks.csv'.format(output_directory, sample_base), index=False)
    print('number of non-blacklist peaks', new_data.shape)


def telomere_exclusion(sample_base):
    peak_file = '{}/peaks/{}_blacklist_filtered_peaks.csv'.format(output_directory, sample_base)

    peaks = pd.read_csv(peak_file, comment='#')
    print('total numbered of blacklist filtered hits', peaks.shape)
    peaks['index_name'] = peaks['chr'].map(str) + '_' + peaks['start'].map(str) + '_' + peaks['end'].map(str)

    overlaps = []
    for i in peaks['chr'].unique():
        repeat_index = repeats[repeats['chr'] == i]
        if len(repeat_index) == 0:
            continue
        repeat_ranges = list()
        for ind, row in repeat_index.iterrows():
            rng = set(range(row['start'], row['end']))
            repeat_ranges.append(rng)
        peak_index = peaks[peaks['chr'] == i]
        for ind, row in peak_index.iterrows():
            rng = set(range(row['start'], row['end']))
            for repeat_r in repeat_ranges:
                if len(repeat_r.intersection(rng)) > 0:
                    # print('Intersection', row['index_name'])
                    overlaps.append(row['index_name'])
    new_data = peaks[~peaks['index_name'].isin(overlaps)]

    new_data.to_csv('{}/peaks/{}_telomere_filtered.csv'.format(output_directory, sample_base), index=False)
    print('number of non-telomere peaks', new_data.shape)


def peak_control_filter(treatment_file_name, control_file_name, save_name):
    #  This function removes the peaks that are also present in the control peak file

    treatment_file = '{}/peaks/{}.csv'.format(output_directory, treatment_file_name)
    control_file = '{}/peaks/{}.csv'.format(output_directory, control_file_name)

    ctrl_file = pd.read_csv(control_file, comment='#')
    print('total number of control peaks', ctrl_file.shape)

    treat_file = pd.read_csv(treatment_file, comment='#')
    print('total number of treatment peaks', treat_file.shape)

    treat_file['index_name'] = treat_file['chr'].map(str) + '_' + treat_file['start'].map(str) + '_' + treat_file[
        'end'].map(str)
    ctrl_file['index_name'] = ctrl_file['chr'].map(str) + '_' + ctrl_file['start'].map(str) + '_' + ctrl_file[
        'end'].map(str)

    overlaps = []
    for i in treat_file['chr'].unique():
        ctrl_index = ctrl_file[ctrl_file['chr'] == i]
        if len(ctrl_index) == 0:
            continue
        ctrl_ranges = list()
        for ind, row in ctrl_index.iterrows():
            rng = set(range(row['start'], row['end']))
            ctrl_ranges.append(rng)
        treat_index = treat_file[treat_file['chr'] == i]
        for ind, row in treat_index.iterrows():
            rng = set(range(row['start'], row['end']))
            for ctrl_r in ctrl_ranges:
                if len(ctrl_r.intersection(rng)) > 0:
                    # print('Intersection', row['index_name'])
                    overlaps.append(row['index_name'])
    new_data = treat_file[~treat_file['index_name'].isin(overlaps)]
    over = treat_file[treat_file['index_name'].isin(overlaps)]

    new_data.to_csv('{}/peaks/{}_control_filtered.csv'.format(output_directory, save_name), index=False)
    over.to_csv('{}/peaks/{}_overlapping_with_control_peaks.csv'.format(output_directory, save_name), index=False)
    print('number of peaks not in control', new_data.shape)


def peak_filter_blacklist_telomere(sample_base):
    peak_blacklist_filter(sample_base)
    telomere_exclusion(sample_base)

