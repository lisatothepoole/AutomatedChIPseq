import pandas as pd

output_directory = '/Users/lisapoole/Desktop/E76b_asisi_rpa_chip-seq'
ip_file = '/Users/lisapoole/Desktop/E76b_asisi_rpa_chip-seq/'
species = 'human'

# Step 1: import blacklist files
if species == 'human':
    blacklist_file = '/Users/lisapoole/Sources/hg38.blacklist.bed'
    asisi_cut_sites = '/Users/lisapoole/Sources/asisi_cut_sites.csv'
    repeat_file = '/Users/lisapoole/Sources/hg38_repeat_telomere_centromere.txt'
else:
    print('Species value invalid.')

blacklist = pd.read_table(blacklist_file, usecols=[0, 1, 2], names=['chr', 'start', 'end'])
blacklist['index_name'] = blacklist['chr'].map(str) + '_' + blacklist['start'].map(str) + '_' + blacklist['end'].map(
   str)

known = pd.read_csv(asisi_cut_sites)
known['index_name'] = known['chr'].map(str) + '_' + known['start'].map(str) + '_' + known['end'].map(str)

repeats = pd.read_table(repeat_file, usecols=[0,1,2,3], names=['chr', 'start', 'end', 'repeat_type'])
keep_list = ['telomere']
repeats = repeats[repeats['repeat_type'].isin(keep_list)]

def peak_blacklist_filter(sample_base, save_name):
    #  This function removes the peaks that fall within the blacklist of regions of the genome, compiled by ENCODE
    data_file = '{}/peaks/{}_peaks.xls'.format(output_directory, sample_base)

    peak_file = pd.read_table(data_file, comment='#')
    print('total number of peaks', peak_file.shape)
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

    new_data.to_csv('{}/peaks/{}.csv'.format(output_directory, save_name), index=False)
    print('number of non-blacklist peaks', new_data.shape)

def peak_control_filter(treatment_base, control_base, save_name):
    #  This function removes the peaks that are also present in the etoh control peak file

    treatment_file = '{}/peaks/{}_filtered_peaks.csv'.format(output_directory, treatment_base)
    control_file = '{}/peaks/{}_filtered_peaks.csv'.format(output_directory, control_base)

    ctrl_file = pd.read_csv(control_file, comment='#')
    print('total number of control peaks', ctrl_file.shape)
    #print(ctrl_file.head(10))

    treat_file = pd.read_csv(treatment_file, comment='#')
    print('total number of treatment peaks', treat_file.shape)
    #print(treat_file.head(10))

    treat_file['index_name'] = treat_file['chr'].map(str) + '_' + treat_file['start'].map(str) + '_' + treat_file[
        'end'].map(str)
    ctrl_file['index_name'] = ctrl_file['chr'].map(str) + '_' + ctrl_file['start'].map(str) + '_' + ctrl_file[
        'end'].map(str)
    print(treat_file.head(10))
    print(ctrl_file.head(10))

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

    new_data.to_csv('{}/peaks/{}.csv'.format(output_directory, save_name), index=False)
    print('number of peaks not in control', new_data.shape)

def overlap_asisi_known(sample_base, save_name):
    #  This function will find peaks that overlap with known asisi cut sites
    peak_file = '{}/peaks/{}_peaks.xls'.format(output_directory, sample_base)

    peaks = pd.read_csv(peak_file, comment='#')
    print('total number of filtered peaks', peaks.shape)

    peaks['index_name'] = peaks['chr'].map(str) + '_' + peaks['start'].map(str) + '_' + peaks['end'].map(str)
    overlaps = []
    for i in peaks['chr'].unique():
        known_index = known[known['chr'] == i]
        if len(known_index) == 0:
            continue
        known_ranges = list()
        for ind, row in known_index.iterrows():
            rng = set(range(row['start'], row['end']))
            known_ranges.append(rng)
        peak_index = peaks[peaks['chr'] == i]
        for ind, row in peak_index.iterrows():
            rng = set(range(row['start'], row['end']))
            for known_r in known_ranges:
                if len(known_r.intersection(rng)) > 0:
                    # print('Intersection', row['index_name'])
                    overlaps.append(row['index_name'])
    new_data = peaks[peaks['index_name'].isin(overlaps)]

    new_data.to_csv('{}/peaks/{}.csv'.format(output_directory, save_name), index=False)
    print('number of peaks at known sites', new_data.shape)


def telomere_exclusion(sample_base):
    peak_file = '{}/peaks/{}_filtered_peaks.csv'.format(output_directory, sample_base)

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

#peak_control_filter('asisi_4oht', 'asisi_etoh', 'asisi_control_filtered_peaks')
overlap_asisi_known('asisi_4oht', 'asisi_rpa_ip')
#telomere_exclusion('asisi_4oht')