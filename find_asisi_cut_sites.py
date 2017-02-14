import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


output_directory = '/Users/lisapoole/Desktop/E76b_asisi_rpa_chip-seq'
cutoff_of_reads = '30'

def peak_histogram_pileup_plot(sample_base):
#  This function creates a histogram of reads piled at peak summit
    filename = '{}/peaks/{}_telomere_filtered.csv'.format(output_directory, sample_base)
    data = pd.read_csv(filename, comment='#', header=0)

    bins = range(0, 101, 10)
    data_for_hist = data.copy()
    cutoff = 100
    data_for_hist.loc[data_for_hist['pileup'] > cutoff, 'pileup'] = cutoff + 1

    bins.append(cutoff + 10)
    fig, ax = plt.subplots(figsize=(9, 5))
    plt.hist(data_for_hist['pileup'], bins=bins,
                       color='red')
    bin_labels = np.array(bins, dtype='|S4')
    bin_labels = []
    for i in bins:
        if i == 100:
            bin_labels.append('>100')
        elif i >100:
            continue
        else:
            bin_labels.append('{}-{}'.format(i,i+10))
    bins = np.array(bins)+5
    ax.set_xticks(bins)
    ax.set_xticklabels(bin_labels, fontsize=10)
    plt.xlim(0, 115)
    plt.yscale('log', base=10)
    plt.title('{} Tag Counts'.format(sample_base), fontsize=20)
    plt.ylabel('log$10$(Frequency)')
    plt.xlabel('read_pileup')
    plt.savefig("{}/{}_pileup_plot.png".format(output_directory, sample_base), bbox_inches='tight')

def peak_cutoff_filtering(sample_base):
    # This removes peaks that are below the cutoff for number of reads determined to legitimate
    data_file = '{}/peaks/{}_telomere_filtered.csv'.format(output_directory, sample_base)
    peaks = pd.read_csv(data_file, header=0)
    print('number of peaks', peaks.shape)
    new_peaks = peaks[peaks['pileup'] >= 30]
    print('number of peaks above cutoff', new_peaks.shape)
    new_peaks.to_csv('{}/peaks/{}_cutoff_{}.csv'.format(output_directory, sample_base, cutoff_of_reads), index=False)

peak_cutoff_filtering('asisi_4oht')