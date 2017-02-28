import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
font = {'fontname':'Times New Roman'}

experiment_name = 'E76_RPA_ChIP-seq'
output_directory = '/Users/lisapoole/Desktop/{}'.format(experiment_name)
bin_range = range(0, 101, 10)
bin_cutoff = 100 # If no cutoff is desired, indicate '0'
y_scale = 'log' # Valid options are 'log' or 'normal'
y_log_scale = '10'  # Commong options include 10 or 2

def peak_pileup_histogram(peak_file, title_name, save_name):
#  This function creates a histogram of reads piled at peak summit
    filename = '{}/peaks/{}.csv'.format(output_directory, peak_file)
    data = pd.read_csv(filename, comment='#', header=0)
    bins = '{}'.format(bin_range)
    data_for_hist = data.copy()
    if bin_cutoff == '0':
        cutoff = '{}'.format(bin_cutoff)
        data_for_hist.loc[data_for_hist['pileup'] > cutoff, 'pileup'] = cutoff + 1

        bins.append(cutoff + 10)
        bin_labels = []
        for i in bins:
            if i == 100:
                bin_labels.append('>100')
            elif i > 100:
                continue
            else:
                bin_labels.append('{}-{}'.format(i, i + 10))
        plt.xlim(0, {}).format(bin_cutoff)
    else:
        bin_labels = bins
        plt.xlim(0, 115)

    fig, ax = plt.subplots(figsize=(9, 5))
    plt.hist(data_for_hist['pileup'], bins=bins,
                       color='red')
    bin_labels = np.array(bins, dtype='|S4')
    bins = np.array(bins)+5
    ax.set_xticks(bins)
    ax.set_xticklabels(bin_labels, fontsize=10)
    if y_scale == 'log':
        plt.yscale('log', base={}.__format__(y_log_scale))
        plt.ylabel('$log_{10}$(Frequency)', **font)
    elif y_scale == 'normal':
        plt.ylabel('Frequency', **font)
    else:
        print('Invalid y_scale choice.')

    plt.title('{}'.format(title_name), fontsize=20, **font)
    plt.xlabel('Read Pileup', **font)
    for tick in ax.get_xticklabels():
        tick.set_fontname("Times New Roman")
    for tick in ax.get_yticklabels():
        tick.set_fontname("Times New Roman")
    plt.savefig("/Users/lisapoole/Desktop/{}.png".format(save_name), dpi=300, bbox_inches='tight')


