import pandas as pd
import matplotlib.pyplot as plt

#  Create a dictionary of locations of ends of chromosomes -- do not need to change unless there is an update to hg38
ends = {'chr1': 248956422, 'chr2': 242193529, 'chr3': 198295559, 'chr4': 190214555, 'chr5': 181538259,
        'chr6': 170805979, 'chr7': 159345973, 'chr8': 145138636, 'chr9': 138394717, 'chr10': 133797422,
        'chr11': 135086622, 'chr12': 133275309, 'chr13': 114364328, 'chr14': 107043718, 'chr15': 101991189,
        'chr16': 90338345, 'chr17': 83257441, 'chr18': 80373285, 'chr19': 58617616, 'chr20': 64444167,
        'chr21': 46709983, 'chr22': 50818468, 'chrX': 156040895, 'chrY': 57227415}

hg38_repeat_regions = '/Users/lisapoole/Sources/hg38_repeat_telomere_centromere.txt'

# Read in peak file
def peak_location(sample_base):
    output_directory = '/Users/lisapoole/Desktop/E71_SMARCAL1_ENDSEQ/default_homer_analysis/{}'.format(sample_base)
    peak_file = '{}/{}_filtered_peaks.csv'.format(output_directory, sample_base)

    peaks = pd.read_csv(peak_file, header=0)
    result = []
    distance = []
    for i,row in peaks.iterrows():
        Chr = row['chr']
        End = row['end']
        actual_end = ends[Chr]
        d_to_end = actual_end-End
        d_to_start = row['start']
        final_d = d_to_start
        direction = 'start_closer'
        if d_to_start > d_to_end:
            direction = 'end_closer'
            final_d = d_to_end
        result.append(direction)
        distance.append(final_d)
        #print("peak = {}, distance = {}, distance2={}, direction={}".format(
        #    row['PeakID'],d_to_end, d_to_start, direction
        #))
    peaks['closest'] = result
    peaks['distance'] = distance
    #print(peaks.head(10))
    peaks.to_csv('{}/{}_closer_peaks.csv'.format(output_directory, sample_base), index=False)
    group = peaks.groupby(by='chr')
    for i,j in group:
        #print(i,j)
        x = len(j[j['closest'] == 'end_closer'])
        y = len(j[j['closest'] == 'start_closer'])
        print(i,x,y)
        plt.hist(x, bins=20, label='end_closer')
        #plt.hist(y, bins = 20)
        plt.show()
