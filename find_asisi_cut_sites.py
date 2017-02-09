import pandas as pd

repeat_file = '/Users/lisapoole/Sources/hg38_repeat_telomere_centromere.txt'

repeats = pd.read_table(repeat_file, usecols=[0,1,2,3], names=['chr', 'start', 'end', 'repeat_type'])
keep_list = ['telomere']
repeats = repeats[repeats['repeat_type'].isin(keep_list)]
print(repeats.head(10))