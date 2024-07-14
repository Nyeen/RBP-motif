import pandas as pd
from tqdm import tqdm


data = pd.read_csv('positive.csv')


p1 = data[['chromosome', 'strand', 'rbp1', 'rbp1_start', 'rbp1_end']]
p1.rename(columns={"rbp1": "rbp", "rbp1_start": "start", 'rbp1_end': 'end'}, inplace=True)
p2 = data[['chromosome', 'strand', 'rbp2', 'rbp2_start', 'rbp2_end']]
p2.rename(columns={"rbp2": "rbp", "rbp2_start": "start", 'rbp2_end': 'end'}, inplace=True)
p = pd.concat([p1,p2], axis=0)
p=p.drop_duplicates().reset_index(drop=True)


all_motifs = pd.read_csv('postar3-peaks/human.txt', sep = '\t', names = ['chromosome', 'bs_start', 'bs_end',
                                                          'peak_id', 'strandness', 'rbp',
                                                          'experiment', 'tissue', 'access'])


positive_rbps = list(p.rbp.unique())
for rbp in positive_rbps:
    all_motifs = all_motifs.drop(all_motifs[all_motifs.rbp == rbp].index)

all_motifs['length'] = all_motifs.bs_end - all_motifs.bs_start
all_motifs = all_motifs[all_motifs.length<23]
all_motifs = all_motifs[all_motifs.length>2]

df_cols = ['chromosome', 'strand', 'rbp1', 'rbp1_start', 'rbp1_end', 'rbp2', 'rbp2_start', 'rbp2_end']
negative = pd.DataFrame(columns = df_cols)
df_ix = 0
for i in tqdm (range(len(p)), desc="Building negative data..."):
    tmp_motifs = all_motifs[all_motifs.chromosome == p.loc[i].chromosome]
    tmp_motifs = tmp_motifs[tmp_motifs.strandness == p.loc[i].strand]
    tmp_motifs['distance'] = tmp_motifs.bs_start-p.loc[i].start
    tmp_motifs = tmp_motifs[tmp_motifs.distance > 1500].reset_index(drop=True)
    if not len(tmp_motifs):
        continue
    negative.loc[df_ix] = [p.loc[i]['chromosome'], p.loc[i]['strand'], p.loc[i]['rbp'],
                      p.loc[i]['start'], p.loc[i]['end'], tmp_motifs.loc[0]['rbp'], tmp_motifs.loc[0]['bs_start'],
                      tmp_motifs.loc[0]['bs_end']]
    df_ix += 1


negative = negative.drop_duplicates().reset_index(drop=True)
negative = negative.sample(n=len(data))
negative.to_csv('negative.csv', index = False)
