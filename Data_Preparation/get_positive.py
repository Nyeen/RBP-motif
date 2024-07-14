import pandas as pd
import glob

pairs = glob.glob('data_genes_all/*')

df_cols = ['chromosome', 'strand', 'rbp1', 'rbp1_start', 'rbp1_end', 'rbp2', 'rbp2_start', 'rbp2_end']
positive = pd.DataFrame(columns = df_cols)
for pair in pairs:
    sub_data = pd.read_csv(pair)
    positive = pd.concat([positive, sub_data], axis=0)

positive['bs1_length']=positive.rbp1_end-positive.rbp1_start
positive['bs2_length']=positive.rbp2_end-positive.rbp2_start
positive = positive[positive.bs1_length<23]
positive = positive[positive.bs2_length<23]
positive = positive[positive.bs1_length>2]
positive = positive[positive.bs2_length>2]
positive['dbbs'] = abs(positive.rbp1_start-positive.rbp2_start)
positive = positive[positive.dbbs<20]

positive = positive.drop_duplicates().reset_index(drop=True)
positive.drop(columns=['bs1_length', 'bs2_length', 'dbbs'], inplace = True)
positive.to_csv('positive.csv', index = False)
