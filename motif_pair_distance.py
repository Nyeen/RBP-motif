import pandas as pd
import time
import os
import numpy as np

RBP_human = 'POSTAR3-peaks/human.txt'

all_motifs = pd.read_csv(RBP_human, sep = '\t', names = ['chromosome', 'bs_start', 'bs_end',
                                                          'peak_id', 'strandness', 'rbp',
                                                          'experiment', 'tissue', 'access'])


def calc_dist(df1, df2, df_cols):
    
    if not (len(df1) and len(df2)):
        return pd.DataFrame(columns = df_cols)
        
    start1 = np.transpose(np.tile(np.matrix(df1.bs_start), (len(df2), 1)))
    start2 = np.tile(np.matrix(df2.bs_start), (len(df1), 1))
    dis = abs(start1 - start2)
    index1 = np.transpose(np.tile(np.matrix(df1.index), (len(df2), 1)))
    index2 = np.tile(np.matrix(df2.index), (len(df1), 1))
    
    dis = dis.flatten()
    index1 = index1.flatten()
    index2 = index2.flatten()
    
    r,c=np.shape(dis)
    if r*c == 1:
        df=pd.DataFrame(columns = df_cols)
        df.loc[0]=[index1[0,0], index2[0,0], dis[0,0]]
        return df
    
    return pd.DataFrame(np.stack((index1, index2, dis)).T,
                        columns = df_cols)


def make_motif_pair(motif1, motif2):
    
    data_save_path = 'data/'+ motif1 + '_' + motif2 + '/'
    data_save_path_alt = 'data/'+ motif2 + '_' + motif1 + '/'
    if os.path.isdir(data_save_path) or os.path.isdir(data_save_path_alt):
        return
    else:
        os.mkdir(data_save_path)
    motif_save_path = 'data/motifs/'
    
    
    motif1_df = all_motifs[all_motifs.rbp == motif1]
    motif2_df = all_motifs[all_motifs.rbp == motif2]
    
    
    chromosomes_considered = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chrX',
                              'chr8', 'chr9', 'chr11', 'chr10', 'chr12', 'chr13', 'chr14',
                              'chr15', 'chr16', 'chr17', 'chr18', 'chr20', 'chr19', 'chrY',
                              'chr22', 'chr21']
    
    motif1_df = motif1_df[motif1_df['chromosome'].isin(chromosomes_considered)]
    motif2_df = motif2_df[motif2_df['chromosome'].isin(chromosomes_considered)]
    
    motif1_df.drop_duplicates(subset=['chromosome', 'bs_start', 'strandness'], inplace = True)
    motif1_df = motif1_df.reset_index(drop = True)
    motif1_df['index'] = motif1_df.index
    motif1_df.to_csv(motif_save_path + motif1 + '_sequences.csv', index = False)
    
    motif2_df.drop_duplicates(subset=['chromosome', 'bs_start', 'strandness'], inplace = True)
    motif2_df = motif2_df.reset_index(drop = True)
    motif2_df['index'] = motif2_df.index
    motif2_df.to_csv(motif_save_path + motif2 + '_sequences.csv', index = False)
    
    df_cols = [motif1 + '_ix', motif2 + '_ix', 'dbbs']
    _10k_df = pd.DataFrame(columns = df_cols)
    _summary_df = pd.DataFrame(columns=['_10k', '10k_'])
    m = 0
   
    for chromosome in chromosomes_considered:
        print(motif1, motif2, chromosome)
        partition = 5000
        
        temp1_df = motif1_df[motif1_df['chromosome'] == chromosome]
        temp2_df = motif2_df[motif2_df['chromosome'] == chromosome]
        
        temp1_dfp = temp1_df[temp1_df['strandness'] == '+']
        temp2_dfp = temp2_df[temp2_df['strandness'] == '+']
        
        p1 = len(temp1_dfp)//partition + 1
        p2 = len(temp2_dfp)//partition + 1
        for i in range(p1):
            for j in range(p2):
                if p1-i>1:
                    sdf1 = temp1_dfp[i*partition:(i+1)*partition]
                else:
                    sdf1 = temp1_dfp[i*partition:]
                if p2-j>1:
                    sdf2 = temp2_dfp[j*partition:(j+1)*partition]
                else:
                    sdf2 = temp2_dfp[j*partition:]
                
                sdf = calc_dist(sdf1, sdf2, df_cols)
                _10k_temp_df = sdf[sdf['dbbs'] <= 10000]
                _10k_df = pd.concat([_10k_df,_10k_temp_df], axis=0).reset_index(drop=True)
                
                if len(_10k_df)>1000000:
                    _10k_df.to_csv(data_save_path + '_10k_co_occurrence_'+str(m)+'.csv', index = False)
                    _10k_df = pd.DataFrame(columns = df_cols)
                    m += 1
                
                _summary_temp_df = pd.DataFrame(columns=['_10k', '10k_'])
                _summary_temp_df.loc[0] = [len(_10k_temp_df), len(sdf[sdf['dbbs'] > 10000])]
                _summary_df = pd.concat([_summary_df,_summary_temp_df], axis=0).reset_index(drop=True)
                print('P\t',i+1,'/',p1,'\t',j+1,'/',p2)
        
        temp1_dfn = temp1_df[temp1_df['strandness'] == '-']
        temp2_dfn = temp2_df[temp2_df['strandness'] == '-']
        
        p1 = len(temp1_dfn)//partition + 1
        p2 = len(temp2_dfn)//partition + 1
        for i in range(p1):
            for j in range(p2):
                if p1-i>1:
                    sdf1 = temp1_dfn[i*partition:(i+1)*partition]
                else:
                    sdf1 = temp1_dfn[i*partition:]
                if p2-j>1:
                    sdf2 = temp2_dfn[j*partition:(j+1)*partition]
                else:
                    sdf2 = temp2_dfn[j*partition:]
                
                sdf = calc_dist(sdf1, sdf2, df_cols)
                _10k_temp_df = sdf[sdf['dbbs'] <= 10000]
                _10k_df = pd.concat([_10k_df,_10k_temp_df], axis=0).reset_index(drop=True)
                
                if len(_10k_df)>1000000:
                    _10k_df.to_csv(data_save_path + '_10k_co_occurrence_'+str(m)+'.csv', index = False)
                    _10k_df = pd.DataFrame(columns = df_cols)
                    m += 1
                
                _summary_temp_df = pd.DataFrame(columns=['_10k', '10k_'])
                _summary_temp_df.loc[0] = [len(_10k_temp_df), len(sdf[sdf['dbbs'] > 10000])]
                _summary_df = pd.concat([_summary_df,_summary_temp_df], axis=0).reset_index(drop=True)
                print('N\t',i+1,'/',p1,'\t',j+1,'/',p2)
    
    _10k_df.to_csv(data_save_path + '_10k_co_occurrence_'+str(m)+'.csv', index = False)
    _summary_df.to_csv(data_save_path + 'summary_co_occurrence.csv', index = False)
    
    return


hnRNPs = ['HNRNPA1', 'HNRNPA2B1', 'HNRNPC', 'HNRNPD', 'HNRNPF', 'HNRNPH', 'HNRNPH1',
          'HNRNPK', 'HNRNPL', 'HNRNPM', 'HNRNPU', 'HNRNPUL1']
splicing = ['HNRNPA1', 'HNRNPC', 'HNRNPK', 'HNRNPM', 'HNRNPU', 'PCBP1', 'BCLAF1']
silencing = ['DICER1', 'AGO2']
polyadenylation = ['PABPC4', 'PABPN1', 'LARP4']
granules = ['TIA1', 'G3BP1', 'CAPRIN1', 'EIF3A', 'EIF3B', 'EIF3D', 'EIF3G', 'EIF3H', 'DDX6']

rbps = hnRNPs

start = time.time()
for i in range(len(rbps)-1):
    for j in range(i+1, len(rbps)):
        make_motif_pair(rbps[i], rbps[j])
end = time.time()

print("Elapsed time: ", end-start)

