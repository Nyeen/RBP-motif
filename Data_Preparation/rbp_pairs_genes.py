import pandas as pd
import os
import glob
from tqdm import tqdm





genes = pd.read_csv('rip_seq/genes.csv')

files = glob.glob('data/motifs/*')
motifs = {}
for file in files:
    motifs[file.split('\\')[-1].split('_')[0]] = file




hnRNPs = ['HNRNPA1', 'HNRNPA2B1', 'HNRNPC', 'HNRNPD', 'HNRNPF', 'HNRNPH', 'HNRNPH1',
          'HNRNPK', 'HNRNPL', 'HNRNPM', 'HNRNPU', 'HNRNPUL1']
splicing = ['HNRNPA1', 'HNRNPC', 'HNRNPK', 'HNRNPM', 'HNRNPU', 'PCBP1', 'BCLAF1']
silencing = ['DICER1', 'AGO2']
polyadenylation = ['PABPC4', 'PABPN1', 'LARP4']
granules = ['TIA1', 'G3BP1', 'CAPRIN1', 'EIF3A', 'EIF3B', 'EIF3D', 'EIF3G', 'EIF3H', 'DDX6']



# rbps = polyadenylation
rbps = list(set(granules+hnRNPs+polyadenylation+silencing+splicing))
data_save_path = 'data_genes_all/'
if not os.path.isdir(data_save_path):
    os.mkdir(data_save_path)

pair_files = glob.glob(data_save_path+'*')


for i in range(len(rbps)-1):
    for j in range(i+1, len(rbps)):
        filename = rbps[i]+'_'+rbps[j]+'.csv'
        filename_alt = rbps[j]+'_'+rbps[i]+'.csv'
        try:
            pd.read_csv(data_save_path+filename)
            continue
        except:
            pass
        
        try:
            pd.read_csv(data_save_path+filename_alt)
            continue
        except:
            pass
        
        motif1_df = pd.read_csv(motifs[rbps[i]])
        motif2_df = pd.read_csv(motifs[rbps[j]])
        motif1_df['length'] = motif1_df.bs_end - motif1_df.bs_start
        motif2_df['length'] = motif2_df.bs_end - motif2_df.bs_start
        
        print(rbps[i] + '_' + rbps[j])
        chromosomes_considered = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chrX',
                                  'chr8', 'chr9', 'chr11', 'chr10', 'chr12', 'chr13', 'chr14',
                                  'chr15', 'chr16', 'chr17', 'chr18', 'chr20', 'chr19', 'chrY',
                                  'chr22', 'chr21']
        
        df_cols = ['chromosome', 'strand', 'rbp1', 'rbp1_start', 'rbp1_end', 'rbp2', 'rbp2_start', 'rbp2_end']
        pair_df = pd.DataFrame(columns = df_cols)
        df_ix = 0
        
        '''separate chromosome, strand, then sort based on bs_start'''
        for chromosome in chromosomes_considered:
            temp1_df = motif1_df[motif1_df['chromosome'] == chromosome]
            temp2_df = motif2_df[motif2_df['chromosome'] == chromosome]
            tempgenes_df = genes[genes['chromosome'] == chromosome]
            
            temp1_dfp = temp1_df[temp1_df['strandness'] == '+']
            temp2_dfp = temp2_df[temp2_df['strandness'] == '+']
            tempgenes_dfp = tempgenes_df[tempgenes_df['strand'] == '+']
            
            tempgenes_dfp.reset_index(inplace = True)
            
            for m in tqdm (range (len(tempgenes_dfp)), desc= chromosome + ' +...'):    
                temp1_dfp_ = temp1_dfp[tempgenes_dfp.loc[m]['position_start'] - temp1_dfp['bs_start']  <= 1500]
                temp1_dfp_ = temp1_dfp_[tempgenes_dfp.loc[m]['position_end'] - temp1_dfp_['bs_start']  >= 0]
                temp2_dfp_ = temp2_dfp[tempgenes_dfp.loc[m]['position_start'] - temp2_dfp['bs_start']  <= 1500]
                temp2_dfp_ = temp2_dfp_[tempgenes_dfp.loc[m]['position_end'] - temp2_dfp_['bs_start']  >= 0]
                    
                if len(temp1_dfp_) > 0 and len(temp2_dfp_) > 0:
                    for _, row1 in temp1_dfp_.iterrows():
                        for _, row2 in temp2_dfp_.iterrows():
                            if row2['bs_start']-row1['bs_start'] < 20 and row2['bs_start']-row1['bs_start'] >= 0:
                                pair_df.loc[df_ix] = [row1['chromosome'], row1['strandness'], row1['rbp'],
                                                      row1['bs_start'], row1['bs_end'], row2['rbp'], row2['bs_start'],
                                                      row2['bs_end']]
                                df_ix += 1
                            
                            elif row1['bs_start']-row2['bs_start'] < 20 and row1['bs_start']-row2['bs_start'] > 0:
                                pair_df.loc[df_ix] = [row1['chromosome'], row1['strandness'], row2['rbp'],
                                                      row2['bs_start'], row2['bs_end'], row1['rbp'], row1['bs_start'],
                                                      row1['bs_end']]
                                df_ix += 1
            
            
            temp1_dfn = temp1_df[temp1_df['strandness'] == '-']
            temp2_dfn = temp2_df[temp2_df['strandness'] == '-']
            tempgenes_dfn = tempgenes_df[tempgenes_df['strand'] == '-']
            
            tempgenes_dfn.reset_index(inplace = True)
            
            for m in tqdm (range (len(tempgenes_dfn)), desc= chromosome + ' -...'):    
                temp1_dfn_ = temp1_dfn[tempgenes_dfn.loc[m]['position_start'] - temp1_dfn['bs_start']  <= 1500]
                temp1_dfn_ = temp1_dfn_[tempgenes_dfn.loc[m]['position_end'] - temp1_dfn_['bs_start']  >= 0]
                temp2_dfn_ = temp2_dfn[tempgenes_dfn.loc[m]['position_start'] - temp2_dfn['bs_start']  <= 1500]
                temp2_dfn_ = temp2_dfn_[tempgenes_dfn.loc[m]['position_end'] - temp2_dfn_['bs_start']  >= 0]
                    
                if len(temp1_dfn_) > 0 and len(temp2_dfn_) > 0:
                    for _, row1 in temp1_dfn_.iterrows():
                        for _, row2 in temp2_dfn_.iterrows():
                            if row2['bs_start']-row1['bs_start'] < 20 and row2['bs_start']-row1['bs_start'] >= 0:
                                pair_df.loc[df_ix] = [row1['chromosome'], row1['strandness'], row1['rbp'],
                                                      row1['bs_start'], row1['bs_end'], row2['rbp'], row2['bs_start'],
                                                      row2['bs_end']]
                                df_ix += 1
                            
                            elif row1['bs_start']-row2['bs_start'] < 20 and row1['bs_start']-row2['bs_start'] > 0:
                                pair_df.loc[df_ix] = [row1['chromosome'], row1['strandness'], row2['rbp'],
                                                      row2['bs_start'], row2['bs_end'], row1['rbp'], row1['bs_start'],
                                                      row1['bs_end']]
                                df_ix += 1
                    
        
        pair_df.to_csv(data_save_path + rbps[i] + '_' + rbps[j] + '.csv', index = False)

