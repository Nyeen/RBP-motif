import pandas as pd
from tqdm import tqdm

motif_data = 'RBP_human.txt'

'''motif data has been loaded in a pandas dataframe with proper column label'''
all_motifs = pd.read_csv(motif_data, sep='\t', names=['chromosome', 'bs_start', 'bs_end',
                                                      'sequence', 'strandness', 'rbp',
                                                      'tissue'])



'''The range of motif length has taken arbitrarily'''
motif_length_lower = 5
motif_length_upper = 20


'''Preparing the dataframe for scanning crm data'''
all_motifs = all_motifs.sort_values(by = 'bs_start')
all_motifs['length'] = all_motifs.bs_end-all_motifs.bs_start
all_motifs.drop(all_motifs[all_motifs.length < motif_length_lower].index, inplace = True)
all_motifs.drop(all_motifs[all_motifs.length > motif_length_upper].index, inplace = True)

'''prepare an output file to store all crm data'''
fw = open('crm_data.txt', 'w')
fw.write('chromosome\tbs1_start\tbs2_start\tbs1_end\tbs2_end\tbs1_sequence\tbs2_sequence\t'+
         'bs1_rbp\tbs2_rbp\tbs1_tissue\tbs2_tissue\tcrm\n')


for i in tqdm (range (len(all_motifs)-1), desc="Iterating..."):
    for j in range(i+1, len(all_motifs)):
        
        # check if the pair have the same motif
        if all_motifs.iloc[i].rbp == all_motifs.iloc[j].rbp:
            continue
        
        distance = all_motifs.iloc[j].bs_start - all_motifs.iloc[i].bs_start
        
        # check if the pairing going outside the promoter region
        if distance > 1500:
            break
        
        fw.write(all_motifs.iloc[i].chromosome+'\t'+
                 str(all_motifs.iloc[i].bs_start)+'\t'+str(all_motifs.iloc[j].bs_start)+'\t'+
                 str(all_motifs.iloc[i].bs_end)+'\t'+str(all_motifs.iloc[j].bs_end)+'\t'+
                 all_motifs.iloc[i].sequence+'\t'+all_motifs.iloc[j].sequence+'\t'+
                 all_motifs.iloc[i].rbp+'\t'+all_motifs.iloc[j].rbp+'\t'+
                 all_motifs.iloc[i].tissue+'\t'+all_motifs.iloc[j].tissue+'\t')
        
        
        if distance > 30:
            fw.write('0\n')
        else:
            fw.write('1\n')

fw.close()
