import pandas as pd
import numpy as np
import CNN_utils as cnn_ut
import CRM_utils as crm_ut
from tqdm import tqdm
from tensorflow.keras.models import load_model

loaded_model = load_model('best_model')



flanking_sequence = 0 # True or False
flanking_sequence_length = 30
cut_off = 22 + flanking_sequence*flanking_sequence_length
mode = 'pstack'

if flanking_sequence:
    try:
        df = pd.read_csv('positive_wf_'+str(flanking_sequence_length)+'.csv')
    except:
        df = pd.read_csv('positive.csv')
        df = cnn_ut.add_flanking_sequence(df, length = flanking_sequence_length)
        df.to_csv('positive_wf_'+str(flanking_sequence_length)+'.csv', index = False)

else:
    df = pd.read_csv('positive.csv')

x_positive, y_positive = cnn_ut.get_data_label(df, label = 1, cut_off = cut_off, mode = mode)
x_positive = x_positive.astype('float32')
if len(np.shape(x_positive)) < 4:
    x_positive = x_positive.reshape((x_positive.shape[0], x_positive.shape[1], x_positive.shape[2], 1))

y_prob = loaded_model.predict(x_positive, verbose=0)
y_prob = y_prob[:, 0]
df['score']=y_prob




'''
Scanning for clusters
'''

chromosome, gene_name, cluster_region, cluster_rbps, CRC_score = crm_ut.generate_crc(df, 'crc.bed', 'crc.csv')




'''
Coverage of pre-defined RBPs of each biological category within identified clusters
'''

hnRNPs = ['HNRNPA1', 'HNRNPA2B1', 'HNRNPC', 'HNRNPD', 'HNRNPF', 'HNRNPH', 'HNRNPH1',
          'HNRNPK', 'HNRNPL', 'HNRNPM', 'HNRNPU', 'HNRNPUL1']
splicing = ['HNRNPA1', 'HNRNPC', 'HNRNPK', 'HNRNPM', 'HNRNPU', 'PCBP1', 'BCLAF1']
silencing = ['DICER1', 'AGO2']
polyadenylation = ['PABPC4', 'PABPN1', 'LARP4']
granules = ['TIA1', 'G3BP1', 'CAPRIN1', 'EIF3A', 'EIF3B', 'EIF3D', 'EIF3G', 'EIF3H', 'DDX6']

rbps = []
no_unique_rbps = []
cluster_size = []
cluster_length = []
for i in range(len(gene_name)):
    if len(np.unique(cluster_rbps[i].split('_'))) > 3:
        rbps.append(list(np.unique(cluster_rbps[i].split('_'))))
        no_unique_rbps.append(len(np.unique(cluster_rbps[i].split('_'))))
        cluster_size.append(len(cluster_rbps[i]))
        cluster_length.append(int(cluster_region[i].split('_')[1])-int(cluster_region[i].split('_')[0]))

        



hnRNPs_i, found_hnRNPs = crm_ut.get_crm_indexes(rbps, hnRNPs)
splicing_i, found_splicing = crm_ut.get_crm_indexes(rbps, splicing)
silencing_i, found_silencing = crm_ut.get_crm_indexes(rbps, silencing)
polyadenylation_i, found_polyadenylation = crm_ut.get_crm_indexes(rbps, polyadenylation)
granules_i, found_granules = crm_ut.get_crm_indexes(rbps, granules)




acc_h = len(found_hnRNPs)/len(hnRNPs)
acc_sp = len(found_splicing)/len(splicing)
acc_si = len(found_silencing)/len(silencing)
acc_p = len(found_polyadenylation)/len(polyadenylation)
acc_g = len(found_granules)/len(granules)




print(acc_h)
print(acc_sp)
print(acc_si)
print(acc_p)
print(acc_g)

