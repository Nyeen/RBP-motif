import pandas as pd
import numpy as np
from random import randint
from tqdm import tqdm


def scan_cluster(df):

    clusters = []
    cluster_lengths = []
    
    while len(df):
        
        sorted_df = df.sort_values(by=['rbp1_start']).reset_index(drop=True)
        cluster_RBPs = [sorted_df.loc[0].rbp1]
        cluster_span = [sorted_df.loc[0].rbp1_start, sorted_df.loc[0].rbp1_start]
        
        while True:
            tmp_df = sorted_df[sorted_df['rbp1'] == cluster_RBPs[-1]]
            tmp_df = tmp_df[tmp_df['rbp1_start'] == cluster_span[1]]
            if len(tmp_df) and max(tmp_df['score']) > 0.5:
                selected_row = tmp_df[tmp_df['score'] == max(tmp_df['score'])].iloc[0]
                index = tmp_df[tmp_df['score'] == max(tmp_df['score'])].index[0]
                sorted_df.drop(index, inplace = True)
                cluster_RBPs.append(selected_row.rbp2)
                cluster_span[1] = selected_row.rbp2_start
                span_end = selected_row.rbp2_end
            else:
                break
            
        df = df[df['rbp1_start'] > cluster_span[1]]
        
        if len(cluster_RBPs)>2:
            clusters.append(cluster_RBPs)
            cluster_span[1] = span_end
            cluster_lengths.append(cluster_span)
            
        
    
    return clusters, cluster_lengths



def get_crm_indexes(rbps, known_cluster, used_i = [], unique_cluster = False):
    cluster_i = []
    found_rbps = []
    for i in range(len(rbps)):
        if i in used_i and unique_cluster:
            continue
        count = 0
        for rbp in rbps[i]:
            if rbp in known_cluster:
                count += 1
        if count>1:
            cluster_i.append(i)
            if unique_cluster:
                used_i.append(i)
            for rbp in rbps[i]:
                if rbp in known_cluster and rbp not in found_rbps:
                    found_rbps.append(rbp)
        
        if set(found_rbps) == set(known_cluster):
            break
    
    if unique_cluster:
        return cluster_i, known_cluster, used_i
    else:
        return cluster_i, found_rbps


######################################################################################################
######################################################################################################

gene_data_df = pd.read_csv('Data_Preparation/genes.csv')









def get_rightCRC(motif, rbp, df1):
    df1 = df1.reset_index(drop = True)
    df1 = df1[df1['rbp1_start'] == motif[0]]
    df1 = df1[df1['rbp1_end'] == motif[1]]
    df1 = df1[df1['rbp1'] == rbp]
    # df1 = df1[df1['bs1'] == motif]
    if not df1.shape[0]:
        CRC = []
        rbps = []
        score = []
        return CRC, rbps, score
    scanned = []
    CRC = []
    rbps = []
    score = []
    for ix, row in df1.iterrows():
        motif1 = tuple(row[['rbp1_start','rbp1_end']])
        rbp1 = row.rbp1
        motif2 = tuple(row[['rbp2_start','rbp2_end']])
        rbp2 = row.rbp2

        if motif1 not in scanned:
            
            df2 = df1[df1['rbp1_start'] == motif1[0]]
            df2 = df2[df2['rbp1_end'] == motif1[1]]
            df2 = df2[df2['rbp1'] == rbp1]
            
            max_score = max(df2['score'])
            selected_motif = tuple(df2[df2['score'] == max_score].iloc[0][['rbp2_start','rbp2_end']])
            selected_rbp = df2[df2['score'] == max_score].iloc[0]['rbp2']
    #         print(df2)
            CRC.append(selected_motif)
            rbps.append(selected_rbp)
            score.append(max_score)
            next_search_index = max(list(df2.index)) + 1
            if df1[next_search_index:].shape[0]:
                new_CRC, new_rbps, new_score = get_rightCRC(selected_motif, selected_rbp, df1[next_search_index:])
                CRC.extend(new_CRC)
                rbps.extend(new_rbps)
                score.extend(new_score)
            else:
                CRC.extend([])
                rbps.extend([])
                score.extend([])
            scanned.append(motif1)
    
        return CRC, rbps, score
    
def get_leftCRC(motif1, rbp1, df2, df1):
    CRC = []
    rbps = []
    score = []
    max_score = max(df2['score'])
    selected_motif = tuple(df2[df2['score'] == max_score].iloc[0][['rbp1_start','rbp1_end']])
    selected_rbp = df2[df2['score'] == max_score].iloc[0]['rbp1']
    
    CRC.append(selected_motif)
    rbps.append(selected_rbp)
    score.append(max_score)
    
    
    df2 = df1[df1['rbp2_start'] == selected_motif[0]]
    df2 = df2[df2['rbp2_end'] == selected_motif[1]]
    df2 = df2[df2['rbp2'] == selected_rbp]
#     print(selected_motif)
#     print(df2.shape[0])
    if df2.shape[0]:
        new_CRC, new_rbps, new_score = get_leftCRC(selected_motif, selected_rbp, df2.copy(), df1.copy())
        new_CRC.extend(CRC)
        new_rbps.extend(rbps)
        new_score.extend(score)
        CRC = new_CRC.copy()
        rbps = new_rbps.copy()
        score = new_score.copy()
    else:
        return CRC, rbps, score
    
    return CRC, rbps, score
    

def get_CRCs(df1):
    
#     scanned = []
    all_CRCs = []
    all_rbps = []
    all_scores = []
    
    for ix, row in df1.iterrows():
        # print(ix)
        motif1 = tuple(row[['rbp1_start','rbp1_end']])
        rbp1 = row.rbp1
        motif2 = tuple(row[['rbp2_start','rbp2_end']])
        rbp2 = row.rbp2
        co_score = row['score']

#         if motif1 not in scanned:
#         print(motif1)
        CRC = []
        rbps = []
        score = []
        df2 = df1[df1['rbp1_start'] == motif1[0]]
        df2 = df2[df2['rbp1_end'] == motif1[1]]
        df2 = df2[df2['rbp1'] == rbp1]
        max_score = max(df2['score'])
        selected_motif = tuple(df2[df2['score'] == max_score].iloc[0][['rbp2_start','rbp2_end']])
        selected_rbp = df2[df2['score'] == max_score].iloc[0]['rbp2']
#         print(df2)
        CRC.append(motif1)
        rbps.append(rbp1)
        CRC.append(motif2)
        rbps.append(rbp2)
        score.append(co_score)
        next_search_index = max(list(df2.index)) + 1
        # df2 = df1[df1['bs2'] == motif1]
        df2 = df1[df1['rbp2_start'] == row.rbp1_start]
        df2 = df2[df2['rbp2_end'] == row.rbp1_end]
        if df2.shape[0]:
#             print(df1)
            new_CRC, new_rbps, new_score = get_leftCRC(motif1, rbp1, df2.copy(), df1.copy())
            new_CRC.extend(CRC)
            new_rbps.extend(rbps)
            new_score.extend(score)
            CRC = new_CRC.copy()
            rbps = new_rbps.copy()
            score = new_score.copy()

        if df1[next_search_index:].shape[0]:
#             print(motif1)
#             print(df1[next_search_index:])
            new_CRC, new_rbps, new_score = get_rightCRC(selected_motif, selected_rbp, df1[next_search_index:].copy())
            CRC.extend(new_CRC)
            rbps.extend(new_rbps)
            score.extend(new_score)
        else:
            CRC.extend([])
            rbps.extend([])
            score.extend([])
#         scanned.append(motif1)
        all_CRCs.append(CRC)
        all_rbps.append(rbps)
        all_scores.append(round(np.mean(score),3))
    return all_CRCs, all_rbps, all_scores

def remove_duplicate_etc(CRCs, rbps, scores, first = 1):
    new_CRCs = []
    new_rbps = []
    new_scores = []
    for i in range(len(CRCs)):
        cs=0
        ce=0
        for j in range(len(CRCs[i])):
            if cs==0 or cs>CRCs[i][j][0]:
                cs=CRCs[i][j][0]
            if CRCs[i][j][1]>ce:
                ce = CRCs[i][j][1]
        
        cluster = str(cs)+'_'+str(ce)
        proteins = '_'.join(rbps[i])
        if cluster not in new_CRCs and proteins not in new_rbps:
            new_CRCs.append(cluster)
            new_rbps.append(proteins)
            new_scores.append(scores[i])
    
    return new_CRCs, new_rbps, new_scores


def generate_crc(predictions_df, output_filepath_bed, output_filepath_csv):

    putative_CRC = pd.DataFrame(columns = ['chromosome', 'gene_name', 'cluster_region', 'cluster_rbps', 'CRC_score'])

    with open(output_filepath_bed, 'w') as fp:
        
        ngs = []
        chrs = []
        gene_CRCs = []
        gene_rbps = []
        probs = []
        
        for i in tqdm (range (len(gene_data_df)), desc= 'Genes...'):
            # print(i)
            temp_df = predictions_df[predictions_df['chromosome'] == gene_data_df.loc[i]['chromosome']]
            temp_df = temp_df[temp_df['strand'] == gene_data_df.loc[i]['strand']]
            temp_df = temp_df[gene_data_df.loc[i]['position_start'] - temp_df['rbp1_start']  <= 1500]
            temp_df = temp_df[gene_data_df.loc[i]['position_end'] - temp_df['rbp1_start']  >= 0]
            
            gene = gene_data_df.loc[i]['gene_name']
            chromosome = gene_data_df.loc[i]['chromosome']
        
            df1 = temp_df.reset_index(drop = True)
            if not len(df1):
                continue

            CRCs, rbps, scores = get_CRCs(df1.copy())
            CRCs, rbps, scores = remove_duplicate_etc(CRCs, rbps, scores)

            ngs.extend([gene] * len(CRCs))
            chrs.extend([chromosome] * len(CRCs))
            gene_CRCs.extend(CRCs)
            gene_rbps.extend(rbps)
            probs.extend(scores)
            
        
        putative_CRC['chromosome'] = chrs
        putative_CRC['gene_name'] = ngs
        putative_CRC['cluster_region'] = gene_CRCs
        putative_CRC['cluster_rbps'] = gene_rbps
        putative_CRC['CRC_score'] = probs

        putative_CRC.to_csv(output_filepath_csv, index = False)


    return chrs, ngs, gene_CRCs, gene_rbps, probs



