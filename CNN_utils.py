import pandas as pd
import numpy as np
from tqdm import tqdm
from Bio import SeqIO
import random


def add_flanking_sequence(df, length = 30):
    fasta_file = 'rip_seq/hg38.fa'
    chrom = 'chr'
    bs1_list = []
    bs2_list = []
    for i in tqdm (range(len(df)), desc="Adding flanking sequence..."):
        
        if df.loc[i]['chromosome'] != chrom:
            fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')
            for fasta in fasta_sequences:
                if fasta.id == df.loc[i]['chromosome']:
                    chrom = fasta.id
                    break
        
        bs1 = str(fasta.seq).upper()[df.loc[i]['rbp1_start']-length:df.loc[i]['rbp1_end']]
        bs2 = str(fasta.seq).upper()[df.loc[i]['rbp2_start']:df.loc[i]['rbp2_end']+length]
        
        bs1_list.append(bs1)
        bs2_list.append(bs2)
    
    df['bs1'] = bs1_list
    df['bs2'] = bs2_list
    
    return df


def data_encoder(bs1, bs2, mode):
    nucleotide = {'A': np.array([0,0,0,1]),
                  'C': np.array([0,0,1,0]),
                  'G': np.array([0,1,0,0]),
                  'T': np.array([1,0,0,0]),
                  'X': np.array([0,0,0,0])}
    
    bs = []
    for l in bs1:
        bs.append(nucleotide[l])
    bs1 = np.array(bs)
    
    bs = []
    for l in bs2:
        bs.append(nucleotide[l])
    bs2 = np.array(bs)
    
    if mode == 'group':
        return np.array([bs1, bs2])
    
    elif mode == 'pstack':
        bs1 = bs1.transpose()
        bs2 = bs2.transpose()
        return np.concatenate((bs1, bs2), axis=0)
    
    elif mode == 'sstack':
        return np.concatenate((bs1, bs2), axis=0).transpose()


def get_data_label(df, label, cut_off, mode = 'group'):
    data = []
    for index, row in df.iterrows():
        bs1 = row.bs1+'X'*(cut_off-len(row.bs1))
        bs2 = row.bs2+'X'*(cut_off-len(row.bs2))
        data.append(data_encoder(bs1, bs2, mode = mode))
    x = np.array(data)
    y = np.array([label]*len(x))
    return x, y

def prepare_train_test(positive, negative, split = 0.2):
    train_positive = positive[:int((1-split)*len(positive))]
    train_negative = negative[:int((1-split)*len(negative))]
    train = train_positive+train_negative
    random.shuffle(train)
    x_train, y_train = zip(*train)
    x_train = np.array(x_train)
    y_train = np.array(y_train)

    test_positive = positive[int((1-split)*len(positive)):]
    test_negative = negative[int((1-split)*len(negative)):]
    test = test_positive+test_negative
    random.shuffle(test)
    x_test, y_test = zip(*test)
    x_test = np.array(x_test)
    y_test = np.array(y_test)
    
    return x_train, y_train, x_test, y_test
