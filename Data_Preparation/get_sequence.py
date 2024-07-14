import pandas as pd
from tqdm import tqdm
from Bio import SeqIO

positive = pd.read_csv('positive.csv')
negative = pd.read_csv('negative.csv')

positive = positive.sort_values(by=['chromosome']).reset_index(drop=True)
negative = negative.sort_values(by=['chromosome']).reset_index(drop=True)

fasta_file = 'rip_seq/hg38.fa'

chrom = 'chr'
positive_bs1 = []
positive_bs2 = []
for i in tqdm (range(len(positive)), desc="Positive..."):
    
    if positive.loc[i]['chromosome'] != chrom:
        fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')
        for fasta in fasta_sequences:
            if fasta.id == positive.loc[i]['chromosome']:
                chrom = fasta.id
                break

    bs1 = str(fasta.seq).upper()[positive.loc[i]['rbp1_start']:positive.loc[i]['rbp1_end']]
    bs2 = str(fasta.seq).upper()[positive.loc[i]['rbp2_start']:positive.loc[i]['rbp2_end']]
    
    positive_bs1.append(bs1)
    positive_bs2.append(bs2)

positive['bs1'] = positive_bs1
positive['bs2'] = positive_bs2
positive.to_csv('positive.csv', index = False)





negative_bs1 = []
negative_bs2 = []
for i in tqdm (range(len(negative)), desc="negative..."):
    
    if negative.loc[i]['chromosome'] != chrom:
        fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')
        for fasta in fasta_sequences:
            if fasta.id == negative.loc[i]['chromosome']:
                chrom = fasta.id
                break

    bs1 = str(fasta.seq).upper()[negative.loc[i]['rbp1_start']:negative.loc[i]['rbp1_end']]
    bs2 = str(fasta.seq).upper()[negative.loc[i]['rbp2_start']:negative.loc[i]['rbp2_end']]
    
    negative_bs1.append(bs1)
    negative_bs2.append(bs2)

negative['bs1'] = negative_bs1
negative['bs2'] = negative_bs2
negative.to_csv('negative.csv', index = False)
    
