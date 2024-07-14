import pandas as pd
from tqdm import tqdm


chromosomes_considered = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chrX',
                          'chr8', 'chr9', 'chr11', 'chr10', 'chr12', 'chr13', 'chr14',
                          'chr15', 'chr16', 'chr17', 'chr18', 'chr20', 'chr19', 'chrY',
                          'chr22', 'chr21']
    

genes = pd.DataFrame(columns=['chromosome', 'gene_name', 'gene_type', 'position_start', 'position_end', 'gene_source',
                              'strand'])
df_ix = 0


f=open('rip_seq/gencode.v38lift37.annotation.gtf', 'r')

line = f.readline()

while line:
    line_elements = line.split('\t')
    if len(line_elements) == 1:
        pass
    
    elif line_elements[0] in chromosomes_considered and line_elements[2] == 'gene':
        genes.loc[df_ix] = [line_elements[0], line_elements[8].split(';')[2].split('"')[1],
                            line_elements[8].split(';')[1].split('"')[1], int(line_elements[3]),
                            int(line_elements[4]), line_elements[1], line_elements[6]]
        df_ix += 1
        print(df_ix)
    
    line = f.readline()

genes.sort_values(by=['chromosome', 'position_start'])
genes.to_csv('genes.csv', index = False)

f.close()
