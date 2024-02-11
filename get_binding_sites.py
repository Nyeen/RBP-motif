from Bio import SeqIO
from tqdm import tqdm

fasta_file = 'RIP-seq/hg38.fa'
RBP_human = 'POSTAR3-peaks/human.txt'


fr = open(RBP_human, 'r')
fw = open('RBP_human.txt', 'w')

chrom = 'chr'

for i in tqdm (range (1000000), desc="Converting..."):
    
    line = fr.readline()        
    if not line:
        break

    bs_info = line.split('\t')
    
    if bs_info[0] != chrom:
        fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')
        for fasta in fasta_sequences:
            if fasta.id == bs_info[0]:
                chrom = fasta.id
                break

    sequence = str(fasta.seq).upper()[int(bs_info[1]):int(bs_info[2])]
    fw.write(bs_info[0]+'\t'+bs_info[1]+'\t'+bs_info[2]+'\t'+sequence+'\t'+
             bs_info[4]+'\t'+bs_info[5]+'\t'+bs_info[7]+'\n')

fr.close()
fw.close()
