from Bio import SeqIO
import pandas as pd
import glob


def bed_to_csv(bed_file, fasta_file):
    
    with open(bed_file, 'r') as fp:
        all_motifs = fp.readlines()
    
    df = pd.DataFrame(columns = ['chromosome', 'position_start', 'position_end', 'sequence'])
    chromosomes = []
    positions_start = []
    positions_end = []
    sequences = []
    
    for motif in all_motifs:
        info = motif.split('\t')
        chromosome = info[0]
        position_start = int(info[1])
        position_end = int(info[2])
        fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')
        for fasta in fasta_sequences:
            if fasta.id == chromosome:
                sequence = str(fasta.seq).upper()[position_start:position_end]
                break
    
        chromosomes.append(chromosome)
        positions_start.append(position_start)
        positions_end.append(position_end)
        sequences.append(sequence)
    
    df = pd.DataFrame({'chromosome': chromosomes,
                       'position_start': positions_start,
                       'position_end': positions_end,
                       'sequence': sequences})
    
    df.to_csv(bed_file.split('.')[0]+'.csv', index=False)
    
    return

def main():
    fasta_file = 'RIP-seq/hg19.fa'
    bed_files = glob.glob('bed_files/*.bed')
    for i in tqdm (range (len(bed_files)), desc="Converting..."):
    bed_file = bed_files[i]
    bed_to_csv(bed_file, fasta_file)

if __name__ == "__main__":
    main()
    
