
import os
from Bio import SeqIO
from tqdm import tqdm

domains = ['Archaea', 'Bacteria', 'Eukaryota', 'Viruses']
amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
win_size = 20

def main():
    
    # CALCULATES THE AMINO ACID FREQUENCIES IN ALL PROTEOMES, WITH "_additional" ISOFORM FILES EXCLUDED, ONLY CANONICAL AMINO ACIDS COUNTED, AND A MINIMUM WINDOW SIZE OF 20 AAs.
    
    output = open('Background_AAfrequencies_AllProteomes.tsv', 'w')
    output.write('\t'.join(['Domain of Life', 'Proteome'] + list(amino_acids)) + '\n')
    
    output2 = open('Background_AA-RANKS_AllProteomes.tsv', 'w')
    output2.write('\t'.join(['Domain of Life', 'Proteome'] + list(amino_acids)) + '\n')
    for domain in domains:
        for file in tqdm(os.listdir('./' + domain)):
            if not file.endswith('.fasta') or '_DNA' in file or '_additional' in file:
                continue
                
            proteome, ext = file.split('.')
                
            df = {aa:0 for aa in amino_acids}
            
            h = open(os.path.join('.', domain, file))
            total_aas = 0
            for seq_record in SeqIO.parse(h, 'fasta'):
                seq = str(seq_record.seq)
                if len(seq) < win_size:
                    continue
    
                for aa in amino_acids:
                    count = seq.count(aa)
                    df[aa] += count
                    total_aas += count
            h.close()

            output.write('\t'.join([domain, proteome] + [str(df[aa]) for aa in amino_acids]) + '\n')
            
            counts = [df[aa] for aa in amino_acids]
            sorted_counts, sorted_aas = zip(*sorted(zip(counts, list(amino_acids)), reverse=True))
            ranks = [sorted_aas.index(aa) + 1 for aa in amino_acids]
            
            output2.write('\t'.join([domain, proteome] + [str(ranks[i]) for i, aa in enumerate(amino_acids)]) + '\n')
            
    output.close()
    output2.close()

if __name__ == '__main__':
    main()