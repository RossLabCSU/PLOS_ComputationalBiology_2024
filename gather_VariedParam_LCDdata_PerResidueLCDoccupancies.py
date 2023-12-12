
from Bio import SeqIO
import os
from tqdm import tqdm
cwd = os.getcwd()
domains = ['Archaea', 'Bacteria', 'Eukaryota', 'Viruses']

def main():

    cwd = os.getcwd()
    proteome_size_df = get_proteome_sizes()
    
    output = open('LCDcomposer_VariedParameters_Results_PerResidueLCDoccupancies.tsv', 'w')
    output.write('\t'.join(['Domain', 'Organism ID', 'Window Size', 'Primary Amino Acid', 'Secondary Amino Acid', 'Primary Composition Threshold', 'Secondary Composition Threshold', 'Per-Residue LCD Occupancy']) + '\n')
    
    data_types = ['VariedCompThresholds', 'VariedWindowSize']
    for i, data_type in enumerate(data_types):
        for domain in domains:
            dirname = domain + '_' + data_type
            base_path = os.path.join(cwd, 'RandomlySelectedOrganisms', dirname)
        
            files = os.listdir(os.path.join(base_path))
            files = [x for x in files if x.endswith('RESULTS.tsv')]
            for file in tqdm(files):
                fileinfo = file.split('_')
                proteome = fileinfo[0] + '_' + fileinfo[1]
                lcd_class = fileinfo[3]
                if data_type == 'VariedCompThresholds':
                    window_size = '20'
                    if len(lcd_class) == 1:
                        prim_aa = lcd_class
                        sec_aa = 'none'
                        prim_threshold = fileinfo[4]
                        sec_threshold = 'none'
                    else:
                        prim_aa, sec_aa = lcd_class.split('-')
                        prim_threshold, sec_threshold = fileinfo[4].split('-')
                        if int(sec_threshold) > int(prim_threshold):
                            prim_aa, sec_aa = sec_aa, prim_aa
                            prim_threshold, sec_threshold = sec_threshold, prim_threshold
                else:
                    window_size = fileinfo[4].replace('windowSize', '')
                    if window_size == '20':
                        continue
                    if len(lcd_class) == 1:
                        prim_aa = lcd_class
                        sec_aa = 'none'
                        prim_threshold = '40'
                        sec_threshold = 'none'
                    else:
                        prim_aa, sec_aa = lcd_class.split('-')
                        prim_threshold = '40'
                        sec_threshold = '20'
                per_res_occup = calculate_per_res_occupancy(file, base_path, proteome, proteome_size_df)
                output.write('\t'.join([domain, proteome, window_size, prim_aa, sec_aa, prim_threshold, sec_threshold, str(per_res_occup)]) + '\n')
                if prim_threshold == sec_threshold and prim_aa != sec_aa:
                    output.write('\t'.join([domain, proteome, window_size, sec_aa, prim_aa, sec_threshold, prim_threshold, str(per_res_occup)]) + '\n')
    output.close()
    

def calculate_per_res_occupancy(file, base_path, proteome, proteome_size_df):

    lcd_aas = 0
    h = open(os.path.join(base_path, file))
    for i in range(11):
        h.readline()
    for line in h:
        items = line.rstrip().split('\t')
        lcd_len = len(items[2])
        lcd_aas += lcd_len
    h.close()
    
    perres_occup = lcd_aas / proteome_size_df[proteome] * 100

    return perres_occup
    
    
def get_proteome_sizes():

    df = {}
    for domain in domains:
        files = os.listdir(domain)
        fasta_files = [file for file in files if file.endswith('.fasta')]
        for file in fasta_files:
            h = open(os.path.join(domain, file))
            total_aas = 0
            for seq_record in SeqIO.parse(h, 'fasta'):
                seq = str(seq_record.seq)
                if seq[-1] == '*':
                    seq = seq[:-1]
                total_aas += len(seq)
                
            h.close()
            
            proteome = file.replace('.fasta', '')
            df[proteome] = total_aas
    
    return df


if __name__ == '__main__':
    main()