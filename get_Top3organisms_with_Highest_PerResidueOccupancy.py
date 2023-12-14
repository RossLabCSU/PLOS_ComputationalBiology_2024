
from tqdm import tqdm
import os
from Bio import SeqIO
import statistics

domains = ['Archaea', 'Bacteria', 'Eukaryota', 'Viruses']

amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
aa_strings = [aa for aa in amino_acids]
for res1 in amino_acids:
    for res2 in amino_acids:
        if res1 == res2:
            continue
        aa_strings.append(res1+res2)
        
def main():

    proteome_length_df = {}
    total_prots_df = {}
    total_aas_df = {}
    window_size = 20
    top_x_threshold = 3
    
    perc_in_lcds_df, proteomes_df, org_info = get_data('LCDfrequency_data_PerResidue_LCDoccupancy.tsv')
    mean_percs_df = get_domain_wide_mean_percs('LCDfrequency_data_PerResidue_LCDoccupancy.tsv')
    
    output = open('TableS11_PerResidueLevelStatistics_Top' + str(top_x_threshold) + 'values_for_LCDclasses.tsv', 'w')
    header = ['LCD Class', 'Domain of Life', 'Maximum Value', 'Value Category', 'Proteome', 'Scientific Name', 'Common Name (when available)', 'Total Proteins in Proteome', 'Total Number of AAs in Proteome (only proteins at least ' + str(window_size) + 'aa in length)', 'Domain-Wide Mean Percentage of Proteome in LCDs']
    output.write('\t'.join(header) + '\n')
    for domain in domains:
        
        proteomes = proteomes_df[domain]
        for proteome in tqdm(proteomes):
            proteome_length_df, total_prots_df, total_aas = calc_proteome_stats(domain, proteome+'.fasta', proteome, proteome_length_df, total_prots_df, window_size)
            total_aas_df[proteome] = total_aas
            
        for aa_string in aa_strings:
            perc_in_lcds = perc_in_lcds_df[domain][aa_string]
            
            domain_wide_mean_perc = mean_percs_df[domain][aa_string]
            
            sorted_percs, perc_sorted_proteomes = zip(*sorted(zip(perc_in_lcds, proteomes), reverse=True))
            best_proteomes = perc_sorted_proteomes[:top_x_threshold]
            proteome_indices = [proteomes.index(proteome) for proteome in best_proteomes]
            temp_total_aas = [total_aas_df[proteome] for proteome in best_proteomes]
            if sorted_percs[0] == 0.0:
                output.write('\t'.join([aa_string, domain, '; '.join([str(x) for x in sorted_percs[:top_x_threshold]]), 'Percentage of Proteome in LCDs', '; '.join(['N/A']*3)] + ['; '.join(['N/A']*3), '; '.join(['N/A']*3), '; '.join(['N/A']*3)] + ['; '.join(['N/A']*3)] + [str(domain_wide_mean_perc)]) + '\n')
            else:
                output.write('\t'.join([aa_string, domain, '; '.join([str(x) for x in sorted_percs[:top_x_threshold]]), 'Percentage of Proteome in LCDs', '; '.join(best_proteomes)] + ['; '.join([org_info[proteome][0] for proteome in best_proteomes]), '; '.join([org_info[proteome][1] for proteome in best_proteomes]), '; '.join([org_info[proteome][2] for proteome in best_proteomes])] + ['; '.join([str(x) for x in temp_total_aas])] + [str(domain_wide_mean_perc)]) + '\n')

    output.close()
    
    
def calc_proteome_stats(domain, proteome, proteome_name, proteome_length_df, total_prots_df, window_size):
    
    total_aas = 0
    h = open(os.path.join('.', domain, proteome))
    for seq_record in SeqIO.parse(h, 'fasta'):
        seq = str(seq_record.seq)
        if seq[-1] == '*':
            seq = seq[:-1]
            
        if len(seq) < window_size:
            continue
            
        total_aas += len(seq)
        
        total_prots_df[proteome_name] = total_prots_df.get(proteome_name, 0) + 1
        
    h.close()
    
    proteome_length_df[proteome_name] = total_aas
    
    return proteome_length_df, total_prots_df, total_aas
    
    
def get_data(file):

    df = {domain:{aa_string:[] for aa_string in aa_strings} for domain in domains}
    org_info = {}
    
    h = open(file)
    header = h.readline()
    proteomes = {domain:[] for domain in domains}
    for line in h:
        proteome, taxid, domain, sci_name, common_name, total_prots, *data = line.rstrip().split('\t')
        if proteome == 'UP000597762_158019':    # EXCLUDE SEPIA PHARAONIS
            continue
        data = [float(x) if x != 'N/A' else 0 for x in data]
        for i, aa_string in enumerate(aa_strings):
            val = data[i]
            df[domain][aa_string].append(val)
        proteomes[domain].append(proteome)
        org_info[proteome] = [sci_name, common_name, total_prots]
        
    h.close()
    
    return df, proteomes, org_info
    
    
def get_domain_wide_mean_percs(file):

    h = open(file)
    header = h.readline()
    
    df = {domain:{aa_string:[] for aa_string in aa_strings} for domain in domains}
    
    for line in h:
        proteome, taxid, domain, sci_name, common_name, total_prots, *data = line.rstrip().split('\t')
        data = [float(x) for x in data]
        
        for i, aa_string in enumerate(aa_strings):
            df[domain][aa_string].append(data[i])
    h.close()
            
    means_df = {domain:{} for domain in domains}
    for domain in df:
        for aa_string in aa_strings:
            ave_perc = statistics.mean(df[domain][aa_string])
            means_df[domain][aa_string] = ave_perc
            
    return means_df
    
    
if __name__ == '__main__':
    main()
