
amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
import os
from Bio import SeqIO
import statistics
from tqdm import tqdm
from scipy import stats
domains = ['Archaea', 'Bacteria', 'Eukaryota', 'Viruses']


def main(args):

    run_scrambled_analysis = args.scrambled_analysis

    aa_strings = [aa for aa in amino_acids]
    for res1 in amino_acids:
        for res2 in amino_acids:
            if res1 == res2:
                continue
            aa_strings.append(res1+res2)
            
    organism_order = []
    labels = []
    organism_labels_df, proteome_to_info = get_organisms(run_scrambled_analysis)

    lcd_secondary_percentofproteins_df = {aa_string:[] for aa_string in aa_strings}
    domain_specific_secondary_percentofproteins_df = {domain:{aa_string:[] for aa_string in aa_strings} for domain in domains}

    lcd_secondaryprots_df = {}
    
    lcds_per_residue_perc_of_proteome_df = {}
    number_of_lcds_df = {}
            
    secondary_lcds_percentofproteins_df = {aa_string:[] for aa_string in aa_strings}
    lcd_counts_df = {}
    proteome_length_df = {}
    total_prots_df = {}
    window_size = 20
    counter = 0
    initiate = False
    num_skipped = 0
    for domain in domains:
        if run_scrambled_analysis:
            files = list(os.listdir('./' + domain + '_SCRAMBLED'))
        else:
            files = list(os.listdir('./' + domain))
        files = [x for x in files if x.startswith('UP0') and x.endswith('.fasta') and '_DNA' not in x and '_additional' not in x]
        for file in tqdm(files):
            counter += 1
            proteome_name, ext = file.split('.')

            proteome_length_df, total_prots_df, total_aas = calc_proteome_stats(domain, file, proteome_name, proteome_length_df, total_prots_df, window_size, run_scrambled_analysis)
            
            lcd_counts_df, lcd_secondaryprots_df, lcds_per_residue_perc_of_proteome_df, number_of_lcds_df = get_LCD_data(domain, proteome_name, lcd_counts_df, lcd_secondaryprots_df, aa_strings, lcds_per_residue_perc_of_proteome_df, number_of_lcds_df, proteome_length_df, run_scrambled_analysis)

            domain_label = organism_labels_df[proteome_name]
            for aa_string in aa_strings:
                 perc_of_protein = len(lcd_secondaryprots_df[proteome_name][aa_string]) / total_prots_df[proteome_name] *100
                 lcd_secondary_percentofproteins_df[aa_string].append(perc_of_protein)
                 domain_specific_secondary_percentofproteins_df[domain_label][aa_string].append(perc_of_protein)
                 
            organism_order.append(proteome_name)
            labels.append(domain_label)

    if run_scrambled_analysis:
        output = open('TableS1_LCDfrequency_data_NumberOfProtsWithLCDs_SCRAMBLED.tsv', 'w')
        output2 = open('TableS9_LCDfrequency_data_PercentageOfProtsWithLCDs_SCRAMBLED.tsv', 'w')
        output3 = open('LCDfrequency_data_PERCENTILE_SCORE_PercProtsWithLCDs_SCRAMBLED.tsv', 'w')
        output4 = open('LCDfrequency_data_MEAN_DIFFERENCE_PercProtsWithLCDs_SCRAMBLED.tsv', 'w')
        output5 = open('LCDfrequency_data_ProtIdentitiesWithLCDs_SCRAMBLED.tsv', 'w')
        output6 = open('LCDfrequency_data_PerResidue_LCDoccupancy_SCRAMBLED.tsv', 'w')
        output7 = open('LCDfrequency_data_TotalNumberOfLCDs_SCRAMBLED.tsv', 'w')
    else:
        output = open('TableS1_LCDfrequency_data_NumberOfProtsWithLCDs.tsv', 'w')
        output2 = open('TableS9_LCDfrequency_data_PercentageOfProtsWithLCDs.tsv', 'w')
        output3 = open('LCDfrequency_data_PERCENTILE_SCORE_PercProtsWithLCDs.tsv', 'w')
        output4 = open('LCDfrequency_data_MEAN_DIFFERENCE_PercProtsWithLCDs.tsv', 'w')
        output5 = open('LCDfrequency_data_ProtIdentitiesWithLCDs.tsv', 'w')
        output6 = open('LCDfrequency_data_PerResidue_LCDoccupancy.tsv', 'w')
        output7 = open('LCDfrequency_data_TotalNumberOfLCDs.tsv', 'w')

    header = ['Proteome', 'Taxon ID', 'Domain', 'Scientific Name', 'Common Name (when available)', 'Total Proteins'] + aa_strings
    output.write('\t'.join(header) + '\n')
    output2.write('\t'.join(header) + '\n')
    output3.write('\t'.join(header) + '\n')
    output4.write('\t'.join(header) + '\n')
    output5.write('\t'.join(['Proteome', 'Taxon ID', 'Domain', 'Scientific Name', 'Common Name (when available)', 'Total Proteins'] + aa_strings[20:]) + '\n')
    output6.write('\t'.join(header) + '\n')
    output7.write('\t'.join(header) + '\n')
    
    for i, proteome in tqdm(enumerate(organism_order)):
        row_data = [proteome] + proteome_to_info[proteome] + [str(total_prots_df[proteome])] + [len(lcd_secondaryprots_df[proteome][aa_string]) for aa_string in aa_strings]
        output.write('\t'.join([str(x) for x in row_data]) + '\n')
        
        row_data = [proteome] + proteome_to_info[proteome] + [str(total_prots_df[proteome])] + [lcd_secondary_percentofproteins_df[aa_string][i] for aa_string in aa_strings]
        output2.write('\t'.join([str(x) for x in row_data]) + '\n')
        
        domain_label = organism_labels_df[proteome]
        row_data = [proteome] + proteome_to_info[proteome] + [str(total_prots_df[proteome])]
        row_data2 = [proteome] + proteome_to_info[proteome] + [str(total_prots_df[proteome])]
        for aa_string in aa_strings:
            scores = domain_specific_secondary_percentofproteins_df[domain_label][aa_string]
            percentile_val = str(stats.percentileofscore(scores, lcd_secondary_percentofproteins_df[aa_string][i]))
            row_data += [percentile_val]
            
            class_ave = statistics.mean(scores)
            mean_diff = str( lcd_secondary_percentofproteins_df[aa_string][i] - class_ave )
            row_data2 += [mean_diff]
            
        output3.write('\t'.join([str(x) for x in row_data]) + '\n')
        output4.write('\t'.join([str(x) for x in row_data2]) + '\n')
        
        row_data = [proteome] + proteome_to_info[proteome] + [str(total_prots_df[proteome])] + ['_'.join(lcd_secondaryprots_df[proteome][aa_string]) if len(lcd_secondaryprots_df[proteome][aa_string]) > 0 else 'N/A' for aa_string in aa_strings[20:]]
        output5.write('\t'.join([str(x) for x in row_data]) + '\n')
        
        
        row_data = [proteome] + proteome_to_info[proteome] + [str(total_prots_df[proteome])] + [lcds_per_residue_perc_of_proteome_df[proteome][aa_string] for aa_string in aa_strings]
        output6.write('\t'.join([str(x) for x in row_data]) + '\n')

        row_data3 = [proteome] + proteome_to_info[proteome] + [str(total_prots_df[proteome])] + [number_of_lcds_df[proteome][aa_string] for aa_string in aa_strings]
        output7.write('\t'.join([str(x) for x in row_data3]) + '\n')

    output.close()
    output2.close()
    output3.close()
    output4.close()
    output5.close()
    output6.close()
    output7.close()
   
    
def get_organisms(run_scrambled_analysis):

    h = open('OrganismKey_TabDelimited_2022.tsv')
    header = h.readline()
    df = {}
    proteome_to_info = {}
    for line in h:
        proteome_name, taxonid, domain, scientific_name, common_name, synonyms, organism_code = line.rstrip().split('\t')
        if run_scrambled_analysis:
            proteome_name += '_SCRAMBLED'
        df[proteome_name] = domain
        proteome_to_info[proteome_name] = [taxonid, domain, scientific_name, common_name]
        
    h.close()
    
    return df, proteome_to_info
    
        
def calc_proteome_stats(domain, proteome, proteome_name, proteome_length_df, total_prots_df, window_size, run_scrambled_analysis):
    
    total_aas = 0
    if run_scrambled_analysis:
        h = open(os.path.join('.', domain + '_SCRAMBLED', proteome))
    else:
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
    
    
def get_LCD_data(domain, proteome_name, lcd_counts_df, lcd_secondaryprots_df, aa_strings, lcds_per_residue_perc_of_proteome_df, number_of_lcds_df, proteome_length_df, run_scrambled_analysis):
    
    lcd_counts_df[proteome_name] = {aa:0 for aa in amino_acids}
    lcd_secondaryprots_df[proteome_name] = {aa_string:[] for aa_string in aa_strings}
    lengths = {aa:0 for aa in amino_acids}

    lcds_df = {aa_string:[] for aa_string in aa_strings+list(amino_acids)}
    lcds_per_residue_perc_of_proteome_df[proteome_name] = {aa_string:[] for aa_string in aa_strings}
    number_of_lcds_df[proteome_name] = {aa_string:0 for aa_string in aa_strings}
    
    if run_scrambled_analysis:
        h = open(os.path.join('.', domain + '_SCRAMBLED', proteome_name + '_LCDcomposer_RESULTS.tsv'))
    else:
        h = open(os.path.join('.', domain, proteome_name + '_LCDcomposer_RESULTS.tsv'))
    for i in range(11):
        header = h.readline()

    lcd_secondary_prots = {aa:set() for aa in aa_strings}
    
    for line in h:
        prot, uniprot_id, seq, bounds, lcd_class, final_comp, final_disp, *comps = line.rstrip().split('\t')
        if not run_scrambled_analysis:
            junk, prot, *junk = prot.split('|')
            
        lcd_secondary_prots[lcd_class].add(prot)
        lengths[lcd_class] += len(seq)
        lcd_counts_df[proteome_name][lcd_class] += 1
        
        lcds_df[lcd_class].append(seq)
        number_of_lcds_df[proteome_name][lcd_class] += 1
        
    h.close()
    
    if run_scrambled_analysis:
        h = open(os.path.join('.', domain + '_SCRAMBLED', proteome_name + '_LCDcomposer_SecondaryLCDs_RESULTS.tsv'))
    else:
        h = open(os.path.join('.', domain, proteome_name + '_LCDcomposer_SecondaryLCDs_RESULTS.tsv'))

    for i in range(11):
        header = h.readline()
    
    for line in h:
        prot, uniprot_id, seq, bounds, lcd_class, final_comp, final_disp, *comps = line.rstrip().split('\t')
        if not run_scrambled_analysis:
            junk, prot, *junk = prot.split('|')
        res1, res2 = lcd_class

        lcd_secondary_prots[lcd_class].add(prot)
        lcds_df[lcd_class].append(seq)
        number_of_lcds_df[proteome_name][lcd_class] += 1
        
    h.close()

    for res1 in amino_acids:

        lcd_secondaryprots_df[proteome_name][res1] = lcd_secondary_prots[res1]
        for res2 in amino_acids:
            if res1 == res2:
                continue

            lcd_secondaryprots_df[proteome_name][res1+res2] = lcd_secondary_prots[res1+res2]
            
    for aa_string in aa_strings:
        lcds = lcds_df[aa_string]
        lcd_lengths = [len(lcd) for lcd in lcds]
        combined_lcd_length = sum(lcd_lengths)
        perc_proteome_as_lcd = combined_lcd_length / proteome_length_df[proteome_name] * 100
        
        lcds_per_residue_perc_of_proteome_df[proteome_name][aa_string] = perc_proteome_as_lcd
        
    
    return lcd_counts_df, lcd_secondaryprots_df, lcds_per_residue_perc_of_proteome_df, number_of_lcds_df
    
    
def get_args(arguments):

    parser = argparse.ArgumentParser(description='Run main calculations for all proteome LCD statistics', prog='LCD-Composer stats')

    parser.add_argument('-s', '--scrambled_analysis', action='store_true',
                        help="""When this flag is used, calculations are performed on the scrambled proteomes rather than the native proteomes.""")
         
    args = parser.parse_args(arguments)
    
    return args

if __name__ == '__main__':
    import sys, argparse
    args = get_args(sys.argv[1:])
    main(args)
