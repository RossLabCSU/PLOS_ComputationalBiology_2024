
amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
import os
from tqdm import tqdm
domains = ['Archaea', 'Bacteria', 'Eukaryota', 'Viruses']

def main():
            
    organism_order = []
    labels = []
    organism_labels_df, proteome_to_info = get_organisms()
    total_organisms_df = {domain:0 for domain in domains}
    lcd_class_organismcount_df = {domain:{aa:0 for aa in amino_acids} for domain in domains}
    
    window_size = 20
    counter = 0
    secondary_threshold = 20
    
    for domain in domains:
        files = list(os.listdir('./' + domain))
        files = [x for x in files if x.startswith('UP0') and x.endswith('.fasta') and '_DNA' not in x and '_additional' not in x]
        for file in tqdm(files):
            counter += 1
            proteome_name, ext = file.split('.')

            domain_label = organism_labels_df[proteome_name]
            total_organisms_df[domain_label] += 1
            
            lcd_class_organismcount_df = get_LCD_data(proteome_name, domain_label, lcd_class_organismcount_df)

    output = open('TableS6_PrimaryLCDs_NumberOfOrganisms.tsv', 'w')
    output.write('\t'.join(['Domain of Life', 'LCD Class', 'Number of Organisms with LCD', 'Total # of Organisms in Domain of Life', 'Percentage of Organisms in Domain of Life with LCD']) + '\n')
    for domain in domains:
        for aa in amino_acids:
            num_orgs_with_lcd = lcd_class_organismcount_df[domain][aa]
            total_orgs = total_organisms_df[domain]
            perc_with_lcd = num_orgs_with_lcd / total_orgs *100
            output.write('\t'.join([domain, aa, str(num_orgs_with_lcd), str(total_orgs), str(perc_with_lcd)]) + '\n')
            
    output.close()
    
def get_organisms():

    h = open('OrganismKey_TabDelimited_2022.tsv')
    header = h.readline()
    df = {}
    proteome_to_info = {}
    for line in h:
        proteome_name, taxonid, domain, scientific_name, common_name, synonyms, organism_code = line.rstrip().split('\t')
        df[proteome_name] = domain
        proteome_to_info[proteome_name] = [taxonid, domain, scientific_name, common_name]
        
    h.close()
    
    return df, proteome_to_info


def get_LCD_data(proteome_name, domain_label, lcd_class_organismcount_df):

    h = open(os.path.join('.', domain_label, proteome_name + '_LCDcomposer_RESULTS.tsv'))
    for i in range(11):
        header = h.readline()

    lcd_class_set = set()
    
    for line in h:
        prot, uniprot_id, seq, bounds, lcd_class, final_comp, final_disp, *comps = line.rstrip().split('\t')
        
        lcd_class_set.add(lcd_class)
        
    h.close()
    
    for res in lcd_class_set:
        lcd_class_organismcount_df[domain_label][res] += 1

    return lcd_class_organismcount_df
    
    
if __name__ == '__main__':
    main()