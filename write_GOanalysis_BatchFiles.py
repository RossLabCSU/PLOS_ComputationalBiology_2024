
import os
amino_acids = 'ACDEFGHIKLMNPQRSTVWY'

def main():
    
    domains = ['Archaea', 'Bacteria', 'Eukaryota', 'Viruses']
    
    for domain in domains:
        proteomes = get_proteomes(domain)
        output = open('Run_' + domain + '_GOanalyses.bat', 'w')
        
        for proteome in proteomes:
            output.write('python run_GOanalyses_by_Organism_All_LCD_Classes.py ' + proteome + '\n')
                
        output.close()


def get_proteomes(domain_of_interest):
    
    h = open('GOAfiles_to_Proteomes_MAP_WITH_SCIENTIFIC_NAMES.tsv')
    
    header = h.readline()
    proteomes = []
    for line in h:
        items = line.rstrip().split('\t')
        domain = items[0]
        if domain != domain_of_interest:
            continue
        proteome = items[3]
        proteomes.append(proteome.replace('.fasta', ''))
        
    h.close()
    
    return proteomes
    

if __name__ == '__main__':
    main()