
all_proteomes = ['UP000006548_3702', 'UP000009136_9913', 'UP000001940_6239', 'UP000002254_9615', 'UP000000437_7955', 'UP000002195_44689', 'UP000000803_7227', 'UP000000539_9031', 'UP000005640_9606', 'UP000000589_10090', 'UP000002311_559292', 'UP000008227_9823', 'UP000002494_10116']

all_abbrevs = ['Athaliana', 'Btaurus', 'Celegans', 'Clupusfamiliaris', 'Drerio', 'Ddiscoideum', 'Dmelanogaster', 'Ggallusdomesticus', 'Hsapiens', 'Mmusculus', 'Scerevisiae', 'Sscrofadomesticus', 'Rnorvegicus']
proteomes_to_abbrevs = {all_proteomes[i]:all_abbrevs[i] for i in range(len(all_proteomes))}

amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
aa_strings = []
for res1 in amino_acids:
    for res2 in amino_acids:
        if res1 == res2:
            continue
        aa_strings.append(res1+res2)

lcd_classes = ['C' + aa for aa in amino_acids if aa != 'C']
lcd_classes += [aa + 'C' for aa in amino_acids if aa != 'C']

def main():

    num_prots_threshold = 5
    output = open('TableS6_ModelOrganisms_KeratinAssociated_CrichSubset.tsv', 'w')
    output.write('\t'.join(['Proteome', 'Scientific Name', 'Common Name', 'LCD Class', 'Proteins']) + '\n')
    h = open('ModelOrganisms_CrichLCDproteins.tsv')
    header = h.readline()

    for line in h:
        proteome, taxid, domain, sci_name, common_name, total_prots, *data = line.rstrip().split('\t')
        freqs = data[:len(lcd_classes)]
        prots = data[len(lcd_classes):]
        prots = [prots[i].split('_') for i in range(len(prots))]
        for i, lcd_class in enumerate(lcd_classes):
            freq = int(freqs[i])
            if freq >= num_prots_threshold:
                keratin_prots = get_keratin_prots(proteome, lcd_class)
                if len(keratin_prots) > 0:
                    output.write('\t'.join([proteome, sci_name, common_name, lcd_class, ', '.join(keratin_prots)]) + '\n')

    h.close()
    output.close()
    
    
def get_keratin_prots(proteome, lcd_class):
        
    index = all_proteomes.index(proteome)
    abbrev = all_abbrevs[index]
    
    h = open(proteome + '_' + abbrev + '_' + lcd_class + '_GOterm_RESULTS.tsv')
    header = h.readline()
    
    all_prots = []
    
    for line in h:
        items = line.rstrip().split('\t')
        goterm = items[3].strip()
        if goterm in ['keratin filament', 'keratinization']:
            prots = items[-1].split(', ')
            if prots != ['1']:
                all_prots += prots
                
    h.close()
    
    all_prots = set(all_prots)
    
    return all_prots
        
        
if __name__ == '__main__':
    main()