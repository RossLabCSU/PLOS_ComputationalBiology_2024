
import os

all_proteomes = ['UP000006548_3702', 'UP000009136_9913', 'UP000001940_6239', 'UP000002254_9615', 'UP000000437_7955', 'UP000002195_44689', 'UP000000803_7227', 'UP000000539_9031', 'UP000005640_9606', 'UP000000589_10090', 'UP000002311_559292', 'UP000008227_9823', 'UP000002494_10116']
all_abbrevs = ['Athaliana', 'Btaurus', 'Celegans', 'Clupusfamiliaris', 'Drerio', 'Ddiscoideum', 'Dmelanogaster', 'Ggallusdomesticus', 'Hsapiens', 'Mmusculus', 'Scerevisiae', 'Sscrofadomesticus', 'Rnorvegicus']

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

    num_prots_threshold = 1
    output = open('TableS7_ModelOrganisms_KeratinAssociated_CrichSubset.tsv', 'w')
    output.write('\t'.join(['Proteome', 'Scientific Name', 'Common Name', 'LCD Class', 'Proteins']) + '\n')

    uniprot_to_sciname, uniprot_to_commonname = get_proteome_map()

    for proteome in all_proteomes:
        for i, lcd_class in enumerate(lcd_classes):
            keratin_prots = get_keratin_prots(proteome, lcd_class)
            if len(keratin_prots) > 0:
                sci_name = uniprot_to_sciname[proteome]
                common_name = uniprot_to_commonname[proteome]
                output.write('\t'.join([proteome, sci_name, common_name, lcd_class, ', '.join(keratin_prots)]) + '\n')

    output.close()
    
    
def get_keratin_prots(proteome, lcd_class):
        
    index = all_proteomes.index(proteome)
    abbrev = all_abbrevs[index]
    
    try:
        h = open(os.path.join(os.getcwd(), 'Eukaryota_GOresults', proteome + '_' + lcd_class + 'class_LCDprots_GOterm_RESULTS.tsv'))
    except:

        return []
    header = h.readline()
    
    all_prots = []
    
    for line in h:
        items = line.rstrip().split('\t')
        goterm = items[3].strip()
        if goterm in ['keratin filament', 'keratinization']:
            prots = items[-1].strip().split(', ')
            if prots != ['1']:
                all_prots += prots
                
    h.close()
    
    all_prots = set(all_prots)
    
    return all_prots
        
    
def get_proteome_map():

    h = open('TableS1_LCDfrequency_data_NumberOfProtsWithLCDs.tsv')
    header = h.readline()
    uniprot_to_sciname = {}
    uniprot_to_commonname = {}
    for line in h:
        proteome, taxid, domain, sciname, commonname, *remainder = line.rstrip().split('\t')
        uniprot_to_sciname[proteome] = sciname
        uniprot_to_commonname[proteome] = commonname
        
    h.close()
    
    return uniprot_to_sciname, uniprot_to_commonname


if __name__ == '__main__':
    main()