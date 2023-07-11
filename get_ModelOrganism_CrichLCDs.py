
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
class_indices = [aa_strings.index(lcd_class) for lcd_class in lcd_classes]
        
def main():

    h = open('LCDfrequency_data_ProtIdentitiesWithLCDs.tsv')
    header = h.readline()
        
    output = open('ModelOrganisms_CrichLCDproteins.tsv', 'w')
    output.write('\t'.join(['Proteome', 'Taxon ID', 'Domain', 'Scientific Name', 'Common Name (when available)', 'Total Proteins'] + lcd_classes + [c + ' Protein IDs' for c in lcd_classes]) + '\n')
    
    for line in h:
        proteome, taxid, domain, sci_name, common_name, total_prots, *protid_lists = line.rstrip().split('\t')
        if proteome in all_proteomes:
            class_protids = [protid_lists[class_index].split('_') for class_index in class_indices]
            class_protids = [x if x != ['N/A'] else [] for x in class_protids]
            class_protids_strs = [protid_lists[class_index] for class_index in class_indices]
            class_freqs = [str(len(protids)) for protids in class_protids]

            output.write('\t'.join([proteome, taxid, domain, sci_name, common_name, total_prots] + class_freqs + class_protids_strs) + '\n')

    h.close()
    output.close()
    

if __name__ == '__main__':
    main()