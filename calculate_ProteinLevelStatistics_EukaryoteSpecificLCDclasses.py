
amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
domains = ['Archaea','Bacteria','Eukaryota','Viruses']

aa_strings = []
for res1 in amino_acids:
    for res2 in amino_acids:
        if res1 == res2:
            continue
        aa_strings.append(res1+res2)
        
def main():

    num_prots_threshold = 10

    h = open('DomainSpecific_Common_LCDclasses.txt')

    rare_classes_df = {}
    rare_class_indices_df = {}
    for i in range(4):
        domain = domains[i]
        items = h.readline().rstrip().split(': ')
        if len(items) > 1:
            junk, rare_classes = items
            rare_classes = rare_classes.split(', ')
        else:
            rare_classes = []
        rare_classes_df[domain] = rare_classes
        indices = [aa_strings.index(rare_class) for rare_class in rare_classes]
        rare_class_indices_df[domain] = indices
    h.close()

    h = open('TableS1_LCDfrequency_data_NumberOfProtsWithLCDs.tsv')
    header = h.readline()
        
    lcdprot_count_df = {}
    lcdprot_perc_df = {}
    domains_df = {domain:set() for domain in domains}
    proteome_to_orginfo = {}
    for line in h:
        proteome, taxid, domain, sci_name, common_name, total_prots, *data = line.rstrip().split('\t')
        data = data[20:]
        total_prots = int(total_prots)
        
        for i, lcd_class in enumerate(rare_classes_df[domain]):
            class_index = rare_class_indices_df[domain][i]
            num_lcd_prots = int(data[class_index])
            
            lcdprot_count_df[proteome] = lcdprot_count_df.get(proteome, {})
            lcdprot_count_df[proteome][lcd_class] = num_lcd_prots
            
            lcdprot_perc_df[proteome] = lcdprot_perc_df.get(proteome, {})
            lcdprot_perc_df[proteome][lcd_class] = num_lcd_prots / total_prots * 100

            domains_df[domain].add( proteome )
            proteome_to_orginfo[proteome] = (sci_name, common_name)


    domain = 'Eukaryota'
    output = open('ProteinLevelStatistics_EukaryoteSpecificLCDclasses.tsv', 'w')
    output.write('\t'.join(['Rare LCD Class', 'Number of Organisms Passing Minimum-Protein-Number Threshold', 'Maximum Number of Proteins in Single Organism', 'Organism(s) with Maximum Number of Proteins', 'Maximum Percentage of Proteins with LCD in Single Organism', 'Organism(s) with Maximum Percentage of Proteins', 'All Organisms Passing Minimum Number of Proteins Threshold', 'Protein Counts for Each Organism Passing Minimum Number of Proteins Threshold', 'Protein Count, all organisms']) + '\n')
    proteome_list = sorted(list(domains_df[domain]))
    
    counts = []
    for rare_lcd_class in rare_classes_df[domain]:
        num_organisms_passing_threshold = sum([1 for proteome in proteome_list if rare_lcd_class in lcdprot_count_df[proteome] and lcdprot_count_df[proteome][rare_lcd_class] >= num_prots_threshold])
        all_organisms_passing_threshold = [proteome for proteome in proteome_list if rare_lcd_class in lcdprot_count_df[proteome] and lcdprot_count_df[proteome][rare_lcd_class] >= num_prots_threshold]
        protcounts_all_organisms_passing_threshold = [lcdprot_count_df[proteome][rare_lcd_class] for proteome in proteome_list if rare_lcd_class in lcdprot_count_df[proteome] and lcdprot_count_df[proteome][rare_lcd_class] >= num_prots_threshold]
        all_protcounts = [lcdprot_count_df[proteome][rare_lcd_class] for proteome in proteome_list if rare_lcd_class in lcdprot_count_df[proteome]]
        counts = [lcdprot_count_df[proteome][rare_lcd_class] if rare_lcd_class in lcdprot_count_df[proteome] else 0 for proteome in proteome_list]
        max_prots = max(counts)
        max_prot_locs = [1 if count == max_prots else 0 for count in counts]
        max_prot_proteomes = [proteome_list[i] for i in range(len(proteome_list)) if max_prot_locs[i] == 1]
        
        percs = [lcdprot_perc_df[proteome][rare_lcd_class] if rare_lcd_class in lcdprot_perc_df[proteome] else 0 for proteome in proteome_list]
        max_perc = max(percs)
        max_perc_locs = [1 if perc == max_perc else 0 for perc in percs]
        max_perc_proteomes = [proteome_list[i] for i in range(len(proteome_list)) if max_perc_locs[i] == 1]
        
        output.write('\t'.join([rare_lcd_class, str(num_organisms_passing_threshold), str(max_prots), '; '.join(['(' + x + ', ' + proteome_to_orginfo[x][0] + ', ' + proteome_to_orginfo[x][1] + ')' for x in max_prot_proteomes]), str(max_perc), '; '.join(['(' + x + ', ' + proteome_to_orginfo[x][0] + ', ' + proteome_to_orginfo[x][1] + ')' for x in max_perc_proteomes]), '; '.join(['(' + x + ', ' + proteome_to_orginfo[x][0] + ', ' + proteome_to_orginfo[x][1] + ')' for x in all_organisms_passing_threshold]), '; '.join([str(x) for x in protcounts_all_organisms_passing_threshold]), '; '.join([str(x) for x in all_protcounts])]) + '\n')
        
    output.close()
            
if __name__ == '__main__':
    main()