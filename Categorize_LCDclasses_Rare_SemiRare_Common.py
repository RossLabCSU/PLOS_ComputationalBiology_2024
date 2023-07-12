
domains = ['Archaea','Bacteria','Eukaryota','Viruses']
amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
aa_strings = [aa for aa in amino_acids]
for res1 in amino_acids:
    for res2 in amino_acids:
        if res1 == res2:
            continue
        aa_strings.append(res1+res2)
        
def main():

    labels = ['Archaea', 'Bacteria', 'Eukaryota', 'Viruses']
    data = [set(), set(), set(), set()]
    
    perc_threshold = 5
    common_perc_threshold = 15
    extremely_rare_perc_threshold = 2
    
    df, organism_info_lol, organism_data_lol = get_LCD_data_DOMAIN_SPECIFIC()
    rare_classes, rare_class_indicator, rare_classes_no_instances, rare_noinstance_combined, common_classes, remaining_classes, extremely_rare_classes = get_common_LCDclasses_DOMAIN_SPECIFIC(df, perc_threshold, organism_data_lol, common_perc_threshold, extremely_rare_perc_threshold)
    
    output = open('Common_LCDclasses.txt', 'w')
    for i, domain in enumerate(labels):
        lcd_classes = common_classes[domain]
        output.write(domain + ' Common LCD Classes: ' + ', '.join(lcd_classes) + '\n')
    output.close()
    
    output = open('SemiRare_LCDclasses.txt', 'w')
    for i, domain in enumerate(labels):
        lcd_classes = remaining_classes[domain]
        output.write(domain + ' SemiRare LCD Classes: ' + ', '.join(lcd_classes) + '\n')
    output.close()
    
    output = open('Rare_LCDclasses.txt', 'w')
    for i, domain in enumerate(labels):
        lcd_classes = rare_classes[domain]
        output.write(domain + ' Rare LCD Classes: ' + ', '.join(lcd_classes) + '\n')
    output.close()
    
    output = open('ZeroInstance_LCDclasses.txt', 'w')
    for i, domain in enumerate(labels):
        lcd_classes = rare_classes_no_instances[domain]
        output.write(domain + ' ZeroInstance LCD Classes: ' + ', '.join(lcd_classes) + '\n')
    output.close()
                
    output = open('DomainSpecific_Common_LCDclasses.txt', 'w')
    for i, domain in enumerate(labels):
        lcd_classes = common_classes[domain]
        other_domains = [x for x in domains if x != domain]
        rare_classes_other_domains = [extremely_rare_classes[d] for d in other_domains]
        filtered_lcd_classes = []
        for c in lcd_classes:
            checks_other_domains = [1 if c in domain_classes else 0 for domain_classes in rare_classes_other_domains]
            if sum(checks_other_domains) == len(other_domains):
                filtered_lcd_classes.append(c)
        output.write(domain + ' DomainSpecific LCD Classes: ' + ', '.join(filtered_lcd_classes) + '\n')
    output.close()
    
    
def get_common_LCDclasses_DOMAIN_SPECIFIC(df, perc_threshold, organism_data_lol, common_perc_threshold, extremely_rare_perc_threshold):
    
    rare_classes = {domain:[] for domain in domains}
    rare_class_indicators = {domain:{} for domain in domains}
    rare_classes_no_instances = {domain:[] for domain in domains}
    rare_noinstance_combined = {domain:[] for domain in domains}
    common_classes = {domain:[] for domain in domains}
    remaining_classes = {domain:[] for domain in domains}
    extremely_rare_classes = {domain:[] for domain in domains}
    for domain in df:
        total_organisms = len(organism_data_lol[domain])
        num_orgs_threshold = total_organisms * perc_threshold/100 +0.000001
        extremely_rare_num_orgs_threshold = total_organisms * extremely_rare_perc_threshold/100 -0.000001
        common_num_orgs_threshold = total_organisms * common_perc_threshold/100 -0.000001
        for aa_string in aa_strings:
            nonzeros_indicator = [1 if val != 0 else 0 for val in df[domain][aa_string]]
            nonzeros_count = sum(nonzeros_indicator)

            # STORE COMMON LCD CLASSES
            if nonzeros_count >= common_num_orgs_threshold:
                common_classes[domain].append(aa_string)
                
            if num_orgs_threshold<nonzeros_count<common_num_orgs_threshold:
                remaining_classes[domain].append(aa_string)
                
            # STORE RARE LCD CLASSES
            if 0 < nonzeros_count <= num_orgs_threshold:
                rare_classes[domain].append(aa_string)
                rare_class_indicators[domain][aa_string] = nonzeros_indicator
                rare_noinstance_combined[domain].append(aa_string)
                
            # STORE EXTREMELY RARE LCD CLASSES
            if nonzeros_count < extremely_rare_num_orgs_threshold:
                extremely_rare_classes[domain].append(aa_string)
                
            # STORE ZERO-INSTANCE LCD CLASSES
            if nonzeros_count == 0:
                rare_classes_no_instances[domain].append(aa_string)
                rare_noinstance_combined[domain].append(aa_string)
            
    return rare_classes, rare_class_indicators, rare_classes_no_instances, rare_noinstance_combined, common_classes, remaining_classes, extremely_rare_classes

    
def get_LCD_data_DOMAIN_SPECIFIC():

    h = open('TableS1_LCDfrequency_data_NumberOfProtsWithLCDs.tsv')
    header = h.readline()
    df = {domain:{aa_string:[] for aa_string in aa_strings} for domain in domains}
    organism_info = {domain:[] for domain in domains}
    organism_data = {domain:[] for domain in domains}
    
    for line in h:
        items = line.rstrip().split('\t')
        domain = items[2]
        info = items[:6]
        data = [int(x) for x in items[6:]]     
        for i, val in enumerate(data):
            df[domain][aa_strings[i]].append(val)

        organism_info[domain].append(info)
        organism_data[domain].append(data)
            
    h.close()
    
    return df, organism_info, organism_data
    

if __name__ == '__main__':
    main()