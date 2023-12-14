
amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
ordered_aas_df = {'Archaea':'LVEAGISDRTKPFNYQMHWC',
                'Bacteria':'LAGVEISDRTKPFQNYMHWC',
                'Eukaryota':'LSAEVGKRTPDINQFYHMCW',
                'Viruses':'LASVGETKDIRNPFQYMHWC'}
domains = ['Archaea','Bacteria','Eukaryota','Viruses']

def main():

    category = 'PerResidueLCDoccupancy'
    for domain in domains:
        ordered_aas = ordered_aas_df[domain]
        file = 'TableS11_PerResidueLevelStatistics_Top3values_for_LCDclasses.tsv'
        
        df, orgs_df = get_data(file, category, domain)
        mat, orgs_mat = construct_matrix(df, orgs_df, category, ordered_aas)
        output_matrix(domain, category, mat, 'ActualValues', ordered_aas)
        output_matrix(domain, category, orgs_mat, 'CorrespondingOrganisms', ordered_aas)
            
def construct_matrix(df, orgs_df, category, ordered_aas):
    
    mat = []
    orgs_mat = []
    for res2 in ordered_aas:
        row = []
        orgs_row = []
        for res1 in ordered_aas:
            val = df[res1][res2]
            org = orgs_df[res1][res2]
            row.append(val)
            orgs_row.append(org)
            
        mat.append(row)
        orgs_mat.append(orgs_row)
        
    return mat, orgs_mat
    
    
def output_matrix(domain, category, mat, filetype, ordered_aas):
    
    output = open(domain + '_' + category.replace(' ', '') + '_Matrix_' + filetype + '.tsv', 'w')
    output.write('\t'.join(['Secondary LCD Class'] + list(ordered_aas)) + '\n')
    for i, row in enumerate(mat):
        output.write('\t'.join([ordered_aas[i]] + [str(x) for x in row]) + '\n')
    output.close()
  
    
def get_data(file, category, domain_of_interest):
    
    h = open(file)
    header = h.readline()
    
    df = {aa:{res:0 for res in amino_acids} for aa in amino_acids}
    orgs_df = {aa:{res:0 for res in amino_acids} for aa in amino_acids}
    organisms = []
    for line in h:
        lcd_class, domain, main_val, val_cat, proteome, sci_name, common_name, total_prots, total_aas, mean_lcd_length = line.rstrip().split('\t')
        
        if domain != domain_of_interest:
            continue
            
        main_val_list = main_val.split('; ')
        main_val = float(main_val_list[0])
        
        sciname_list = sci_name.split('; ')
        sci_name = sciname_list[0]
        if sci_name == 'Sepia pharaonis':
            sci_name = sciname_list[1]
            main_val = float(main_val_list[1])
        
        if len(lcd_class) == 1:
            res1 = lcd_class
            res2 = lcd_class
        else:
            res1 = lcd_class[0]
            res2 = lcd_class[1]
            
        df[res1][res2] = main_val
        orgs_df[res1][res2] = sci_name
        
    h.close()
    
    return df, orgs_df

    
if __name__ == '__main__':
    main()
