

def main():

    proteomes_to_analyze = get_proteomes_for_GOanalysis()
    proteomes_with_goa_files = get_proteomes_with_goa_files()
    
    output = open('RUN_Eukaryotes_HQprots_with_HrichPrimary_Excluded.bat', 'w')
    total_analyses = 0
    for proteome in proteomes_to_analyze:
        if proteome not in proteomes_with_goa_files:
            continue
        
        output.write('python run_GOanalyses_by_Eukaryotes_HrichPrimary_excluded.py ' + proteome + '\n')
        total_analyses += 1
        
    output.close()

    
def get_proteomes_with_goa_files():

    h = open('GOAfiles_to_Proteomes_MAP_WITH_SCIENTIFIC_NAMES.tsv')
    header = h.readline()
    goa_proteomes = []
    for line in h:
        items = line.rstrip().split('\t')
        if items[0] != 'Eukaryota':
            continue
            
        proteome = items[3]
        goa_proteomes.append(proteome.replace('.fasta', ''))
        
    h.close()
    
    return goa_proteomes
    

def get_proteomes_for_GOanalysis():

    h = open('HQprots_Without_HrichLCDs.tsv')
    header = h.readline()
    
    proteomes_to_analyze = []
    
    for line in h:
        items = line.rstrip().split('\t')
        if len(items) < 2:
            continue
        proteome = items[0]
        proteomes_to_analyze.append(proteome)
        
    h.close()
        
    return proteomes_to_analyze
    

if __name__ == '__main__':
    main()