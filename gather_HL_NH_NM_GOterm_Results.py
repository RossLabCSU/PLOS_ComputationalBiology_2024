
import os

def main():

    goterm_freq_threshold = 1
    lcd_classes = ['NM', 'NH', 'HL']
    
    output = open('TableS8_Summarized_GOtermResults_HL_NH_NM_proteins.tsv', 'w')

    write_header = True # HEADER IS ONLY WRITTEN THE FIRST TIME GO RESULTS ARE GATHERED

    for lcd_class in lcd_classes:

        counts, goterms, output, write_header = get_GO_results(lcd_class + 'prots_GOterm_RESULTS.tsv', lcd_class, goterm_freq_threshold, output, write_header)
        
    output.close()
    

def get_GO_results(file_ending, lcd_class, goterm_freq_threshold, output, write_header):

    minimum_depth = 0
    all_enriched_terms = []
    org_count = 0

    for file in os.listdir():
        if not file.endswith(file_ending):
            continue
            
        file_info = file.split('_')
        uniprot_acc = file_info[0] + '_' + file_info[1]
        organism = file_info[2]
            
        org_count += 1
        h = open(file)
        header = h.readline().rstrip().split('\t')
        header.insert(0, 'Organism')
        header.insert(0, 'UniProt Organism ID')
        header.insert(0, 'LCD Class')
        
        if write_header:
            output.write('\t'.join(header) + '\n')
            write_header = False
        
        enriched_terms = []
        for line in h:
            items = line.rstrip().split('\t')
            e_or_p = items[2]
            goterm = items[3].rstrip()
            depth = int(items[12])
            if depth < minimum_depth:
                continue
            sidak_pval = float(items[15])
            if sidak_pval < 0.05:
                enriched_terms.append(goterm)
                output.write('\t'.join([lcd_class, uniprot_acc, organism] + items) + '\n')
        h.close()
        
        all_enriched_terms += enriched_terms
                
        fileinfo = file.replace(lcd_class + 'prots_GOterm_RESULTS.tsv', '')
        
    counts = []
    goterms = []
    for goterm in set(all_enriched_terms):
        counts.append(all_enriched_terms.count(goterm))
        goterms.append(goterm)
        
    counts, goterms = zip(*sorted(zip(counts, goterms), reverse=True))
    
    return counts, goterms, output, write_header
    

if __name__ == '__main__':
    main()