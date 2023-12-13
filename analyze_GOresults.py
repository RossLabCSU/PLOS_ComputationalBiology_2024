
import os
from tqdm import tqdm
domains = ['Archaea', 'Bacteria', 'Eukaryota', 'Viruses']     
amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
all_lcd_classes = list(amino_acids)
for res1 in amino_acids:
    for res2 in amino_acids.replace(res1, ''):
        all_lcd_classes.append( res1+res2 )
        
cwd = os.getcwd()

def main():

    df = initialize_df()
    org_count = {'Archaea':294,
                 'Bacteria':7533,
                 'Eukaryota':2061,
                 'Viruses':6134}

    for domain in domains:
        dirname = os.path.join(cwd, domain + '_GOresults')
        
        for file in tqdm(os.listdir(dirname)):
            if not file.endswith('RESULTS.tsv'):
                continue
                
            file_info = file.split('_')
            org_id = file_info[0] + '_' + file_info[1]

            lcd_class = file_info[2].replace('class', '')
                
            h = open(os.path.join(dirname, file))
            header = h.readline()
            
            for line in h:
                items = line.rstrip().split('\t')
                go_id = items[0]
                go_category = items[1]
                enr_or_pur = items[2]
                go_desc = items[3].strip()
                hits_in_pop = int(items[7])
                totalprots = int(items[8])
                depth = int(items[12])
                perc_hits_in_total = str(round(hits_in_pop/totalprots*100, 3))
                if enr_or_pur != 'e' or depth < 4:
                    continue
                p_sidak = float(items[15])
                go_info = '\t'.join([go_id, go_category, go_desc]) 
                if p_sidak < 0.05:
                    df[domain][lcd_class][go_info] = df[domain][lcd_class].get(go_info, [])
                    df[domain][lcd_class][go_info].append( (p_sidak, org_id, perc_hits_in_total) )

            h.close()

    output_data(df, org_count, domain)
    

def output_data(df, org_count, domain):

    output = open('TableS7_GOterm_Summary_AllOrganisms_AllClasses.tsv', 'w')
    output.write('\t'.join(['Domain of Life', 'LCD Class', 'GO term ID', 'GO term Category', 'GO term Description', 'Number of Organisms with GO term Significantly Enriched', 'Total Number of Organisms Analyzed', 'Percentage of Organisms with GO term Significantly Enriched', 'Sidak-corrected p-values (semicolon-delimited)', 'Corresponding Organism UniProt IDs (semicolon-delimited)', 'Percentage of GO term-Associated Proteins in Background Proteome']) + '\n')
    for domain in domains:
        total_orgs = org_count[domain]
        for lcd_class in all_lcd_classes:
            for go_info in df[domain][lcd_class]:
                data = df[domain][lcd_class][go_info]
                data = sorted(data)
                pvals = [str(tup[0]) for tup in data]
                org_ids = ';'.join([tup[1] for tup in data])
                perc_hits = ';'.join([tup[2] for tup in data])
                num_sig = len(pvals)
                pvals = ';'.join(pvals)
                output.write('\t'.join([domain, lcd_class, go_info, str(num_sig), str(total_orgs), str(num_sig/total_orgs*100), pvals, org_ids, perc_hits]) + '\n')
    output.close()

        
def initialize_df():
    df = {}
    for domain in domains:
        df[domain] = {}
        for res1 in amino_acids:
            for res2 in amino_acids:
                if res1 == res2:
                    df[domain][res1] = {}
                else:
                    df[domain][res1+res2] = {}
                    
    return df
    

if __name__ == '__main__':
    main()
