
import math
import pandas as pd
from tqdm import tqdm
import pickle

domains = ['Archaea', 'Bacteria', 'Eukaryota']

def main():

    method = 'Manhattan'
    
    data_df, map_df, all_proteomes = get_lcd_arrays()
    lineage_df, taxid_to_common, taxid_to_sci, all_clades, all_domains = get_taxonomy_lineages()

    proteomes = []
    for proteome in all_proteomes:
        leader, taxid = proteome.split('_')
        if taxid in lineage_df:
            proteomes.append(proteome)

    domains_df = {domain1:{domain2:{'dist sum':0, 'n':0} for domain2 in domains} for domain1 in domains}
    clades_df = {clade1:{clade2:{'dist sum':0, 'n':0} for clade2 in all_clades} for clade1 in all_clades}
    
    df = {'Proteome1':[],
        'Proteome2':[],
        'Distance':[],
        'Domain1':[],
        'Domain2':[],
        'Clade1':[],
        'Clade2':[]}
        
    clade_tracker = set()

    output = open('AllOrganisms_SecondaryLCD_PerResidueOccupancy_Pairwise' + method + 'Distances_FULL_CALCULATIONS.tsv', 'w')
    output2 = open('AllOrganisms_SecondaryLCD_PerResidueOccupancy_Pairwise' + method + 'Distances_LOGSCALE_FULL_CALCULATIONS.tsv', 'w')
    output_longformat = open('AllOrganisms_SecondaryLCD_PerResidueOccupancy_Pairwise' + method + 'Distances_LONG_FORMAT_FULL_CALCULATIONS.tsv', 'w')
    
    output.write('\t'.join([''] + proteomes) + '\n')
    output2.write('\t'.join([''] + proteomes) + '\n')
    output_longformat.write('\t'.join(['Proteome 1', 'Domain of Life 1', 'Clade 1', 'Proteome 2', 'Domain of Life 2', 'Clade 2', 'Distance (' + method + ')']) + '\n')
    
    missing_count = 0
    i = 0
    total_output_longformat = 0
    
    for proteome1 in tqdm(proteomes):
        leader, taxid = proteome1.split('_')
        lineage1 = lineage_df[taxid]
        domain1 = lineage1[0]
        clade1 = lineage1[1]
        clade_tracker.add(clade1)
        
        row = []
        row2 = []
        for j, proteome2 in enumerate(proteomes):
            leader, taxid2 = proteome2.split('_')
            lineage2 = lineage_df[taxid2]
            domain2 = lineage2[0]
            clade2 = lineage2[1]
            vals1 = data_df[proteome1]
            vals2 = data_df[proteome2]
            dist = calc_Distance(vals1, vals2)
            row.append(str(dist))
            if dist == 0:
                row2.append(str(-1000))
            else:
                row2.append(str(math.log10(dist)))

            total_output_longformat += 1
            output_longformat.write('\t'.join([proteome1, domain1, clade1, proteome2, domain2, clade2, str(dist)]) + '\n')
            
            domains_df[domain1][domain2]['dist sum'] += dist
            domains_df[domain1][domain2]['n'] += 1

            clades_df[clade1][clade2]['dist sum'] += dist
            clades_df[clade1][clade2]['n'] += 1

        output.write('\t'.join([proteome1] + row) + '\n')
        output2.write('\t'.join([proteome1] + row2) + '\n')
        
        i += 1
        
    output.close()
    output2.close()
    output_longformat.close()

    pf = open('Saved_domains_df_for_MeanPairwise' + method + 'Distance_FULL_CALCULATIONS.dat', 'wb')
    pickle.dump(domains_df, pf)
    pf.close()
    
    pf = open('Saved_clades_df_for_MeanPairwise' + method + 'Distance_FULL_CALCULATIONS.dat', 'wb')
    pickle.dump(clades_df, pf)
    pf.close()
    
    calc_mean_dists(domains_df, domains, method, 'Domain', all_domains)
    calc_mean_dists(clades_df, all_clades, method, 'Clade', all_domains)
    
    
def calc_mean_dists(df, categories, method, data_type, all_domains):
    
    df = pd.DataFrame.from_dict(df)
    output = open(data_type + 's_MeanPairwise' + method + 'Distance_FULL_CALCULATIONS.tsv', 'w')
    if data_type == 'Clade':
        output.write('\t'.join(['', ''] + all_domains) + '\n')
    output.write('\t'.join(['', ''] + categories) + '\n')
    for i, cat1 in enumerate(categories):
        if data_type == 'Clade':
            row = [all_domains[i], cat1]
        else:
            row = [cat1]
        for cat2 in categories:
            dist_sum = df[cat1][cat2]['dist sum']
            n = df[cat1][cat2]['n']
            if n == 0:
                row.append(str(-1000))
            else:
                row.append(str(dist_sum / n))
        output.write('\t'.join(row) + '\n')
    output.close()
    
    
def get_taxonomy_lineages():

    h = open('UniProt_Taxonomy.tsv')
    header = h.readline()
    
    df = {}
    taxid_to_common = {}
    taxid_to_sci = {}
    all_clades = []
    all_domains = []
    for line in h:
        items = line.rstrip().split('\t')
        lineage = items[7].replace('"', '').split(', ')
        domain = lineage[0]
        if domain == 'Viruses':
            continue
        clade = lineage[1]
        if clade not in all_clades:
            all_clades.append(clade)
            all_domains.append(domain)
        taxid = items[0]
        common = items[1]
        sci = items[2]
        taxid_to_common[taxid] = common
        taxid_to_sci[taxid] = sci
        df[taxid] = lineage
        
    h.close()
    
    return df, taxid_to_common, taxid_to_sci, all_clades, all_domains
    
    

def get_lcd_arrays():

    h = open('LCDfrequency_data_PerResidue_LCDoccupancy.tsv')
    header = h.readline()
    
    df = {}
    map_df = {}
    proteomes = []
    for line in h:
        proteome, taxid, domain, sciname, commonname, totalprots, *data = line.rstrip().split('\t')
        if domain == 'Viruses':
            continue
        df[proteome] = [float(x) for x in data[20:]]    # SECONDARY LCD CLASSES ONLY
        map_df[proteome] = [taxid, domain, sciname, commonname, totalprots]
        proteomes.append(proteome)
        
    h.close()
    
    return df, map_df, proteomes


def calc_Distance(vals1, vals2):

    dist = sum( [ abs( vals1[i] - vals2[i] ) for i in range(len(vals1)) ] )

    return dist

    
if __name__ == '__main__':
    main()