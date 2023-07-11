
import pandas as pd
amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
ordered_aas_df = {'Archaea':'LVEAGISDRTKPFNYQMHWC',
                'Bacteria':'LAGVEISDRTKPFQNYMHWC',
                'Eukaryota':'LSAEVGKRTPDINQFYHMCW',
                'Viruses':'LASVGETKDIRNPFQYMHWC'}

total_orgcount_df = {'Archaea':340,
                    'Bacteria':8623,
                    'Eukaryota':2119,
                    'Viruses':10788}
                    
domains = ['Archaea','Bacteria','Eukaryota','Viruses']
domain_to_tablenum = {'Archaea':'TableS2', 'Bacteria':'TableS3', 'Eukaryota':'TableS4', 'Viruses':'TableS5'}

aa_strings = [aa+aa for aa in amino_acids]
for res1 in amino_acids:
    for res2 in amino_acids:
        if res1 == res2:
            continue
        aa_strings.append(res1+res2)

total_organisms = 21870

def main():

    org_count_df = get_data()

    for domain in domains:
        ordered_aas = ordered_aas_df[domain]
        df = {aa:[] for aa in ordered_aas}
        total_orgcount = total_orgcount_df[domain]
        ordered_aas = ordered_aas_df[domain]

        for res1 in ordered_aas:
            for i, res2 in enumerate(ordered_aas):
                aa_string = res1 + res2
                perc = org_count_df[domain][aa_string] / total_orgcount *100
                df[res1].append(perc)

        df = pd.DataFrame.from_dict(df)
        df['Secondary LCD Class'] = list(ordered_aas)
        df.set_index('Secondary LCD Class', inplace=True)
        
        df.to_csv(domain_to_tablenum[domain] + '_' + domain + '_LCDclassRarity_CategoricalHeatmap.tsv', sep='\t')


def get_data():

    h = open('TableS1_LCDfrequency_data_NumberOfProtsWithLCDs.tsv')
    header = h.readline()
    org_count_df = {domain:{aa_string:0 for aa_string in aa_strings} for domain in domains}
    
    for line in h:
        proteome, taxon, domain, sci_name, common_name, total_prots, *prot_counts = line.rstrip().split('\t')
        prot_counts = [int(x) for x in prot_counts]
        for i, count in enumerate(prot_counts):
            aa_string = aa_strings[i]
            if count > 0:
                org_count_df[domain][aa_string] += 1
    h.close()
    
    return org_count_df
    
        
if __name__ == '__main__':
    main()