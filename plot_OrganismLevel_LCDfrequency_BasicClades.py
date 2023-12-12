
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
domains = ['Archaea', 'Bacteria', 'Eukaryota', 'Viruses']
aa_names = {'A':'Alanine (A)', 'C':'Cysteine (C)', 'D':'Aspartic Acid (D)', 'E':'Glutamic Acid (E)',
        'F':'Phenylalanine (F)', 'G':'Glycine (G)', 'H':'Histidine (H)', 'I':'Isoleucine (I)', 
        'K':'Lysine (K)', 'L':'Leucine (L)', 'M':'Methionine (M)', 'N':'Asparagine (N)',
        'P':'Proline (P)', 'Q':'Glutamine (Q)', 'R':'Arginine (R)', 'S':'Serine (S)',
        'T':'Threonine (T)', 'V':'Valine (V)', 'W':'Tryptophan (W)','Y':'Tyrosine (Y)'}
fig_labels = {'Archaea':'FigS5',
                'Bacteria':'FigS6',
                'Eukaryota':'FigS7',
                'Viruses':'FigS8'}

def main():

    lineage_df, taxid_to_common, taxid_to_sci = get_taxonomy_lineages()
    data_df = get_data(lineage_df)
    
    for domain in domains:
        taxids = [taxid for taxid in lineage_df if domain in lineage_df[taxid]]
        lineages = [lineage_df[taxid] for taxid in taxids]
        lineages, taxids = zip(*sorted(zip(lineages, taxids)))
        second_tax_items = []
        for lineage in lineages:
            items = lineage.split(',')
            sec_item = items[1].strip()
            second_tax_items.append(sec_item)
            
        plotting_df = {'Second Taxonomy Label':[],
                        'Per-Residue Occupancy':[],
                        'LCD Class':[],
                        'Sample Sizes':[]}

        samplesize_df = {}

        for j, aa in enumerate(amino_acids):
            for tax_label in sorted(list(set(second_tax_items))):
                euks = [taxids[i] for i, second_tax_item in enumerate(second_tax_items) if second_tax_item == tax_label]
                vals = [data_df[taxid][j] for taxid in euks]
                percentage_with_lcds = sum(vals) / len(vals) * 100
                label_str = tax_label + ' (' + str(len(vals)) + ')'
                plotting_df['Second Taxonomy Label'].append(label_str)
                plotting_df['Per-Residue Occupancy'].append( percentage_with_lcds )
                plotting_df['LCD Class'].append(aa)
                plotting_df['Sample Sizes'].append(len(vals))
                samplesize_df[label_str] = len(vals)

        plotting_df = pd.DataFrame.from_dict(plotting_df)

        subplots_second_taxonomy_items(plotting_df, domain, samplesize_df)

    
def subplots_second_taxonomy_items(plotting_df, domain, samplesize_df):

    aa_index = 1
    fig, axes = plt.subplots(5,4,sharex=True, sharey=False, figsize=(8,10))
    big_subplot = fig.add_subplot(111)
    big_subplot.spines['top'].set_color('none')
    big_subplot.spines['bottom'].set_color('none')
    big_subplot.spines['left'].set_color('none')
    big_subplot.spines['right'].set_color('none')
    big_subplot.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)
    
    samplesizes = []
    xtick_labels = []
    for label in samplesize_df:
        xtick_labels.append(label)
        samplesizes.append(samplesize_df[label])
        
    samplesizes, xtick_labels = zip(*sorted(zip(samplesizes, xtick_labels), reverse=True))
    retained_labels = xtick_labels[:10]

    for aa in amino_acids:

        ax = plt.subplot(5,4, aa_index)
        temp = plotting_df[plotting_df['LCD Class'] == aa]
        temp2 = temp[temp['Second Taxonomy Label'].isin(retained_labels)]
        sns.barplot(x='Second Taxonomy Label', y='Per-Residue Occupancy', data=temp2)

        overall_max = ax.patches[0].get_height()

        plt.ylabel('')
        plt.xlabel('')
        plt.title(aa_names[aa], fontsize=12, fontname='Arial', pad=0.9)
        plt.xticks(rotation=90)

        if aa_index in range(1,17):
            ax.set_xticklabels([])

        aa_index += 1

    fig.text(-0.025, 0.5, 'Percentage of Organisms with' + r'$\geq$' + '1 LCD', va='center', rotation='vertical', fontname='Arial', fontsize=16)
    fig.text(0.5, -0.02, 'Basic Clade', ha='center', fontname='Arial', fontsize=16)
    plt.tight_layout(pad=0.2)
    plt.savefig(fig_labels[domain] + '_' + domain + '_OrganismLevel_LCDfrequency_BasicClades_Top10cladesOnly.tif', bbox_inches ='tight', dpi=600, pil_kwargs={'compression':'tiff_lzw'})
    plt.close()  
    
    
def get_data(files_df):

    h = open('LCDfrequency_data_TotalNumberOfLCDs.tsv')
    header = h.readline()
    
    df = {}
    for line in h:
        proteome, taxid, domain, sci_name, common_name, total_prots, *data = line.rstrip().split('\t')
        if taxid not in files_df:
            continue
        
        data = [int(x) for x in data]
        data = data[:20]    # PRIMARY LCD CLASSES
        contains_LCD_checks = [1 if x > 0 else 0 for x in data]
        df[taxid] = contains_LCD_checks
        
    h.close()
    
    return df
    
    
def get_taxonomy_lineages():

    h = open('UniProt_Taxonomy.tsv')
    header = h.readline()
    
    df = {}
    taxid_to_common = {}
    taxid_to_sci = {}
    for line in h:
        items = line.rstrip().split('\t')
        lineage = items[7]
        taxid = items[0]
        common = items[1]
        sci = items[2]
        taxid_to_common[taxid] = common
        taxid_to_sci[taxid] = sci
        df[taxid] = lineage
        
    h.close()
    
    return df, taxid_to_common, taxid_to_sci
    
    
if __name__ == '__main__':
    main()