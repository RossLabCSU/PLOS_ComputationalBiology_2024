
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
ordered_aas_df = {'Archaea':'LVEAGISDRTKPFNYQMHWC',
                'Bacteria':'LAGVEISDRTKPFQNYMHWC',
                'Eukaryota':'LSAEVGKRTPDINQFYHMCW',
                'Viruses':'LASVGETKDIRNPFQYMHWC'}

domains = ['Archaea','Bacteria','Eukaryota','Viruses']
domain_to_filenum = {'Archaea':'Fig7A', 'Bacteria':'Fig7B', 'Eukaryota':'Fig7C', 'Viruses':'Fig7D'}

def main():

    category = 'PerResidueLCDoccupancy'

    for domain in domains:
        ordered_aas = ordered_aas_df[domain]
        df = pd.read_csv(domain + '_' + category + '_Matrix_ActualValues.tsv', sep='\t')
        df.set_index('Secondary LCD Class', inplace=True)
        
        max_offdiag, ondiag_vals = get_nondiag_max(df, domain, ordered_aas)

        plot_heatmap(df, domain, category, max_offdiag)
        
        df.to_csv(domain + '_SecondaryLCDs_HeatmapMatrix.tsv', sep='\t')
    

def get_nondiag_max(df, domain, ordered_aas):
    
    offdiag_vals = []
    ondiag_vals = []
    for aa in ordered_aas:
        for i, sec_aa in enumerate(ordered_aas):
            if aa != sec_aa:
                offdiag_vals.append(df[aa][i])
            else:
                ondiag_vals.append(df[aa][i])

    return max(offdiag_vals), ondiag_vals
    

def plot_heatmap(df, domain, category, max_offdiag):

    pal = sns.color_palette('mako', as_cmap=True)
    
    mask = []
    offdiag_vals = []
    for i in range(20):
        row = []
        for j in range(20):
            if i == j:
                row.append(1)
            else:
                row.append(0)
                
        mask.append(row)

    cg = sns.heatmap(df, mask=np.array(mask), vmin=0, vmax=max_offdiag, cmap=pal)
    ax = plt.gca()
    ax.set_facecolor("grey")

    cb=ax.collections[0].colorbar
    cb.set_ticklabels(cb.ax.get_yticklabels(), fontsize=10, fontname='Arial') 
    
    plt.xticks(rotation=0, fontname='Arial', fontsize=10)
    plt.yticks(rotation=0, fontname='Arial', fontsize=10)
    plt.xlabel('Primary Amino Acid', fontname='Arial', fontsize=12)
    plt.ylabel('Secondary Amino Acid', fontname='Arial', fontsize=12)
    
    plt.ylim(20, 0)
    fig = plt.gcf()
    plt.savefig(domain_to_filenum[domain] + '_' + domain + '_SecondaryLCDs_PerResidueOccupancy_Heatmap.tif', bbox_inches='tight', dpi=600, pil_kwargs={'compression':'tiff_lzw'})
    plt.close()

        
if __name__ == '__main__':
    main()