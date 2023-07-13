
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
ordered_aas_df = {'Archaea':'LVEAGISDRTKPFNYQMHWC',
                'Bacteria':'LAGVEISDRTKPFQNYMHWC',
                'Eukaryota':'LSAEVGKRTPDINQFYHMCW',
                'Viruses':'LASVGETKDIRNPFQYMHWC'}
domains = ['Archaea','Bacteria','Eukaryota','Viruses']
domain_to_tablenum = {'Archaea':'TableS2', 'Bacteria':'TableS3', 'Eukaryota':'TableS4', 'Viruses':'TableS5'}
domain_to_fignum = {'Archaea':'FigS1A', 'Bacteria':'FigS1B', 'Eukaryota':'FigS1C', 'Viruses':'FigS1D'}

aa_strings = []
for res1 in amino_acids:
    for res2 in amino_acids:
        if res1 == res2:
            continue
        aa_strings.append(res1+res2)

def main():

    for domain in domains:
        ordered_aas = ordered_aas_df[domain]
        cat_df = get_matrix(domain, ordered_aas)
        
        df = pd.DataFrame.from_dict(cat_df)
        df['Secondary LCD Class'] = list(ordered_aas)
        df.set_index('Secondary LCD Class', inplace=True)
        df = df.transpose()

        plotting(df, domain)
        

def plotting(df, domain):

    pal = sns.color_palette('coolwarm', n_colors=16)
    pal = [pal[i] for i in range(0, 16, 3)]

    cg = sns.heatmap(df, vmin=0, vmax=5, cmap=pal)
    ax = plt.gca()
    colorbar = ax.collections[0].colorbar 
    r = colorbar.vmax - colorbar.vmin 
    colorbar.set_ticks([colorbar.vmin + r / 6 * (0.5 + i) for i in range(6)])
    labels = ['Absent', 'Very Rare', 'Rare', 'Normal', 'Common', 'Very Common']
    colorbar.set_ticklabels(labels) 

    plt.xticks(rotation=0, fontname='Arial', fontsize=10)
    plt.yticks(rotation=0, fontname='Arial', fontsize=10)
    plt.xlabel('Primary LCD Class', fontname='Arial', fontsize=12)
    plt.ylabel('Secondary LCD Class', fontname='Arial', fontsize=12)
    
    plt.ylim(20, -0.5)

    fig = plt.gcf()
    plt.savefig(domain_to_fignum[domain] + '_' + domain + '_OrganismLevel_LCDclassRarity_CategoricalHeatmap.tif', bbox_inches='tight', dpi=600, pil_kwargs={'compression':'tiff_lzw'})
    plt.close()
        

def get_matrix(domain, ordered_aas):

    h = open(domain_to_tablenum[domain] + '_' + domain + '_LCDclassRarity_CategoricalHeatmap.tsv')

    df = {aa:[] for aa in ordered_aas}
    categorical_matrix = []
    header = h.readline()
    index = 0
    for line in h:
        sec_aa, *data = line.rstrip().split('\t')
        data = [float(x) if x != '0.0' else 0 for x in data]
        cat_row = []
        for val in data:
            if val == 0:
                cat = 0
            elif val < 5:
                cat = 1
            elif 5 <= val < 20:
                cat = 2
            elif 20 <= val < 50:
                cat = 3
            elif 50 <= val < 75:
                cat = 4
            elif val >= 75:
                cat = 5
            else:
                print('Unexpected value!', val)
                exit()

            cat_row.append(cat)
            df[ordered_aas[index]].append(cat)
        categorical_matrix.append(cat_row)
        index += 1

    return df
    
        
if __name__ == '__main__':
    main()
