
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import statistics

ordered_aas_df = {'Archaea':'LVEAGISDRTKPFNYQMHWC',
                'Bacteria':'LAGVEISDRTKPFQNYMHWC',
                'Eukaryota':'LSAEVGKRTPDINQFYHMCW',
                'Viruses':'LASVGETKDIRNPFQYMHWC'}
                
amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
aa_strings = [aa for aa in amino_acids]
for res1 in amino_acids:
    for res2 in amino_acids:
        if res1 == res2:
            continue
        aa_strings.append(res1+res2)

domains = ['Archaea','Bacteria','Eukaryota','Viruses']
fig_labels = {'Archaea':'FigS8A', 'Bacteria':'FigS8B', 'Eukaryota':'FigS8C', 'Viruses':'FigS8-NotShown'}

def main():

    lnORs_df = get_data()
    for domain in domains:
        matrix = []
        ordered_aas = ordered_aas_df[domain]
        df = {aa:[] for aa in ordered_aas}
        for res1 in ordered_aas:
            row = []
            for res2 in ordered_aas:
                if res1 == res2:
                    lnORs = lnORs_df[domain][res1]
                else:
                    lnORs = lnORs_df[domain][res1+res2]

                lnORs = [float(x) if x != 'N/A' else 1.0 for x in lnORs]   # CONVERT "N/A" VALUES TO HIGH P-VALUE (1.0) WHICH WON'T PASS THE SIGNIFICANCE THRESHOLD

                if len(lnORs) == 0:
                    med_lnOR = 0
                else:
                    med_lnOR = statistics.median(lnORs)
                row.append(med_lnOR)
                df[res1].append(med_lnOR)
            matrix.append(row)

        df = pd.DataFrame.from_dict(df)
        df['Secondary LCD Class'] = list(ordered_aas)
        df.set_index('Secondary LCD Class', inplace=True)

        plot_heatmap(df, domain)
        df.to_csv(domain + '_SecondaryLCDs_HeatmapMatrix_Median-lnORs.tsv', sep='\t')


def get_data():

    file = 'Observed_vs_Scrambled_FisherExact_Results.tsv'
    h = open(file)
    header = h.readline().rstrip().split('\t')
    df = {domain:{lcd_class:[] for lcd_class in aa_strings} for domain in domains}

    for line in h:
        items = line.rstrip().split('\t')
        if len(items) < 15:
            continue
        domain = items[1]
        lcd_class = items[2]
        pval = items[-1]
        obs_num = int(items[3])
        scr_num = int(items[4])
        total_prots = int(items[5])
        
        lnOR = items[7]
        if lnOR == 'N/A':
            if items[14] == 'N/A':
                continue
            lnOR = float(items[14])     # USE BIASED lnOR IF THE ACTUAL lnOR IS NOT AVAILABLE
        else:
            lnOR = float(lnOR)

        df[domain][lcd_class].append( float(lnOR) )

    h.close()

    return df


def plot_heatmap(df, domain):

    pal = sns.color_palette('gist_heat', as_cmap=True)

    cg = sns.heatmap(df, cmap=pal)
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
    fig_label = fig_labels[domain]
    plt.savefig(fig_label + '_' + domain + '_Median-lnOR_Actual-vs-Scrambled-Proteomes_Heatmap.tif', bbox_inches='tight', dpi=600, pil_kwargs={'compression':'tiff_lzw'})
    plt.close()


if __name__ == '__main__':
    main()