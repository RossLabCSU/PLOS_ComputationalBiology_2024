
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np

ordered_aas_df = {'Archaea':'LVEAGISDRTKPFNYQMHWC',
                'Bacteria':'LAGVEISDRTKPFQNYMHWC',
                'Eukaryota':'LSAEVGKRTPDINQFYHMCW',
                'Viruses':'LASVGETKDIRNPFQYMHWC'}

fig_labels = {('SignificantOnly', '', 'Enriched', 'Archaea'):'Fig6A',
            ('SignificantOnly', '', 'Enriched', 'Bacteria'):'Fig6B',
            ('SignificantOnly', '', 'Enriched', 'Eukaryota'):'Fig6C',
            ('SignificantOnly', '', 'Enriched', 'Viruses'):'Fig6D',
            ('StatisticalSignificanceIgnored', '', 'Enriched', 'Archaea'):'FigS17A',
            ('StatisticalSignificanceIgnored', '', 'Enriched', 'Bacteria'):'FigS17C',
            ('StatisticalSignificanceIgnored', '', 'Enriched', 'Eukaryota'):'FigS17E',
            ('StatisticalSignificanceIgnored', '', 'Enriched', 'Viruses'):'FigS17G',
            ('StatisticalSignificanceIgnored', 'NormalizedScale', 'Depleted', 'Archaea'):'FigS17B',
            ('StatisticalSignificanceIgnored', 'NormalizedScale', 'Depleted', 'Bacteria'):'FigS17D',
            ('StatisticalSignificanceIgnored', 'NormalizedScale', 'Depleted', 'Eukaryota'):'FigS17F',
            ('StatisticalSignificanceIgnored', 'NormalizedScale', 'Depleted', 'Viruses'):'FigS17H',
            ('SignificantOnly', '', 'Depleted', 'Archaea'):'FigS18-NotShown',
            ('SignificantOnly', 'NormalizedScale', 'Depleted', 'Archaea'):'FigS18-NotShown',
            ('SignificantOnly', '', 'Depleted', 'Bacteria'):'FigS18-NotShown',
            ('SignificantOnly', 'NormalizedScale', 'Depleted', 'Bacteria'):'FigS18-NotShown',
            ('SignificantOnly', '', 'Depleted', 'Eukaryota'):'FigS18A',
            ('SignificantOnly', 'NormalizedScale', 'Depleted', 'Eukaryota'):'FigS18B',
            ('SignificantOnly', '', 'Depleted', 'Viruses'):'FigS18-NotShown',
            ('SignificantOnly', 'NormalizedScale', 'Depleted', 'Viruses'):'FigS18-NotShown',
            }
                
amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
aa_strings = [aa for aa in amino_acids]
for res1 in amino_acids:
    for res2 in amino_acids:
        if res1 == res2:
            continue
        aa_strings.append(res1+res2)

domains = ['Archaea','Bacteria','Eukaryota','Viruses']

def main():

    max_df = {}
    
    for use_statistical_significance in ['SignificantOnly', 'StatisticalSignificanceIgnored']:
        for use_normed_scale in ['', 'NormalizedScale']:
            for cat in ['Enriched', 'Depleted']:
                pvals_df = get_data(cat)
                for domain in domains:
                    matrix = []
                    ordered_aas = ordered_aas_df[domain]
                    df = {aa:[] for aa in ordered_aas}
                    for res1 in ordered_aas:
                        row = []
                        for res2 in ordered_aas:
                            if res1 == res2:
                                pvals = pvals_df[domain][res1]
                            else:
                                pvals = pvals_df[domain][res1+res2]

                            pvals = [float(x) if x != 'N/A' else 1.0 for x in pvals]   # CONVERT "N/A" VALUES TO HIGH P-VALUE (1.0) WHICH WON'T PASS THE SIGNIFICANCE THRESHOLD

                            if len(pvals) == 0:
                                perc_sig = 0
                            else:
                                if use_statistical_significance == 'SignificantOnly':
                                    perc_sig = sum([1 for pval in pvals if pval < 0.05]) / len(pvals) *100
                                else:
                                    perc_sig = sum([1 for pval in pvals if pval < 1.0-0.00000001]) / len(pvals) *100
                            row.append(perc_sig)
                            df[res1].append(perc_sig)
                        matrix.append(row)

                    df = pd.DataFrame.from_dict(df)
                    df['Secondary LCD Class'] = list(ordered_aas)
                    df.set_index('Secondary LCD Class', inplace=True)
                    
                    max_offdiag, ondiag_vals = get_nondiag_max(df, domain, ordered_aas)

                    if use_normed_scale == 'NormalizedScale':
                        if cat == 'Enriched':
                            max_df[domain] = max_offdiag
                        else:
                            max_offdiag = max_df[domain]   # OVERRIDES max_offdiag FOR THE "Depleted" CATEGORY SO THAT EACH PAIR OF HEATMAPS IS ON THE SAME SCALE

                    # NO NEED TO PLOT NORMALIZED SCALE FOR ENRICHED CATEGORY. THE NORMALIZATION ONLY APPLIES TO THE 'Depleted' CATEGORY BECAUSE IT IS NORMALIZED TO THE CORRESPONDING 'Enriched' CATEGORY (WHICH ALWAYS HAS HIGHER VALUES). IF PLOTS WERE GENERATED FOR THE 'NormalizedScale' AND 'Enriched' COMBINATION, THEY WOULD BE IDENTICAL TO THOSE GENERATED WITHOUT USING THE NORMALIZED SCALE FOR THE 'Enriched' CATEGORY.
                    if use_normed_scale == 'NormalizedScale' and cat == 'Enriched':
                        continue
                
                    plot_heatmap(df, domain, max_offdiag, cat, use_normed_scale, use_statistical_significance)
                    
                    param_tup = (use_statistical_significance, use_normed_scale, cat, domain)
                    if param_tup in fig_labels:
                        fig_label = fig_labels[param_tup]
                    else:
                        fig_label = 'Extra'
                        
                    df.to_csv(fig_label + '_Data_' + domain + '_' + cat + '_' + use_normed_scale + '_' + use_statistical_significance + '.tsv', sep='\t')


def get_data(cat):

    h = open('Observed_vs_Scrambled_FisherExact_Results.tsv')
    header = h.readline()
    df = {domain:{lcd_class:[] for lcd_class in aa_strings} for domain in domains}
    
    for line in h:
        items = line.rstrip().split('\t')
        domain = items[1]
        lcd_class = items[2]
        pval = items[12]
        obs_num = int(items[3])
        scr_num = int(items[4])

        if cat == 'Enriched':
            if scr_num > obs_num:
                df[domain][lcd_class].append('1.0') # APPEND NON-SIGNIFICANT P-VALUE IF THE CATEGORY IS LOOKING FOR LCDS MORE COMMON IN ACTUAL PROTEOMES COMPARED TO SCRAMBLED PROTEOMES
            else:
                df[domain][lcd_class].append(pval)
        elif cat == 'Depleted':
            if obs_num > scr_num:
                df[domain][lcd_class].append('1.0') # APPEND NON-SIGNIFICANT P-VALUE IF THE CATEGORY IS LOOKING FOR LCDS MORE LESS IN ACTUAL PROTEOMES COMPARED TO SCRAMBLED PROTEOMES
            else:
                df[domain][lcd_class].append(pval)
        
    h.close()
    
    return df
    

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
    

def plot_heatmap(df, domain, max_offdiag, cat, use_normed_scale, use_statistical_significance):

    pal = sns.color_palette('gist_heat', as_cmap=True)

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
    param_tup = (use_statistical_significance, use_normed_scale, cat, domain)
    if param_tup in fig_labels:
        fig_label = fig_labels[param_tup]
    else:
        fig_label = 'Extra'
    plt.savefig(fig_label + '_' + domain + '_Original-vs-Scrambled-Proteomes_Heatmap_' + cat + '_' + use_normed_scale + '_' + use_statistical_significance + '.tif', bbox_inches='tight', dpi=600)
    plt.close()


if __name__ == '__main__':
    main()
