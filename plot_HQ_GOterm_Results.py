
import os
import matplotlib.pyplot as plt
import seaborn as sns

def main():

    goterm_freq_threshold = 30

    output = open('TableS7_Summarized_GOtermResults_HQproteins.tsv', 'w')
    counts, goterms, go_categories, output = get_GO_results('Eukaryota_GOterm_Summary_AllOrganisms_AllClasses.tsv', goterm_freq_threshold, output, 'All HQ Proteins')
    counts_Xexcluded, goterms_Xexcluded, go_categories_Xexcluded, output = get_GO_results('HQ_GOresults_QX-XQ-HX-XH_excluded_GOterm_Summary.tsv', goterm_freq_threshold, output, 'QX-XQ-HX-XH Proteins Excluded')
    counts_allQexcluded, goterms_allQexcluded, go_categories_allQexcluded, output = get_GO_results('HQ_GOresults_QrichPrimary_excluded_GOterm_Summary.tsv', goterm_freq_threshold, output, 'Proteins with Q-rich Primary LCDs Excluded')
    counts_allHexcluded, goterms_allHexcluded, go_categories_allHexcluded, output = get_GO_results('HQ_GOresults_HrichPrimary_excluded_GOterm_Summary.tsv', goterm_freq_threshold, output, 'Proteins with H-rich Primary LCDs Excluded')
    
    output.close()
    
    df = {'Count':[],
        'GO term':[],
        'Category':[]}
        
    df = structure_data(counts, goterms, df, goterm_freq_threshold, 'All HQ Proteins')
    df = structure_data(counts_allHexcluded, goterms_allHexcluded, df, goterm_freq_threshold, 'Proteins Also Containing\nH-rich LCD Excluded')
    df = structure_data(counts_Xexcluded, goterms_Xexcluded, df, goterm_freq_threshold, 'Proteins Also Containing\nQX/XQ/HX/XH LCD Excluded')
    df = structure_data(counts_allQexcluded, goterms_allQexcluded, df, goterm_freq_threshold, 'Proteins Also Containing\nQ-rich LCD Excluded')

    plotting(df, goterm_freq_threshold, go_categories)
    
    
def structure_data(counts, goterms, df, goterm_freq_threshold, category):

    if category == 'All HQ Proteins':
        filtered_goterms = [goterms[i] for i in range(len(goterms)) if counts[i] >= goterm_freq_threshold]
        filtered_counts = [count for count in counts if count >= goterm_freq_threshold]
    else:
        filtered_goterms = [goterm for goterm in goterms if goterm in df['GO term']]
        filtered_counts = [count for i, count in enumerate(counts) if goterms[i] in df['GO term']]

    
    for i, goterm in enumerate(filtered_goterms):
        count = filtered_counts[i]
        df['Count'].append(count)
        df['GO term'].append(goterm)
        df['Category'].append(category)
        
    return df
    
    
def plotting(df, goterm_freq_threshold, go_categories):

    tick_colors = ['#be0119', '#00035b', '#0b5509'] # SCARLET, DARK BLUE, FOREST
    cats = ['BP', 'MF', 'CC']
    tick_cmap = {cat:tick_colors[i] for i, cat in enumerate(cats)}
    
    colors = sns.color_palette('Greens', 10)
    colors = colors[::-1]
    colors = [colors[0], colors[2], colors[4], colors[6]]

    sns.barplot(x='GO term', y='Count', data=df, hue='Category', palette=colors)
    ax = plt.gca()
    orig_labels = ax.get_xticklabels()
    labels = []
    textwrap = 20
    for text_obj in orig_labels:
        term = text_obj.get_text()
        if len(term) > textwrap:
            first_half = term[:textwrap]
            second_half = term[textwrap:]
            if ' ' in second_half:
                ind = second_half.index(' ')
                second_half = second_half[:ind] + '\n' + second_half[ind+1:]
                new_label = first_half + second_half
            else:
                new_label = term
        else:
            new_label = term
        labels.append(new_label)

    plt.xticks([x for x in range(len(labels))], labels=labels, rotation=90, fontname='Arial', fontsize=10)
    plt.yticks(fontname='Arial', fontsize=10)
    
    ax = plt.gca()
    xticklabels = ax.get_xticklabels()
    for i, label in enumerate(go_categories[:len(xticklabels)]):
        xticklabels[i].set_color(tick_cmap[label])
    
    plt.xlabel('Enriched GO term', fontname='Arial', fontsize=12)
    plt.ylabel('# of Organisms with GO Term\nSignificantly Enriched', fontname='Arial', fontsize=12)

    leg = plt.legend(handletextpad=0.2, prop={'family':'Arial', 'size':10})

    fig = plt.gcf()
    fig.set_size_inches(12, 4)
    plt.savefig('Fig5B_HQ_GOterm_Results.tif', bbox_inches='tight', dpi=600, pil_kwargs={'compression':'tiff_lzw'})
    plt.close()
    
    
def get_GO_results(goterm_file, goterm_freq_threshold, output, data_label):

    # minimum_depth = 4     # RESTRICTION IMPOSED WHEN GATHERING GO-TERM RESULTS.
    all_enriched_terms = []

    h = open(goterm_file)
    header = h.readline().rstrip().split('\t')
    header.insert(0, 'Protein Set for GO-term Analysis')

    if data_label == 'All HQ Proteins':
        output.write('\t'.join(header) + '\n')

    counts = []
    goterms = []
    go_categories = []
    for line in h:
        items = line.rstrip().split('\t')
        lcd_class = items[1]
        if lcd_class != 'HQ':
            continue
        goterm = items[4].rstrip()
        num_orgs_significant = int(items[5])
        goterms.append(goterm)
        counts.append(num_orgs_significant)
        go_categories.append(items[3])
        output.write('\t'.join([data_label] + items) + '\n')
    h.close()

    counts, goterms, go_categories = zip(*sorted(zip(counts, goterms, go_categories), reverse=True))
    
    return counts, goterms, go_categories, output
    

if __name__ == '__main__':
    main()