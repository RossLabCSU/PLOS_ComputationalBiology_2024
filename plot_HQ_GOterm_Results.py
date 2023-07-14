
import os
import matplotlib.pyplot as plt
import seaborn as sns

all_proteomes = ['UP000006548_3702', 'UP000009136_9913', 'UP000001940_6239', 'UP000002254_9615', 'UP000000437_7955', 'UP000002195_44689', 'UP000000803_7227', 'UP000000539_9031', 'UP000005640_9606', 'UP000000589_10090', 'UP000002311_559292', 'UP000008227_9823', 'UP000002494_10116']
all_abbrevs = ['Athaliana', 'Btaurus', 'Celegans', 'Clupusfamiliaris', 'Drerio', 'Ddiscoideum', 'Dmelanogaster', 'Ggallusdomesticus', 'Hsapiens', 'Mmusculus', 'Scerevisiae', 'Sscrofadomesticus', 'Rnorvegicus']
proteomes_to_abbrevs = {all_proteomes[i]:all_abbrevs[i] for i in range(len(all_proteomes))}

amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
aa_strings = []
for res1 in amino_acids:
    for res2 in amino_acids:
        if res1 == res2:
            continue
        aa_strings.append(res1+res2)
        
def main():

    goterm_freq_threshold = 15

    output = open('TableS7_Summarized_GOtermResults_HQproteins.tsv', 'w')
    counts, goterms, output = get_GO_results('HQprots_GOterm_RESULTS.tsv', goterm_freq_threshold, output, 'All HQ Proteins')
    counts_Xexcluded, goterms_Xexcluded, output = get_GO_results('HQprots-WITHOUT-QX-XQ-HX-XH-LCDs_GOterm_RESULTS.tsv', goterm_freq_threshold, output, 'QX-XQ-HX-XH Proteins Excluded')
    counts_allQexcluded, goterms_allQexcluded, output = get_GO_results('_HQprots-WITHOUT-QrichLCD_GOterm_RESULTS.tsv', goterm_freq_threshold, output, 'Proteins with Q-rich Primary LCDs Excluded')
    
    output.close()
    
    df = {'Count':[],
        'GO term':[],
        'Category':[]}
        
    df = structure_data(counts, goterms, df, goterm_freq_threshold, 'All HQ Proteins')
    df = structure_data(counts_Xexcluded, goterms_Xexcluded, df, goterm_freq_threshold, 'Proteins Also Containing\nQX/XQ/HX/XH LCD Excluded')
    df = structure_data(counts_allQexcluded, goterms_allQexcluded, df, goterm_freq_threshold, 'Proteins Also Containing\nQ-rich LCD Excluded')

    plotting(df, goterm_freq_threshold)
    
    
def structure_data(counts, goterms, df, goterm_freq_threshold, category):

    filtered_goterms = [goterms[i] for i in range(len(goterms)) if counts[i] >= goterm_freq_threshold]
    filtered_counts = [count for count in counts if count >= goterm_freq_threshold]
    
    for i, goterm in enumerate(filtered_goterms):
        count = filtered_counts[i]
        df['Count'].append(count)
        df['GO term'].append(goterm)
        df['Category'].append(category)
        
    return df
    
    
def plotting(df, goterm_freq_threshold):

    colors = ['#014122', '#018141', '#00C060']
    colors = sns.color_palette('Greens', 10)
    colors = colors[::-1]
    colors = [colors[0], colors[3], colors[6]]

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
    plt.xlabel('Enriched GO term', fontname='Arial', fontsize=12)
    plt.ylabel('# of Organisms with GO Term\nSignificantly Enriched', fontname='Arial', fontsize=12)

    leg = plt.legend(handletextpad=0.2, prop={'family':'Arial', 'size':10})

    fig = plt.gcf()
    fig.set_size_inches(12, 4)
    plt.savefig('Fig3B_HQ_GOterm_Results.tif', bbox_inches='tight', dpi=600, pil_kwargs={'compression':'tiff_lzw'})
    plt.close()


def get_GO_results(file_ending, goterm_freq_threshold, output, data_label):

    minimum_depth = 4
    all_enriched_terms = []
    org_count = 0
    
    for file in os.listdir('.'):
        if not file.endswith(file_ending):
            continue
            
        file_info = file.split('_')
        uniprot_acc = file_info[0] + '_' + file_info[1]
        organism = file_info[2]
            
        org_count += 1
        h = open(file)
        header = h.readline().rstrip().split('\t')
        header.insert(0, 'Protein Set for GO-term Analysis')
        header.insert(0, 'Organism')
        header.insert(0, 'UniProt Organism ID')
        if data_label == 'All HQ Proteins' and org_count == 1:
            output.write('\t'.join(header) + '\n')
        
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
                output.write('\t'.join([uniprot_acc, organism, data_label] + items) + '\n')
        h.close()
        
        all_enriched_terms += enriched_terms
                
        fileinfo = file.replace('_HQprots_GOterm_RESULTS.tsv', '')
        
    counts = []
    goterms = []
    for goterm in set(all_enriched_terms):
        counts.append(all_enriched_terms.count(goterm))
        goterms.append(goterm)
        
    counts, goterms = zip(*sorted(zip(counts, goterms), reverse=True))
    
    return counts, goterms, output
    

if __name__ == '__main__':
    main()