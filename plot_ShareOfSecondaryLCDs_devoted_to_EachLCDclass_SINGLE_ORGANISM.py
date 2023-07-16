
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Patch

all_proteomes = ['UP000005640_9606', 'UP000002311_559292']
fig_labels = ['Fig9B', 'Fig9C']

amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
aa_strings = []
for res1 in amino_acids:
    for res2 in amino_acids:
        if res1 == res2:
            continue
        aa_strings.append(res1+res2)
        
ordered_residues = 'AMILVFWYPGCQNSTHDERK'
        
def main():

    lcd_counts_df = get_data()

    for i, proteome in enumerate(all_proteomes):
        fig_label = fig_labels[i]
        shares_df = {res1:{res2:0 for res2 in amino_acids} for res1 in amino_acids}

        output = open('LCDshares_' + proteome + '.tsv', 'w')
        output.write('\t'.join(['Primary AA', 'Secondary AA', 'LCD Share']) + '\n')
        
        for res1 in amino_acids:
            sec_shares_df = {res2:[] for res2 in amino_acids.replace(res1, '')}
            counts = [lcd_counts_df[proteome][res1+res2] for res2 in amino_acids if res2 != res1]
            total = sum(counts)
            if total == 0:
                continue
            for res2 in amino_acids.replace(res1, ''):
                share = lcd_counts_df[proteome][res1+res2] / total * 100
                sec_shares_df[res2].append(share)
                shares_df[res1][res2] = share
                output.write('\t'.join([res1, res2, str(share)]) + '\n')

        output.close()
                
        sample_sizes = []
        for res1 in ordered_residues:
            n = sum([lcd_counts_df[proteome][res1+res2] for res2 in ordered_residues.replace(res1, '')])
            sample_sizes.append(n)
                
        plotting(shares_df, sample_sizes, proteome, fig_label)
            
        
def plotting(shares_df, sample_sizes, proteome, fig_label):
    
    colors = ['#1f77b4', '#aec7e8', '#ff7f0e', '#ffbb78', '#2ca02c', '#98df8a', '#d62728', '#ff9896', '#9467bd', '#c5b0d5', '#8c564b', '#c49c94', '#e377c2', '#f7b6d2', '#7f7f7f', '#c7c7c7', '#bcbd22', '#dbdb8d', '#17becf', '#9edae65']
    hydrophobics = sns.color_palette('Reds', n_colors=5)
    aromatics = sns.color_palette('Oranges', n_colors=3)
    pro = sns.color_palette('Purples', n_colors=1)
    gly = sns.color_palette('Greys', n_colors=1)
    cys = ['#000000']
    polar = sns.color_palette('Blues', n_colors=4)
    charged = sns.color_palette('Greens', n_colors=5)
    
    colors = hydrophobics[::-1] + aromatics[::-1] + pro + gly + cys + polar[::-1] + charged[::-1]

    for xpos, res1 in enumerate(ordered_residues):
        yvals = [shares_df[res1][res2] for res2 in ordered_residues]

        height_sum = 0
        colors_index = 0
        for yval in yvals:
            plt.bar(x=xpos, height=yval, bottom=height_sum, color=colors[colors_index])
            colors_index += 1
            height_sum += yval

    plt.xticks([x for x in range(len(ordered_residues))], labels=ordered_residues, fontname='Arial', fontsize=12)
    plt.yticks([x for x in range(0, 120, 20)], fontname='Arial', fontsize=12)
    
    plt.xlabel('Primary LCD Category', fontname='Arial', fontsize=14)
    plt.ylabel('% Share for Each Secondary LCD Category', fontname='Arial', fontsize=14)
    
    plt.xlim(-0.5, 19.5)
    plt.ylim(0, 120)
    
    for i, samp_size in enumerate(sample_sizes):
        plt.text(i, 102, 'n = ' + str(samp_size), fontname='Arial', fontsize=10, rotation=90, ha='center')

    leg_items = [Patch(facecolor=colors[i], label=ordered_residues[i]) for i in range(len(ordered_residues))]
    leg = plt.legend(handles=leg_items, prop={'family':'Arial', 'size':10}, loc=2, bbox_to_anchor=(1.0, 1.02), title='Secondary AA', handletextpad=0.2)
    leg.set_title('Secondary AA')
    plt.setp(leg.get_title(), multialignment='center', fontname='Arial', fontsize=14)
        
    fig = plt.gcf()
    fig.set_size_inches(7, 5)
    plt.savefig(fig_label + '_' + proteome + '_SecondaryLCDs_SharePlot.tif', bbox_inches='tight', dpi=600, pil_kwargs={'compression':'tiff_lzw'})
    plt.close()


def get_data():

    lcd_counts_df = {proteome:{} for proteome in all_proteomes}
    for res1 in amino_acids:
        for res2 in amino_acids.replace(res1, ''):

            for proteome in all_proteomes:
                lcd_counts_df[proteome][res1+res2] = 0
                
            h = open(res1+res2 + '_LCDs_AllOrganisms.tsv')
            header = h.readline()
            for line in h:
                proteome, *remainder = line.rstrip().split('\t')
                if proteome in all_proteomes:
                    lcd_counts_df[proteome][res1+res2] += 1
                    
            h.close()

    return lcd_counts_df
    
            
if __name__ == '__main__':
    main()