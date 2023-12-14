
import matplotlib.pyplot as plt
import seaborn as sns
import statistics
from matplotlib.patches import Patch
import os

domains = ['Eukaryota']

amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
aa_strings = []
for res1 in amino_acids:
    for res2 in amino_acids:
        if res1 == res2:
            continue
        aa_strings.append(res1+res2)
        
ordered_residues = 'AMILVFWYPGCQNSTHDERK'
        
def main():

    proteomes_df = get_proteomes()

    for domain in domains:

        domain_proteomes = proteomes_df[domain]

        lcd_counts_df = get_data(domain, domain_proteomes)

        shares_df = {res1:{res2:0 for res2 in amino_acids} for res1 in amino_acids}

        output = open('LCDshares_' + domain + '.tsv', 'w')
        output.write('\t'.join(['Primary AA', 'Secondary AA', 'LCD Share']) + '\n')
        for res1 in amino_acids:
            sec_shares_df = {res2:[] for res2 in amino_acids.replace(res1, '')}
            for proteome in domain_proteomes:
                counts = [lcd_counts_df[proteome][res1+res2] for res2 in amino_acids if res2 != res1]
                total = sum(counts)
                if total == 0:
                    continue
                for res2 in amino_acids.replace(res1, ''):
                    share = lcd_counts_df[proteome][res1+res2] / total * 100
                    sec_shares_df[res2].append(share)

            for i, res2 in enumerate(amino_acids.replace(res1, '')):
                if len(sec_shares_df[res2]) > 0:
                    mean_share = statistics.mean(sec_shares_df[res2])
                else:
                    mean_share = 0
                shares_df[res1][res2] = mean_share
                output.write('\t'.join([res1, res2, str(mean_share)]) + '\n')
        output.close()

        sample_sizes = []
        for res1 in ordered_residues:
            n = 0
            for proteome in domain_proteomes:
                n += sum([lcd_counts_df[proteome][res1+res2] for res2 in ordered_residues.replace(res1, '')])
            sample_sizes.append(n)
        
        plotting(shares_df, sample_sizes, domain)
            
        
def plotting(shares_df, sample_sizes, domain):

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

    if domain in ['Bacteria', 'Eukaryota']:
        plt.ylim(0, 128)
    else:
        plt.ylim(0, 122)
    plt.xlim(-0.5, 19.5)
    
    for i, samp_size in enumerate(sample_sizes):
        plt.text(i, 102, 'n = ' + str(samp_size), fontname='Arial', fontsize=10, rotation=90, ha='center')

    leg_items = [Patch(facecolor=colors[i], label=ordered_residues[i]) for i in range(len(ordered_residues))]
    leg = plt.legend(handles=leg_items, prop={'family':'Arial', 'size':10}, loc=2, bbox_to_anchor=(1.0, 1.02), title='Secondary AA', handletextpad=0.2)
    leg.set_title('Secondary AA')
    plt.setp(leg.get_title(), multialignment='center', fontname='Arial', fontsize=14)
        
    fig = plt.gcf()
    fig.set_size_inches(7, 5)
    plt.savefig('Fig12A_' + domain + '_Mean_SecondaryLCD_Shares.tif', bbox_inches='tight', dpi=600)
    plt.close()

        
def get_data(domain_of_interest, domain_proteomes):

    lcd_counts_df = {proteome:{} for proteome in domain_proteomes}
    for res1 in amino_acids:
        for res2 in amino_acids.replace(res1, ''):

            for proteome in domain_proteomes:
                lcd_counts_df[proteome][res1+res2] = 0
                
            h = open(os.path.join('SecondaryLCDs_by_LCDcategory', res1+res2 + '_LCDs_AllOrganisms.tsv'))
            header = h.readline()
            for line in h:
                proteome, *remainder = line.rstrip().split('\t')

                if proteome in domain_proteomes:
                    lcd_counts_df[proteome][res1+res2] += 1
                    
            h.close()

    return lcd_counts_df
    
    
def get_proteomes():
    
    h = open('TableS1_LCDfrequency_data_NumberOfProtsWithLCDs.tsv')
    header = h.readline()
    df = {domain:set() for domain in domains}
    for line in h:
        proteome, taxid, domain, *remainder = line.rstrip().split('\t')
        if domain not in domains:
            continue
        df[domain].add(proteome)
        
    h.close()
    
    return df
    
            
if __name__ == '__main__':
    main()
