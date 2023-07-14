

import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import seaborn as sns
import pandas as pd
import statistics
import os
from Bio import SeqIO
amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
ordered_aas = 'ILVFWYMCQNPGASTDEKRH'

domains = ['Archaea','Bacteria','Eukaryota','Viruses']
domain_to_filenum = {'Archaea':'Fig6A', 'Bacteria':'Fig6B', 'Eukaryota':'Fig6C', 'Viruses':'Fig6D'}

def main():

    category = 'PerResidueLCDoccupancy'
    top_x_threshold = 3
    
    prim_aa_df = {aa:[] for aa in ordered_aas}
    barplot_df = {'Domain':[],
                'LCD Class':[],
                'Value':[]}
                
    output = open('ProteomeSize_Statistics.tsv', 'w')
    output.write('\t'.join(['Domain', 'Number of Organisms', 'Average Number of Proteins per Proteome', 'Standard Deviation in Number of Proteins per Proteome']) + '\n')

    for domain in domains:
    
        df = {'LCD Class':[],
                'Value':[],
                'Organism Rank':[],
                'Organism Label':[],
                'Total Proteins':[]}
        h = open('TableS9_PerResidueLevelStatistics_Top3values_for_LCDclasses.tsv')
        header = h.readline()
        mean_prots, stdev_prots, total_orgs = get_mean_num_proteins(domain)
        output.write('\t'.join([str(x) for x in (domain, total_orgs, mean_prots, stdev_prots)]) + '\n')
        ordered_numprots = [''] * ((top_x_threshold+1) * 20)

        lcd_class_index = 0
        for line in h:
            items = line.rstrip().split('\t')
            dom = items[1]
            val_cat = items[3]
            if dom != domain or val_cat != 'Percentage of Proteome in LCDs':
                continue
            lcd_class = items[0]
            if len(lcd_class) > 1:
                continue
                
            perc_vals = items[2].split('; ')
            perc_vals = [float(x) for x in perc_vals]
            organisms = items[5].split('; ')
            mean_perc = float(items[-1])
            total_prots = items[7].split('; ')
            total_prots = [int(x) for x in total_prots]

            shifts = [0, 20, 40]
            for i in range(top_x_threshold):
                df['LCD Class'].append(lcd_class)
                df['Value'].append(perc_vals[i])
                df['Organism Rank'].append(i)
                df['Organism Label'].append(organisms[i])
                df['Total Proteins'].append(total_prots[i])
                ordered_numprots[lcd_class_index+shifts[i]] = ', ' + str(total_prots[i])
                
            df['LCD Class'].append(lcd_class)
            df['Value'].append(mean_perc)
            df['Organism Rank'].append(top_x_threshold+1)
            df['Organism Label'].append('Average Among ' + domain)
            if lcd_class == 'A':
                df['Total Proteins'].append(', ' + str(mean_prots) + ', ' + str(stdev_prots))
                ordered_numprots[60] = ', ' + str(round(mean_prots)) + ', ' + str(round(stdev_prots))
            else:
                df['Total Proteins'].append('')
                
            lcd_class_index += 1

        h.close()

        barplot(df, domain, top_x_threshold, ordered_numprots)
        
    output.close()
    
    
def get_mean_num_proteins(domain):
    
    prots_per_organism = []
    for file in os.listdir('./' + domain):
        if not file.endswith('.fasta') or '_DNA' in file or '_additional' in file:
            continue

        h = open('./' + domain + '/' + file)
        num_prots = 0
        for seq_record in SeqIO.parse(h, 'fasta'):
            seq = str(seq_record.seq)
            if seq[-1] == '*':
                seq = seq[:-1]
            
            if len(seq) < 20:
                continue
            
            num_prots += 1
            
        prots_per_organism.append(num_prots)
            
        h.close()
        
    mean_prots = statistics.mean(prots_per_organism)
    stdev_prots = statistics.stdev(prots_per_organism, mean_prots)
    total_orgs = len(prots_per_organism)
    
    return mean_prots, stdev_prots, total_orgs
        
    
def barplot(df, domain, top_x_threshold, ordered_numprots):
    
    colors = sns.color_palette()
    
    ax = sns.barplot(x='LCD Class', y='Value', data=df, hue='Organism Rank', palette=colors)

    heights = []
    for i, p in enumerate(ax.patches):
        height = p.get_height()
        heights.append(height)
        numprots = ordered_numprots[i]
        ax.annotate(str(round(height, 2)) + numprots, (p.get_x() + p.get_width() / 2., height),
                ha='center', va='bottom', xytext =(0, 2), textcoords='offset points', rotation=90, size=6)

    labels = ['1st', '2nd', '3rd']
    labels.append('Average Among\n' + domain)
    legend_elements = [Patch(facecolor=colors[i], label=labels[i]) for i in range(len(labels))]

    ax.legend(handles=legend_elements, loc=1, title='Organism Rank', handletextpad=0.3, prop={'family':'Arial', 'size':10})

    max_height = max(heights)
    if domain == 'Eukaryota':
        plt.ylim(0, max_height+max_height*.25)
    else:
        plt.ylim(0, max_height+max_height*.2)
    
    plt.xticks(fontname='Arial', fontsize=12)
    plt.yticks(fontname='Arial', fontsize=12)
    plt.xlabel('LCD Class', fontname='Arial', fontsize=14)
    plt.ylabel('Whole-proteome % of\nResidues within LCDs', fontname='Arial', fontsize=14)
    fig = plt.gcf()
    fig.set_size_inches(14, 4)
    plt.savefig(domain_to_filenum[domain] + '_' + domain + '_Top3organisms_PrimaryLCDs_PerResidueOccupancy_Barplot.tif', bbox_inches='tight', dpi=600, pil_kwargs={'compression':'tiff_lzw'})
    plt.close()

        
if __name__ == '__main__':
    main()