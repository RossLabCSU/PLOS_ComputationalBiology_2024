
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Patch

def main():

    proteome_of_interest = 'UP000005640_9606'
    organism = 'Hsapiens'

    min_thresholds = [0, 15]
    max_thresholds = [15, 1000000]
    plot_labels = ['FewerThan15', '15orMore']
    
    for i, min_threshold in enumerate(min_thresholds):
        max_threshold = max_thresholds[i]
        plot_label = plot_labels[i]
        h = open('All_LCDclasses_PfamClanMapping.tsv')
        header = h.readline()
        
        df = {'LCD Class':[],
            'Max Count':[],
            'Max Clan':[],
            'Number of LCD Prots':[]}
        
        for line in h:
            domain, proteome, lcd_class, clan_ids, clan_set, perc_vals, max_perc, clan_with_max, num_with_max_clan, num_lcd_prots = line.rstrip().split('\t')
            if proteome != proteome_of_interest:
                continue
                
            if num_with_max_clan == 'N/A' or not min_threshold <= int(num_lcd_prots) < max_threshold:
                continue
                
            clan_with_max = clan_with_max.split(';')
            clan_with_max = clan_with_max[:2]
                
            df['LCD Class'].append(lcd_class)
            df['Max Count'].append(int(num_with_max_clan))
            df['Max Clan'].append(', '.join(clan_with_max))
            df['Number of LCD Prots'].append(int(num_lcd_prots))
            
        h.close()

        plotting(df, organism, plot_label)
    
    
def plotting(df, organism, plot_label):
    
    colors = sns.color_palette()
    
    max_counts = df['Max Count']
    lcd_classes = df['LCD Class']
    max_clans = df['Max Clan']
    num_prots = df['Number of LCD Prots']
    
    num_prots, max_counts, lcd_classes, max_clans = zip(*sorted(zip(num_prots, max_counts, lcd_classes, max_clans), reverse=True))

    plt.bar(lcd_classes, num_prots)
    plt.bar(lcd_classes, max_counts, color=colors[3])

    legend_elements = [Patch(facecolor=colors[0], label='Total # of Proteins with LCD'),
                    Patch(facecolor=colors[3], label='Maximum # of Proteins with\nSingle Pfam Clan Annotation')]

    ax = plt.gca()
    ax.legend(handles=legend_elements, loc='upper right', handletextpad=0.3)    

    plt.xticks(rotation=90, fontname='Arial', fontsize=10)
    plt.xlabel('LCD Class', fontname='Arial', fontsize=12)
    plt.ylabel('Number of Proteins', fontname='Arial', fontsize=12)
    plt.xlim(-1, len(max_counts))
    fig = plt.gcf()
    fig.set_size_inches(20, 5)
    plt.savefig('FigS5A_' + organism + '_PfamData_Barplot_' + plot_label + '.tif', bbox_inches='tight', dpi=600, pil_kwargs={'compression': 'tiff_lzw'})
    plt.close()
    

if __name__ == '__main__':
    main()