
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import math
from matplotlib.patches import Patch
import textwrap as tw

amino_acids = 'ACDEFGHIKLMNPQRSTVWY'

def main():

    plotting_category = 'SecondaryLCDs'
    num_terms_plotting = 25
    domain_of_interest = 'Eukaryota'
    fig_labels = ['FigS16A', 'FigS16B', 'FigS16C']

    master_df = {'Domain':[],
                'Percentage of Organisms\nwith Significant Enrichment':[]}

    for i, lcd_class_of_interest in enumerate(['HL', 'NH', 'NM']):
        fig_label = fig_labels[i]

        h = open('TableS6_GOterm_Summary_AllOrganisms_AllClasses.tsv')
            
        goterms = []
        perc_orgs = []
        lcd_classes = []
        go_categories = []

        header = h.readline()
        
        for line in h:
            items = line.rstrip().split('\t')
            domain = items[0]
            if domain != domain_of_interest:
                continue
                
            lcd_class = items[1]
            if len(lcd_class) == 1 and plotting_category == 'SecondaryLCDs':
                continue
            elif len(lcd_class) == 2 and plotting_category == 'PrimaryLCDs':
                continue
            elif lcd_class[:len(lcd_class_of_interest)] != lcd_class_of_interest:
                continue
                
            perc_enr = float(items[7])
            
            goterms.append(items[4])
            go_categories.append(items[3])
            perc_orgs.append(perc_enr)
            lcd_classes.append(lcd_class)
            
        h.close()
    
        if len(lcd_classes) == 0:
            continue
        perc_orgs, goterms, lcd_classes, go_categories = zip(*sorted(zip(perc_orgs, goterms, lcd_classes, go_categories), reverse=True))

        master_df['Domain'] += [domain_of_interest for x in range(len(perc_orgs))]
        master_df['Percentage of Organisms\nwith Significant Enrichment'] += perc_orgs
        goterms = np.array(goterms)
        perc_orgs = np.array(perc_orgs)
        
        plotting(domain_of_interest, goterms, go_categories, perc_orgs, lcd_classes, plotting_category, num_terms_plotting, lcd_class_of_interest, fig_label)


def plotting(domain, goterms, go_categories, perc_orgs, lcd_classes, plotting_category, num_terms_plotting, lcd_class_of_interest, fig_label):
    
    tick_colors = ['#be0119', '#00035b', '#0b5509'] # SCARLET, DARK BLUE, FOREST
    cats = ['BP', 'MF', 'CC']
    tick_cmap = {cat:tick_colors[i] for i, cat in enumerate(cats)}
    
    colors = list(sns.color_palette(n_colors=50))
    
    xtick_labels = ['\n'.join(tw.wrap(x, width=35)) for x in goterms[:num_terms_plotting]]
    perc_orgs = perc_orgs[:num_terms_plotting]
    lcd_classes = lcd_classes[:num_terms_plotting]

    encountered_lcd_classes = set()
    nonredundant_lcd_classes = []
    for lcd_class in lcd_classes:
        if lcd_class not in encountered_lcd_classes:
            reverse_class = lcd_class[1] + lcd_class[0]
            if reverse_class in nonredundant_lcd_classes:
                loc = nonredundant_lcd_classes.index(reverse_class)
                new_class = reverse_class + '/' + lcd_class
                nonredundant_lcd_classes[loc] = new_class
            else:
                nonredundant_lcd_classes.append(lcd_class)
        encountered_lcd_classes.add(lcd_class)

    color_map = {}
    color_map_nonredundant = {}
    for i, lcd_class in enumerate(nonredundant_lcd_classes):
        color_map_nonredundant[lcd_class] = colors[i]
        classes = lcd_class.split('/')
        for c in classes:
            color_map[c] = colors[i]

    df = {'GO term':xtick_labels,
            'Percentage of Organisms\nwith Significant Enrichment':perc_orgs,
            'Color':[color_map[lcd_class] for lcd_class in lcd_classes]}
    plt.bar([x for x in range(len(xtick_labels))], perc_orgs[:num_terms_plotting], color=[color_map[lcd_class] for lcd_class in lcd_classes])

    legend_elements = [Patch(facecolor=color_map_nonredundant[lcd_class], label=lcd_class) for lcd_class in color_map_nonredundant]
                         
    plt.legend(handles=legend_elements, handletextpad=0.3, prop={'family':'Arial', 'size':14})
    
    plt.xticks([x for x in range(len(xtick_labels))], labels=xtick_labels, rotation=90, va='top', fontname='Arial', fontsize=16)
    plt.yticks(fontname='Arial', fontsize=18)
    
    ax = plt.gca()
    xticklabels = ax.get_xticklabels()
    for i, label in enumerate(go_categories[:num_terms_plotting]):
        xticklabels[i].set_color(tick_cmap[label])
    
    plt.xlabel('GO term', fontname='Arial', fontsize=18)
    plt.ylabel('Percentage of Organisms\nwith Significant Enrichment', fontname='Arial', fontsize=18)

    fig = plt.gcf()

    fig.set_size_inches(20, 5)
    plt.xlim(-1, num_terms_plotting)
    plt.savefig(fig_label + '_' + lcd_class_of_interest + '-LCDs_' + domain + '_Top' + str(num_terms_plotting) + '_GOterms.tif', bbox_inches='tight', dpi=600)
    plt.close()
    

if __name__ == '__main__':
    main()
    
