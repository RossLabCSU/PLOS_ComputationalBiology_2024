
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.patches import Patch

amino_acids = 'ACDEFGHIKLMNPQRSTVWY'

domains = ['Archaea', 'Bacteria', 'Eukaryota', 'Viruses']

def main():

    lcd_classes = ['C', 'C', 'H', 'H', 'H']
    domains_of_interest = ['Eukaryota', 'Eukaryota', 'Archaea', 'Bacteria', 'Eukaryota']
    include_reciprocals = [False, True, True, True, True]
    fig_labels = ['Fig3A', 'FigS10A', 'FigS15B', 'FigS15A', 'Fig5A']

    plotting_category = 'SecondaryLCDs'

    num_terms_plotting = 50
    show_bar_labels = True
    
    master_df = {'Domain':[],
                'Percentage of Organisms\nwith Significant Enrichment':[]}
    
    for i, lcd_class_of_interest in enumerate(lcd_classes):
        domain_of_interest = domains_of_interest[i]
        include_recip = include_reciprocals[i]
        fig_label = fig_labels[i]

        goterms = []
        perc_orgs = []
        lcd_classes = []
        go_categories = []
        
        h = open('GOterm_Summary_AllOrganisms_AllClasses_PercentsOfOrganismsFIXED.tsv')
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
                
            elif include_recip and lcd_class_of_interest not in lcd_class:
                continue
            elif not include_recip and lcd_class[:len(lcd_class_of_interest)] != lcd_class_of_interest:
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
        
        plotting(domain_of_interest, goterms, go_categories, perc_orgs, lcd_classes, plotting_category, num_terms_plotting, lcd_class_of_interest, show_bar_labels, fig_label)


def plotting(domain, goterms, go_categories, perc_orgs, lcd_classes, plotting_category, num_terms_plotting, lcd_class_of_interest, show_bar_labels, fig_label):
    
    tick_colors = ['#be0119', '#00035b', '#0b5509'] # SCARLET, DARK BLUE, FOREST
    cats = ['BP', 'MF', 'CC']
    tick_cmap = {cat:tick_colors[i] for i, cat in enumerate(cats)}
    colors = list(sns.color_palette(n_colors=50))
    
    xtick_labels = goterms[:num_terms_plotting]
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

    if show_bar_labels:
        for i, val in enumerate(perc_orgs):
            offset = max(perc_orgs) * 0.01
            plt.text(i, val+offset, lcd_classes[i], ha='center', fontname='Arial', fontsize=14)

    legend_elements = [Patch(facecolor=color_map_nonredundant[lcd_class], label=lcd_class) for lcd_class in color_map_nonredundant]
                         
    plt.legend(handles=legend_elements, handletextpad=0.3, prop={'family':'Arial', 'size':14})
    
    ax = plt.gca()
    plt.xticks([x for x in range(len(xtick_labels))], labels=xtick_labels, rotation=90, va='top', fontname='Arial', fontsize=14)
    plt.yticks(fontname='Arial', fontsize=14)
    
    xticklabels = ax.get_xticklabels()
    for i, label in enumerate(go_categories[:num_terms_plotting]):
        xticklabels[i].set_color(tick_cmap[label])
        
    plt.xlabel('GO term', fontname='Arial', fontsize=16)
    plt.ylabel('Percentage of Organisms\nwith Significant Enrichment', fontname='Arial', fontsize=16)

    fig = plt.gcf()
    fig.set_size_inches(20, 6.5)
    plt.xlim(-1, num_terms_plotting)
    plt.savefig(fig_label + '_' + lcd_class_of_interest + '-LCDs_' + domain + '_Top' + str(num_terms_plotting) + '_GOterms_' + plotting_category + '_PercentsOfOrganismsFIXED.tif', bbox_inches='tight', dpi=600, pil_kwargs={'compression':'tiff_lzw'})
    plt.close()

        
if __name__ == '__main__':
    main()