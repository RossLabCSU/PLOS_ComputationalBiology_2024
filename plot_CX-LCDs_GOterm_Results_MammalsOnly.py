
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import math
from matplotlib.patches import Patch
amino_acids = 'ACDEFGHIKLMNPQRSTVWY'

def main():

    plotting_category = 'SecondaryLCDs'
    num_terms_plotting = 50
    lcd_class_of_interest = 'C'
    domain_of_interest = 'Eukaryota'

    mammal_taxids, taxid_to_sciname = get_mammals()
    
    master_df = {'Domain':[],
                'Percentage of Organisms\nwith Significant Enrichment':[]}
                
    output = open('Mammals_' + lcd_class_of_interest + 'richLCDs_from_AllEukaryotesGOresults.tsv', 'w')

    h = open('TableS6_GOterm_Summary_AllOrganisms_AllClasses.tsv')
        
    goterms = []
    perc_orgs = []
    lcd_classes = []
    go_categories = []

    header = h.readline().rstrip().split('\t')
    output.write('\t'.join(header[:-1] + ['Corresponding Organism Scientific Name (semicolon-delimited)', header[-1]]) + '\n')
    
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
        
        sig_organisms = items[9].split(';')
        sig_taxids = []
        for org in sig_organisms:
            leader, taxid = org.split('_')
            sig_taxids.append(taxid)
            
        checks = [1 if taxid in mammal_taxids else 0 for org in sig_taxids]
        if sum(checks) == 0:
            continue
            
        model_org_indices = [i for i, check_val in enumerate(checks) if check_val == 1]
        filtered_perc_orgs = len(model_org_indices) / len(mammal_taxids) * 100
        pvals = items[8].split(';')
        filtered_pvals = [pvals[i] for i in model_org_indices]
        filtered_orgs = [sig_organisms[i] for i in model_org_indices]
        filtered_taxids = [taxid for taxid in sig_taxids if taxid in mammal_taxids]
        filtered_org_scinames = [taxid_to_sciname[taxid] for taxid in filtered_taxids]
        
        perc_prots_in_background = items[10].split(';')
        filtered_perc_prots = [perc_prots_in_background[i] for i in model_org_indices]
        
        output.write('\t'.join(items[:5] + [str(len(filtered_orgs)), str(len(mammal_taxids)), str(filtered_perc_orgs), ';'.join(filtered_pvals), ';'.join(filtered_orgs), ';'.join(filtered_org_scinames), ';'.join(filtered_perc_prots)]) + '\n')

        goterms.append(items[4])
        go_categories.append(items[3])
        perc_orgs.append(filtered_perc_orgs)
        lcd_classes.append(lcd_class)
        
    h.close()

    perc_orgs, goterms, lcd_classes, go_categories = zip(*sorted(zip(perc_orgs, goterms, lcd_classes, go_categories), reverse=True))

    master_df['Domain'] += [domain_of_interest for x in range(len(perc_orgs))]
    master_df['Percentage of Organisms\nwith Significant Enrichment'] += perc_orgs
    goterms = np.array(goterms)
    perc_orgs = np.array(perc_orgs)
    
    plotting(domain_of_interest, goterms, go_categories, perc_orgs, lcd_classes, plotting_category, num_terms_plotting, lcd_class_of_interest)

    
def get_mammals():

    h = open('UniProt_Taxonomy.tsv')
    header = h.readline()
    
    mammal_taxids = []
    taxid_to_sciname = {}
    
    for line in h:
        items = line.rstrip().split('\t')
        lineage = items[7].split(', ')
        if 'Mammalia' in lineage:
            mammal_taxids.append(items[0])
            taxid_to_sciname[items[0]] = items[2]
            
    h.close()
    
    return mammal_taxids, taxid_to_sciname
    

def plotting(domain, goterms, go_categories, perc_orgs, lcd_classes, plotting_category, num_terms_plotting, lcd_class_of_interest):
    
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
    for i, val in enumerate(perc_orgs):
        offset = max(perc_orgs) * 0.01
        plt.text(i, val+offset, lcd_classes[i], ha='center', fontname='Arial', fontsize=14)

    legend_elements = [Patch(facecolor=color_map_nonredundant[lcd_class], label=lcd_class) for lcd_class in color_map_nonredundant]
                         
    plt.legend(handles=legend_elements, handletextpad=0.3, prop={'family':'Arial', 'size':14})
    
    plt.xticks([x for x in range(len(xtick_labels))], labels=xtick_labels, rotation=90, va='top', fontname='Arial', fontsize=14)
    plt.yticks(fontname='Arial', fontsize=14)
    
    ax = plt.gca()
    xticklabels = ax.get_xticklabels()
    for i, label in enumerate(go_categories[:num_terms_plotting]):
        xticklabels[i].set_color(tick_cmap[label])
    
    plt.xlabel('GO term', fontname='Arial', fontsize=16)
    plt.ylabel('Percentage of Organisms\nwith Significant Enrichment', fontname='Arial', fontsize=16)
    fig = plt.gcf()
    fig.set_size_inches(20, 6.5)
    plt.xlim(-1, num_terms_plotting)
    plt.savefig('FigS10B_CX-LCDs_Top' + str(num_terms_plotting) + '_GOterms_MammalsOnly.tif', bbox_inches='tight', dpi=600, pil_kwargs={'compression':'tiff_lzw'})
    plt.close()
    

if __name__ == '__main__':
    main()
    