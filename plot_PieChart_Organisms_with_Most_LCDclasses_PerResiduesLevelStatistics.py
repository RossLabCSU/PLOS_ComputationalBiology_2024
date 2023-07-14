
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Patch
domains = ['Archaea', 'Bacteria', 'Eukaryota', 'Viruses']

amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
aa_strings = [aa for aa in amino_acids]
for res1 in amino_acids:
    for res2 in amino_acids:
        if res1 == res2:
            continue
        aa_strings.append(res1+res2)
        
excluded_orgs = ['Sepia pharaonis', 'Spodoptera litura']
        
def main():

    for domain_of_interest in domains:
        file = 'TableS9_PerResidueLevelStatistics_Top3values_for_LCDclasses.tsv'

        h = open(file)
        header = h.readline()
        
        organisms = []
        for line in h:
            lcd_class, domain, main_vals, val_cat, proteomes, sci_names, common_names, total_prots, total_aas, mean_perc_in_lcd = line.rstrip().split('\t')
                
            if domain != domain_of_interest:
                continue
                
            # SKIPS LCD CLASSES ABSENT FROM THE GIVEN DOMAIN OF LIFE
            main_vals = main_vals.split('; ')
            if main_vals[0] == '0.0':
                continue
                
            sci_names = sci_names.split('; ')
            sci_names = [x for x in sci_names if x not in excluded_orgs]
            sci_name = sci_names[0]

            organisms.append(sci_name)
            
        h.close()
        
        freqs, annots = calc_freqs(organisms)
        
        print('\n', domain_of_interest, freqs[:5], annots[:5])
        print(domain_of_interest, 'Summed percent share of top5:', sum(freqs[:5]) / 400 *100)
            
        piechart(freqs, annots, domain_of_interest)

    
def calc_freqs(organisms):

    annots = []
    freqs = []
    for organism in set(organisms):
        freq = organisms.count(organism)
        annots.append(organism)
        freqs.append(freq)

    freqs, annots = zip(*sorted(zip(freqs, annots), reverse=True))

    return freqs, annots
    
    
def piechart(freqs, pfams, domain):

    # DEFAULT COLOR PALETTE
    colors = sns.color_palette("Set2")
    colors += sns.color_palette('pastel')

    total = sum(freqs)
    reverse = freqs[::-1]
    smallest_val_of_interest = min(freqs[:10])
    rev_index = reverse.index(smallest_val_of_interest)
    end_of_small_val = len(freqs) - rev_index-1
    start_of_small_val = freqs.index(smallest_val_of_interest)
    other = sum(freqs[end_of_small_val+1:])
    labels = list(pfams[:end_of_small_val+1]) + ['Other']
    sizes = list(freqs)
    
    wedges = plt.pie(sizes, startangle=20, textprops={'fontfamily':'Arial', 'fontsize':12}, colors=colors, counterclock=False)

    ax = plt.gca()
    ax.axis('equal')  # EQUAL ASPECT RATIO ENSURES THAT PIE IS DRAWN AS A CIRCLE

    leg_cutoff = 5  # MAXIMUM NUMBER OF ORGANISMS TO INCLUDE IN THE LEGEND

    for i, label in enumerate(labels):
        if label == 'candidate division MSBL1 archaeon SCGC-AAA385M02':
            labels[i] = 'candidate division MSBL1\narchaeon SCGC-AAA385M02'
        elif label == 'Methanobrevibacter arboriphilus JCM 13429 = DSM 1125':
            labels[i] = 'Methanobrevibacter arboriphilus'
        elif label == 'White spot syndrome virus (isolate Shrimp/China/Tongan/1996)':
            labels[i] = 'White spot syndrome virus'
            
    legend_elements = [Patch(facecolor=colors[i], label=labels[i]) for i in range(leg_cutoff)]

    legend = ax.legend(handles=legend_elements, bbox_to_anchor=(0.86, 0, 0.5, 1), loc='center left', handletextpad=0.3, title='Organism', prop={'family':'Arial', 'size':12})

    leaders = {'Archaea':'FigS9A', 'Bacteria':'FigS9C', 'Eukaryota':'Fig8A', 'Viruses':'FigS9E'}
    plt.setp(legend.get_texts(), fontsize=8, va='bottom')
    fig = plt.gcf()
    fig.set_size_inches(6, 4)
    plt.savefig(leaders[domain] + '_' + domain + '_PieChart_NumberPerResidueOccupancyCategories.tif', bbox_inches='tight', dpi=600, pil_kwargs={'compression':'tiff_lzw'})
    plt.close()

    
if __name__ == '__main__':
    main()