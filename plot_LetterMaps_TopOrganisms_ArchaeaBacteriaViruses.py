
import seaborn as sns
import matplotlib.pyplot as plt

# NUMBER OF ORGANISMS FOR WHICH TO PLOT THE LETTER PLOTS
# THESE VALUES WERE CHOSEN BASED ON THE PIE CHARTS OF NUMBER OF CLASSES FOR WHICH EACH ORGANISM HAD THE HIGHEST PERCENT IN LCDs VALUE
num_of_orgs_df = {'Archaea':5,
                'Bacteria':5,
                'Viruses':5}
# ====================================================== #
                
amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
ordered_aas = 'ILVFWYMCQNPGASTDEKRH'

domains = ['Archaea','Bacteria','Viruses']
domain_to_filenum = {'Archaea':'FigS9B', 'Bacteria':'FigS9D', 'Viruses':'FigS9F'}
aa_strings = []
for res1 in amino_acids:
    for res2 in amino_acids:
        if res1 == res2:
            continue
        aa_strings.append(res1+res2)
        
        
def main():

    category_of_interest = 'Percentage of Proteome in LCDs'
    
    colors = sns.color_palette('Set2')

    for domain in domains:
        file = 'TableS9_PerResidueLevelStatistics_Top3values_for_LCDclasses.tsv'
        top_x_orgs = num_of_orgs_df[domain]
        top_x_counts, top_x_scinames, top_x_commonnames = get_frequent_orgs(file, top_x_orgs, domain, category_of_interest)

        index = 0
        plot_index = 0
        for color_index, sci_name in enumerate(top_x_scinames):
            df, prim_classes = get_data(file, sci_name, category_of_interest)
            pos1_vals = []
            pos2_vals = []
            for i, res1 in enumerate(ordered_aas):
                for j, res2 in enumerate(ordered_aas):
                    val = df[res1][res2]
                    if val == 1:
                        pos1_vals.append(ordered_aas.index(res1)+1) # THESE VALUES ARE SHIFTED BY 1 SINCE THE LETTERS ARE PLOTTED STARTING AT X-VAL==1. MANUALLY VERIFIED THAT THE PLOTTED LCD CLASSES MATCH THE DATA FILE.
                        pos2_vals.append(ordered_aas.index(res2)+1) # THESE VALUES ARE SHIFTED BY 1 SINCE THE LETTERS ARE PLOTTED STARTING AT X-VAL==1. MANUALLY VERIFIED THAT THE PLOTTED LCD CLASSES MATCH THE DATA FILE.

            plot_xlinks(pos1_vals, pos2_vals, sci_name, colors[color_index], prim_classes, domain, index)
            if (index == 1 and domain in ['Archaea', 'Bacteria']) or (index==4 and domain in ['Archaea', 'Bacteria', 'Viruses']):
                plt.savefig(domain_to_filenum[domain] + '_' + domain + str(plot_index) + '_SharedLetterMap.tif', bbox_inches='tight', dpi=600)
                plt.close()
                plot_index += 1
            index += 1
        
        plt.close()
        
        
def get_data(file, sciname_of_interest, category_of_interest):

    h = open(file)
    header = h.readline()
    prim_classes = []
    df = {aa:{res:0 for res in amino_acids} for aa in amino_acids}
    for line in h:
        lcd_class, domain, max_vals, val_cat, proteomes, sci_names, common_names, *remainder = line.rstrip().split('\t')
        
        sci_names = sci_names.split('; ')
        common_names = common_names.split('; ')
        sci_name = sci_names[0]
        common_name = common_names[0]
            
        if sci_name != sciname_of_interest:
            continue
            
        max_vals = max_vals.split('; ')
        if max_vals[0] == '0.0':
            continue
            
        if len(lcd_class) == 1:
            res1 = lcd_class[:]
            res2 = lcd_class[:]
            df[res1][res2] = 1
        else:
            res1, res2 = lcd_class
            df[res1][res2] = 1
            
    h.close()
    
    return df, prim_classes
    

def get_frequent_orgs(file, top_x_orgs, domain_of_interest, category_of_interest):

    h = open(file)
    header = h.readline()
    sci_names = []
    common_names = []
    scinames_to_commonnames = {}
    for line in h:
        lcd_class, domain, max_vals, val_cat, proteomes, sci_names_top3, common_names_top3, *remainder = line.rstrip().split('\t')
        if domain != domain_of_interest:
            continue
            
        max_vals = max_vals.split('; ')
        if max_vals[0] == '0.0':
            continue

        sci_names_top3 = sci_names_top3.split('; ')
        common_names_top3 = common_names_top3.split('; ')
        sci_name = sci_names_top3[0]
        common_name = common_names_top3[0]

        sci_names.append(sci_name)
        common_names.append(common_name)
        scinames_to_commonnames[sci_name] = common_name
        
    h.close()
    
    org_counts = []
    sci_names_set = []
    common_names_set = []
    for sci_name in set(sci_names):
        count = sci_names.count(sci_name)
        org_counts.append(count)
        sci_names_set.append(sci_name)
        common_names_set.append(scinames_to_commonnames[sci_name])
        
    org_counts, sci_names_set, common_names_set = zip(*sorted(zip(org_counts, sci_names_set, common_names_set), reverse=True))
    
    top_x_counts = org_counts[:top_x_orgs]
    top_x_scinames = sci_names_set[:top_x_orgs]
    top_x_commonnames = common_names_set[:top_x_orgs]
    
    return top_x_counts, top_x_scinames, top_x_commonnames


def plot_xlinks(pos1_vals, pos2_vals, dataset, line_color, prim_aas, domain, index):

    domain_ypos2 = 10
    domain_height = 2
    
    for i in range(len(pos1_vals)):
        positions = (pos1_vals[i], pos2_vals[i])
        plt.plot( [pos1_vals[i], pos2_vals[i]], [domain_ypos2, 1+domain_height], color=line_color, alpha=0.5, zorder=0)

    for i, letter in enumerate(ordered_aas):
        plt.text(i+1, 1+0.7, letter, ha='center', fontname='Arial', fontsize=12)
        
    for i, letter in enumerate(ordered_aas):
        if letter in prim_aas:
            plt.text(i+1, domain_ypos2+0.25, letter, ha='center', fontname='Arial', fontsize=12, color='#d62728')
        else:
            plt.text(i+1, domain_ypos2+0.25, letter, ha='center', fontname='Arial', fontsize=12)

    fig = plt.gcf()
    fig.set_size_inches(10, domain_ypos2 / 4.75)
    plt.ylim(0.001, domain_ypos2+domain_height+0.1)
    plt.xlim(0.001, 21)

    ax = plt.gca()
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    plt.yticks([])
    plt.xticks([])
    sci_name, *junk = dataset.split(' (')
    fig.set_size_inches(5.5, 1.75)


if __name__ == '__main__':
    main()
