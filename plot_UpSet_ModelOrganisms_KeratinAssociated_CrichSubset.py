
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "Arial"
plt.rcParams["font.size"] = 14
import seaborn as sns
import upsetplot as usp

domains = ['Archaea','Bacteria','Eukaryota','Viruses']
organism_to_filenum = {'Human':'Fig3B', 'Bovine':'Fig3C', 'Dog':'Fig3D', 'Mouse':'Fig3E', 'Rat':'Fig3F', 'Pig':'Fig3G'}

def main():

    h = open('TableS7_ModelOrganisms_KeratinAssociated_CrichSubset.tsv')
    header = h.readline()
    all_classes = {}

    df = {}
    for line in h:
        proteome, sciname, common_name, lcd_class, prots = line.rstrip().split('\t')
        prots = prots.split(', ')
        df[common_name] = df.get(common_name, {})
        df[common_name][lcd_class] = set(prots)
        all_classes[common_name] = all_classes.get(common_name, [])
        all_classes[common_name].append(lcd_class)
        
    h.close()
    
    df = combine_like_classes(df)

    for organism in df:
        prot_sets = df[organism]
        if len(list(prot_sets.keys())) < 2:
            continue
        plotting_df = usp.from_contents(prot_sets)
        usp.plot(plotting_df, show_counts=True, sort_by='cardinality', facecolor='#87ae73')
        plt.ylabel('Intersection Size', fontname='Arial', fontsize=16)
        if organism == 'Human' or organism == 'Rat':
            plt.yticks([0, 2, 4, 6, 8, 10])

        plt.title(organism)
        plt.savefig(organism_to_filenum[organism] + '_' + organism + '_KeratinAssociatedProteins_UpSetPlot.tif', bbox_inches='tight', dpi=600, pil_kwargs={'compression':'tiff_lzw'})
        plt.close()


def combine_like_classes(df):

    new_df = {organism:{} for organism in df}
    for organism in df:
        classes = sorted(list(df[organism].keys()))
        for c in classes:
            recip = c[1] + c[0]
            if recip in classes:
                combined_prots = df[organism][c] | df[organism][recip]
                if c[0] == 'C':
                    new_str = c + '/' + recip
                    new_df[organism][new_str] = combined_prots
            else:
                new_df[organism][c] = df[organism][c]
                
    return new_df
    
    
if __name__ == '__main__':
    main()