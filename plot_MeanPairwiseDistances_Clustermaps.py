
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

domains = ['Archaea', 'Bacteria', 'Eukaryota']

def main():

    method = 'Manhattan'

    matrix, domains_of_life = get_matrix('Clades', method)
    plot_clustermap_clades(matrix, domains_of_life, method)
    
    matrix, domains_of_life = get_matrix('Domains', method)
    plot_clustermap_domains(matrix, domains_of_life, method)
    
    
def plot_clustermap_domains(matrix, domains_of_life, method):

    cg = sns.clustermap(matrix, xticklabels=True, yticklabels=True)
    cg.ax_heatmap.set_xticklabels(cg.ax_heatmap.get_xmajorticklabels(), fontsize=14, fontname='Arial', rotation=0)
    cg.ax_heatmap.set_yticklabels(cg.ax_heatmap.get_ymajorticklabels(), fontsize=14, fontname='Arial', rotation=0)
    
    ax = cg.ax_heatmap
    ax.set_ylabel('')

    fig = plt.gcf()
    fig.set_size_inches(4, 3)
    
    plt.savefig('FigS22A_Domains_PairwiseMeanDistance_' + method + '_Clustermap.tif', bbox_inches='tight', dpi=600)
    plt.close()
    
    
def plot_clustermap_clades(matrix, domains_of_life, method):
    
    colors = sns.color_palette()
    colors = colors[:4]
    colormap = {domain:colors[i] for i, domain in enumerate(domains)}
    row_colors = [colormap[domain] for domain in domains_of_life]
    
    cg = sns.clustermap(matrix, xticklabels=True, yticklabels=True, row_colors=row_colors, method='complete', metric='cityblock')
    cg.ax_heatmap.set_xticklabels(cg.ax_heatmap.get_xmajorticklabels(), fontsize=10, fontname='Arial')
    cg.ax_heatmap.set_yticklabels(cg.ax_heatmap.get_ymajorticklabels(), fontsize=10, fontname='Arial')
    
    ax = cg.ax_heatmap
    ax.set_ylabel('')

    fig = plt.gcf()
    fig.set_size_inches(10, 9)
    
    plt.savefig('FigS22B_Clades_PairwiseMeanDistance_' + method + '_Clustermap.tif', bbox_inches='tight', dpi=600)
    plt.close()


def get_matrix(data_type, method):
    
    h = open(data_type + '_MeanPairwise' + method + 'Distance_FULL_CALCULATIONS.tsv')
    header = h.readline().rstrip().split('\t')
    header = header[2:]
    if data_type == 'Clades':
        header = h.readline().rstrip().split('\t')
        header = header[2:]
        
    matrix = []
    categories = []
    domains_of_life = []    # THIS IS ONLY USED FOR THE "Clade" DATA TYPE. FOR "Domains", THE DOMAINS SERVE AS THE CATEGORIES
    
    for line in h:
        items = line.rstrip().split('\t')
        if data_type == 'Clades':
            domain, cat = items[:2]
            data = [cat] + [float(x) for x in items[2:]]
        else:
            domain = ''
            cat = items[0]
            data = [cat] + [float(x) for x in items[1:]]
        categories.append(cat)
        domains_of_life.append(domain)
        matrix.append(data)

    matrix = pd.DataFrame(matrix, columns=['CATEGORY']+header)
    matrix.set_index('CATEGORY', inplace=True)

    return matrix, domains_of_life
    

if __name__ == '__main__':
    main()
    
