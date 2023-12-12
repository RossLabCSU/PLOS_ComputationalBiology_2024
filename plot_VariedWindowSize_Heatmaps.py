
import os
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np

domains = ['Archaea', 'Bacteria', 'Eukaryota', 'Viruses']
window_sizes = [x for x in range(20, 70, 10)]
amino_acids = 'ACDEFGHIKLMNPQRSTVWY'


def main():

    # PLOT ACTUAL VALUES======
    filename = 'AllDomains_VariedWindowSize_OrganismLevel_LCDfrequencies.tsv'
    matrix = pd.read_csv(filename, sep='\t')
    matrix.set_index('Secondary LCD Classes', inplace=True)
    df = pd.DataFrame.to_dict(matrix)
    matrix = matrix.to_numpy()
    
    plot_matrix(matrix, 'OrganismLevel_LCDfrequencies')
    #=========================
    
    # PLOT CATEGORICAL "WINNER"======
    cat_matrix = []
    for sec_aa in amino_acids:
        row = []
        for prim_aa in amino_acids:
            for window_size in window_sizes:
                key1 = prim_aa + '_' + str(window_size) + 'aaWindowSize'
                vals = [df[key1][domain+'_'+sec_aa] for domain in domains]

                max_val = max(vals)
                if max_val == 0:
                    row.append(5)
                elif vals.count(max_val) > 1 and max_val > -1000:
                    row.append(4)
                elif vals.count(max_val) > 1 and max_val == -1000:
                    row.append(6)
                else:
                    max_index = vals.index(max(vals))
                    row.append(max_index)
        cat_matrix.append(row)
    
    plot_categorical_matrix(cat_matrix, 'OrganismLevel_LCDfrequencies')
    

def plot_matrix(matrix, data_type):

    mask = np.array([np.array([1 if x == -1000 else 0 for x in matrix[i]]) for i in range(len(matrix))])
    
    matrix = np.array([np.array(row) for row in matrix])
    
    sns.heatmap(matrix, mask=mask)
    ax = plt.gca()
    ax.set_facecolor('0.8')

    colorbar = ax.collections[0].colorbar
    colorbar.set_ticklabels([round(x) for x in colorbar.get_ticks()], fontname='Arial', fontsize=16)
    
    plt.xticks([x+2.5 for x in range(0, 100, 5)], labels=list(amino_acids), rotation=0, fontname='Arial', fontsize=16)
    plt.yticks([x+2 for x in range(0, 80, 4)], labels=list(amino_acids), rotation=0, fontname='Arial', fontsize=16)
    plt.xlabel('Primary Amino Acid', fontname='Arial', fontsize=16)
    plt.ylabel('Secondary Amino Acid', fontname='Arial', fontsize=16)
    for i in range(5, 100, 5):
        plt.plot((i, i), (0, 80), color='grey')
    for i in range(4, 80, 4):
        plt.plot((0, 100), (i, i), color='grey')
    fig = plt.gcf()
    fig.set_size_inches(14, 8)
    plt.savefig('FigS2A_VariedWindowSize_' + data_type + '.tif', bbox_inches='tight', dpi=600)
    plt.close()
    
    
def plot_categorical_matrix(matrix, data_type):
    
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#17becf', '#7f7f7f']
        
    matrix = np.array([np.array(row) for row in matrix])

    mask = [[1 if matrix[i][j] == -1000 else 0 for j in range(len(matrix[i]))] for i in range(len(matrix))]
    mask = np.array([np.array(row) for row in mask])

    labels = domains + ['Tied', 'Zero']
    
    sns.heatmap(matrix, cmap=colors)
    ax = plt.gca()
    colorbar = ax.collections[0].colorbar
    colorbar.set_ticks(np.linspace(0.4, 4.6, num=6))
    colorbar.set_ticklabels(labels, fontname='Arial', fontsize=16)

    ax.set_facecolor('grey')
    plt.xticks([x+2.5 for x in range(0, 100, 5)], labels=list(amino_acids), rotation=0, fontname='Arial', fontsize=16)
    plt.yticks([x+0.5 for x in range(20)], labels=list(amino_acids), rotation=0, fontname='Arial', fontsize=16)
    plt.xlabel('Primary Amino Acid', fontname='Arial', fontsize=16)
    plt.ylabel('Secondary Amino Acid', fontname='Arial', fontsize=16)
    for i in range(5, 100, 5):
        plt.plot((i, i), (0, 20), color='white')
    for i in range(100):
        if i%5 == 0:
            continue
        plt.plot((i, i), (0, 20), color='grey', linewidth=0.5)
    for i in range(1, 20):
        plt.plot((0, 100), (i, i), color='white')
    fig = plt.gcf()
    fig.set_size_inches(14, 8)
    plt.savefig('FigS2B_VariedWindowSize_' + data_type + '_CategoricalWinner.tif', bbox_inches='tight', dpi=600)
    plt.close()


if __name__ == '__main__':
    main()
