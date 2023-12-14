
import os
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np

domains = ['Archaea', 'Bacteria', 'Eukaryota', 'Viruses']
prim_thresholds = [x for x in range(30, 70, 10)]
sec_thresholds = [x for x in range(20, 70, 10)]
amino_acids = 'ACDEFGHIKLMNPQRSTVWY'

def main():

    all_matrices = {}

    for domain in domains:
        filename = domain + '_VariedCompThresholds_PerResidueLCDoccupancies.tsv'
        matrix = pd.read_csv(filename, sep='\t')
        matrix.set_index('Secondary LCD Classes', inplace=True)
        df = pd.DataFrame.to_dict(matrix)
        matrix = matrix.to_numpy()
        
        all_matrices[domain] = matrix

    # PLOT CATEGORICAL "WINNER"======
    archaea_matrix = all_matrices['Archaea']
    cat_matrix = []
    for i in range(len(archaea_matrix)):
        row = []
        for j in range(len(archaea_matrix[i])):
            vals = [all_matrices[domain][i][j] for domain in domains]
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
    plot_categorical_matrix(cat_matrix, 'PerResidueLCDoccupancies')


def plot_categorical_matrix(matrix, data_type):
    
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#17becf', '#7f7f7f', '0.9']
        
    matrix = np.array([np.array(row) for row in matrix])

    mask = [[1 if matrix[i][j] == -1000 else 0 for j in range(len(matrix[i]))] for i in range(len(matrix))]
    mask = np.array([np.array(row) for row in mask])

    labels = domains + ['Tied', 'Zero', 'Invalid\nthreshold']
    
    sns.heatmap(matrix, cmap=colors)
    ax = plt.gca()
    colorbar = ax.collections[0].colorbar
    colorbar.set_ticks(np.linspace(0.4, 5.6, num=7))
    colorbar.set_ticklabels(labels, fontname='Arial', fontsize=16)
    
    ax.set_facecolor('grey')
    plt.xticks([x+2 for x in range(0, 80, 4)], labels=list(amino_acids), rotation=0, fontname='Arial', fontsize=16)
    plt.yticks([x+2.5 for x in range(0, 100, 5)], labels=list(amino_acids), rotation=0, fontname='Arial', fontsize=16)
    plt.xlabel('Amino Acid $\it{A}$', fontname='Arial', fontsize=16)
    plt.ylabel('Amino Acid $\it{B}$', fontname='Arial', fontsize=16)
    for i in range(4, 100, 4):
        plt.plot((i, i), (0, 100), color='white')
    for i in range(5, 100, 5):
        plt.plot((0, 80), (i, i), color='white')
    fig = plt.gcf()
    fig.set_size_inches(14, 10)
    plt.savefig('FigS20B_VariedCompThreshold_' + data_type + '_CategoricalWinner_PerResidueLCDoccupancies.tif', bbox_inches='tight', dpi=600, pil_kwargs={'compression':'tiff_lzw'})
    plt.close()


if __name__ == '__main__':
    main()