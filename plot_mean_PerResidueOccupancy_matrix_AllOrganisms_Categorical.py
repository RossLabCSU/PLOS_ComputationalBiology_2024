
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

domains = ['Archaea', 'Bacteria', 'Eukaryota', 'Viruses']
amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
all_sec_classes = []
for aa1 in amino_acids:
    for aa2 in amino_acids:
        if aa1 == aa2:
            continue
        all_sec_classes.append(aa1+aa2)

def main():

    df = get_data()
    cat_matrix = []
    for i, sec_aa in enumerate(amino_acids):
        row = []
        for j, prim_aa in enumerate(amino_acids):
            sum_vals = [df[domain][prim_aa][sec_aa]['sum'] for domain in domains]
            n_vals = [df[domain][prim_aa][sec_aa]['n'] for domain in domains]
            vals = [sum_vals[k] / n_vals[k] for k in range(len(sum_vals))]
            
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

    plot_categorical_matrix(cat_matrix, 'AllOrganismsMean')
    
    
def plot_categorical_matrix(matrix, data_type):
    
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
        
    matrix = np.array([np.array(row) for row in matrix])

    mask = [[1 if matrix[i][j] == -1000 else 0 for j in range(len(matrix[i]))] for i in range(len(matrix))]
    mask = np.array([np.array(row) for row in mask])

    labels = domains[:]
    
    sns.heatmap(matrix, cmap=colors)
    ax = plt.gca()
    colorbar = ax.collections[0].colorbar
    colorbar.set_ticks(np.linspace(0.4, 2.6, num=4))
    colorbar.set_ticklabels(labels, fontname='Arial', fontsize=16)
    
    ax.set_facecolor('grey')
    plt.xticks([x+0.5 for x in range(20)], labels=list(amino_acids), rotation=0, fontname='Arial', fontsize=16)
    plt.yticks([x+0.5 for x in range(20)], labels=list(amino_acids), rotation=0, fontname='Arial', fontsize=16)
    
    plt.xlabel('Primary Amino Acid', fontname='Arial', fontsize=16)
    plt.ylabel('Secondary Amino Acid', fontname='Arial', fontsize=16)

    for i in range(1, 20):
        plt.plot((i, i), (0, 20), color='white')
    for i in range(1, 20):
        plt.plot((0, 20), (i, i), color='white')
    
    fig = plt.gcf()
    fig.set_size_inches(14, 10)
    plt.savefig('Fig8_VariedCompThreshold_' + data_type + '_CategoricalWinner_PerResidueLCDoccupancies.tif', bbox_inches='tight', dpi=600, pil_kwargs={'compression':'tiff_lzw'})
    plt.close()
        
    
def get_data():

    h = open('LCDfrequency_data_PerResidue_LCDoccupancy.tsv')
    header = h.readline()
    
    df = {domain:{prim_aa:{sec_aa:{'sum':0, 'n':0} for sec_aa in amino_acids} for prim_aa in amino_acids} for domain in domains}

    for line in h:
        proteome, taxid, domain, sciname, commonname, total_prots, *data = line.rstrip().split('\t')
        prim_class_data = [float(x) for x in data[:20]]
        sec_class_data = [float(x) for x in data[20:]]
        for i, aa in enumerate(amino_acids):
            val = prim_class_data[i]
            df[domain][aa][aa]['sum'] += val
            df[domain][aa][aa]['n'] += 1
            
        for i, lcd_class in enumerate(all_sec_classes):
            aa1, aa2 = lcd_class
            val = sec_class_data[i]
            df[domain][aa1][aa2]['sum'] += val
            df[domain][aa1][aa2]['n'] += 1
            
    h.close()
    
    return df

        
if __name__ == '__main__':
    main()
    