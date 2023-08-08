
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import statistics
from scipy import stats
amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
aa_list = list(amino_acids)
ordered_aas_df = {'Archaea':'LVEAGISDRTKPFNYQMHWC',
                'Bacteria':'LAGVEISDRTKPFQNYMHWC',
                'Eukaryota':'LSAEVGKRTPDINQFYHMCW',
                'Viruses':'LASVGETKDIRNPFQYMHWC'}
ordered_aas = 'ILVFWYMCQNPGASTDEKRH'

domains = ['Archaea','Bacteria','Eukaryota','Viruses']
aa_strings = []
for res1 in amino_acids:
    for res2 in amino_acids:
        aa_strings.append(res1+res2)

def main():

    proteome1 = 'UP000005640_9606'

    for proteome2 in ['UP000002311_559292', 'AVERAGE']:
        domain_of_interest = get_domain_of_interest(proteome1)
        ordered_aas = ordered_aas_df[domain_of_interest]

        if proteome2.upper() == 'AVERAGE':
            df, val_df = get_data([proteome1])
            mean_matrix = calc_mean_matrix(proteome1)
            df['AVERAGE'] = mean_matrix
        else:
            df, val_df = get_data([proteome1, proteome2])

        diff_matrix = calc_matrix_diffs(df[proteome1], df[proteome2])

        diag_vals = get_diagonal_vals(diff_matrix)
        diag_vals, sorted_prim_aas = zip(*sorted(zip(diag_vals, aa_list), reverse=True))

        diff_matrix = list(np.transpose(diff_matrix))
        diff_matrix = [list(row) for row in diff_matrix]
        
        max_offdiag, min_offdiag, ondiag_vals, offdiag_vals = get_nondiag_max(diff_matrix)
        center = 0
        
        plotting_df = pd.DataFrame(data=diff_matrix, index=aa_list, columns=aa_list)
        plotting_df, sorted_aas = sort_columns_ave_value(plotting_df, aa_list)

        plotting(plotting_df, proteome1, proteome2, max_offdiag, min_offdiag, center, mask_diagonal=True)
        
        proteome1_matrix = np.transpose(df[proteome1])
        proteome2_matrix = np.transpose(df[proteome2])
        
        matrices = [proteome1_matrix, proteome2_matrix, diff_matrix]
        labels = [proteome1, proteome2, 'Difference']
        output_matrices(matrices, labels, sorted_aas)
        
        data_labels = [proteome1 + ' vs ' + proteome2]*20
        
        if proteome2 != 'AVERAGE':
            plot_primAA_barplots(diag_vals, sorted_prim_aas, data_labels, 'Difference', proteome1, proteome2)


def plot_primAA_barplots(vals, sorted_prim_aas, data_labels, yaxis_label, proteome1, proteome2):

    colors = ['#05696b', '#fd411e']
    
    aa_encoding = {aa:i for i, aa in enumerate(sorted_prim_aas)}
    df = {'Value':vals,
        'LCD Class':[aa_encoding[aa] for aa in sorted_prim_aas],
        'Labels':data_labels}
    colors = sns.color_palette()
    df = pd.DataFrame.from_dict(df)
    sns.barplot(x='LCD Class', y='Value', data=df, hue='Labels')
    plt.xticks([x for x in range(20)], labels=sorted_prim_aas[:20], rotation=0, fontname='Arial', fontsize=10)
    plt.yticks(rotation=0, fontname='Arial', fontsize=10)
    plt.xlabel('Primary LCD Class', fontname='Arial', fontsize=12)
    plt.ylabel(yaxis_label, fontname='Arial', fontsize=12)
    if yaxis_label == 'Percentile':
        plt.ylim(0, 105)
        
    ax = plt.gca()
    ax.get_legend().remove()

    fig = plt.gcf()
    fig.set_size_inches(5, 3.25)
    if proteome2 == 'Percentile':
        plt.savefig('Fig10A_' + proteome1 + '_' + yaxis_label + '_PrimaryLCD_PerResidueOccupancy_Barplot.tif', bbox_inches='tight', dpi=600)
    else:
        plt.savefig('Fig10B_' + proteome1 + '_vs_' + proteome2 + '_' + yaxis_label + '_PrimaryLCD_PerResidueOccupancy_Barplot.tif', bbox_inches='tight', dpi=600)
    plt.close()
    
    
def get_diagonal_vals(matrix):

    diag = []
    for i, row in enumerate(matrix):
        for j, col in enumerate(matrix[i]):
            if i == j:
                diag.append(matrix[i][j])
                
    return diag
    
    
def output_matrices(matrices, labels, sorted_aas):
    
    for i, matrix in enumerate(matrices):
        label = labels[i]
        
        df = pd.DataFrame(data=matrix, index=aa_list, columns=aa_list)
        df = df.reindex(index=sorted_aas, columns=sorted_aas)
        
        df.to_csv('Fig10_' + label + '_HeatmapMatrix.csv')
        
    
def sort_columns_ave_value(df, aa_list):

    mean_values = [statistics.mean(df[aa]) for aa in aa_list]
    mean_values, sorted_aas = zip(*sorted(zip(mean_values, aa_list), reverse=True))
    
    df = df.reindex(index=sorted_aas, columns=sorted_aas)
    
    return df, sorted_aas
    
    
def calc_matrix_diffs(mat1, mat2):
    
    diff_mat = []
    for i, row in enumerate(mat1):
        diff_row = []
        for j, col in enumerate(mat1[i]):
            diff = mat1[i][j] - mat2[i][j]
            diff_row.append(diff)
            
        diff_mat.append(diff_row)
        
    return diff_mat
    

def main_percentile():

    proteome1 = 'UP000005640_9606'
    proteome2 = 'UP000002311_559292'

    df, val_df = get_data([proteome1])
    percentile_matrix = calc_percentiles(proteome1, val_df[proteome1])

    percentile_matrix = list(np.transpose(percentile_matrix))
    percentile_matrix = [list(row) for row in percentile_matrix]
    
    diag_vals = get_diagonal_vals(percentile_matrix)

    if proteome2 != 'AVERAGE':
        df2, val_df2 = get_data([proteome2])
        percentile_matrix2 = calc_percentiles(proteome2, val_df2[proteome2])
        diag_vals2 = get_diagonal_vals(percentile_matrix2)
        sorted_prim_aas = aa_list + aa_list
        data_labels = [proteome1]*20 + [proteome2]*20
        diffs = [diag_vals[i] - diag_vals2[i] for i in range(len(diag_vals))]
        sorted_diffs, sorted_prim_aas = zip(*sorted(zip(diffs, aa_list), reverse=True))
        sorted_indices = [amino_acids.index(aa) for aa in sorted_prim_aas]
        sorted_diag_vals = [diag_vals[index] for index in sorted_indices]
        sorted_diag_vals2 = [diag_vals2[index] for index in sorted_indices]
        diag_vals = sorted_diag_vals + sorted_diag_vals2
        sorted_prim_aas += sorted_prim_aas

    else:
        diag_vals, sorted_prim_aas = zip(*sorted(zip(diag_vals, aa_list), reverse=True))
        data_labels = [proteome1]*20
    
    max_offdiag, min_offdiag, ondiag_vals, offdiag_vals = get_nondiag_max(percentile_matrix)
    center = 50
    
    plotting_df = pd.DataFrame(data=percentile_matrix, index=aa_list, columns=aa_list)
    plotting_df, sorted_aas = sort_columns_ave_value(plotting_df, aa_list)

    plotting(plotting_df, proteome1, 'PERCENTILE', 100, 0, center, mask_diagonal=False)
    
    matrices = [percentile_matrix]
    labels = ['Percentile']
    output_matrices(matrices, labels, sorted_aas)

    plot_primAA_barplots(diag_vals, sorted_prim_aas, data_labels, 'Percentile', proteome1, 'Percentile')
    
    
def get_nondiag_max(matrix):
    
    offdiag_vals = []
    ondiag_vals = []
    for i, row in enumerate(matrix):
        for j, col in enumerate(matrix[i]):
            if i == j:
                ondiag_vals.append(matrix[i][j])
            else:
                offdiag_vals.append(matrix[i][j])

    return max(offdiag_vals), min(offdiag_vals), ondiag_vals, offdiag_vals
 
 
def plotting(df, proteome1, proteome2, max_offdiag, min_offdiag, center, mask_diagonal=False):

    mask = []
    offdiag_vals = []
    for i in range(20):
        row = []
        for j in range(20):
            if i == j:
                row.append(1)
            else:
                row.append(0)
                
        mask.append(row)
     
    if mask_diagonal:
        cg = sns.heatmap(df, mask=np.array(mask), vmin=min_offdiag, vmax=max_offdiag, cmap='coolwarm', center=center)
    else:
        max_val = df.max().max()
        cg = sns.heatmap(df, vmin=min_offdiag, vmax=max_val, cmap='coolwarm', center=center)
    ax = plt.gca()
    ax.set_facecolor("grey")
    
    plt.xticks(rotation=0, fontname='Arial', fontsize=10)
    plt.yticks(rotation=0, fontname='Arial', fontsize=10)
    plt.xlabel('Primary Amino Acid', fontname='Arial', fontsize=12)
    plt.ylabel('Secondary Amino Acid', fontname='Arial', fontsize=12)
    plt.ylim(20, 0)
    fig = plt.gcf()
    if proteome2 == 'AVERAGE':
        plt.savefig('Fig10E_' + proteome1 + '_vs_' + proteome2 + '_PerResidueLCDoccupancy_DifferenceMatrix.tif', bbox_inches='tight', dpi=600)
    else:
        if proteome2 == 'PERCENTILE':
            plt.savefig('Fig10D_' + proteome1 + '_vs_' + proteome2 + '_PerResidueLCDoccupancy_DifferenceMatrix.tif', bbox_inches='tight', dpi=600)
        else:
            plt.savefig('Fig10C_' + proteome1 + '_vs_' + proteome2 + '_PerResidueLCDoccupancy_DifferenceMatrix.tif', bbox_inches='tight', dpi=600)
    plt.close()

    
def get_data(proteomes):

    h = open('LCDfrequency_data_PerResidue_LCDoccupancy.tsv')
    header = h.readline()

    df = {}
    val_df = {}
        
    for line in h:
        proteome, taxid, domain, sci_name, common_name, total_prots, *data = line.rstrip().split('\t')

        if proteome not in proteomes:
            continue
            
        primary_lcd_vals = data[:20]
        primary_lcd_vals = [float(x) for x in primary_lcd_vals]
        secondary_lcd_vals = data[20:]
        secondary_lcd_vals = [float(x) for x in secondary_lcd_vals]
        
        matrix = np.reshape(secondary_lcd_vals, (20, 19))
        
        matrix = [list(row) for row in matrix]
        
        for i, aa in enumerate(amino_acids):    
            matrix[i].insert(i, primary_lcd_vals[i])

        df[proteome] = matrix
        
        val_df[proteome] = {}
        for i, res1 in enumerate(amino_acids):
            for j, res2 in enumerate(amino_acids):
                val = matrix[i][j]
                val_df[proteome][res1+res2] = val

    h.close()

    return df, val_df
    
    
def get_domain_of_interest(proteome1):
    
    h = open('LCDfrequency_data_PerResidue_LCDoccupancy.tsv')
    header = h.readline()

    df = {}

    for line in h:
        proteome, taxid, domain, sci_name, common_name, total_prots, *data = line.rstrip().split('\t')
        
        if proteome != proteome1:
            continue
            
        domain_of_interest = domain
    h.close()
    
    return domain_of_interest
    
    
def calc_mean_matrix(proteome1):
    
    domain_of_interest = get_domain_of_interest(proteome1)
    
    h = open('LCDfrequency_data_PerResidue_LCDoccupancy.tsv')
    header = h.readline()

    df = {aa_string:[] for aa_string in aa_strings}

    for line in h:
        proteome, taxid, domain, sci_name, common_name, total_prots, *data = line.rstrip().split('\t')

        if domain != domain_of_interest:
            continue
            
        primary_lcd_vals = data[:20]
        primary_lcd_vals = [float(x) for x in primary_lcd_vals]
        secondary_lcd_vals = data[20:]
        secondary_lcd_vals = [float(x) for x in secondary_lcd_vals]
        
        matrix = np.reshape(secondary_lcd_vals, (20, 19))
        matrix = [list(row) for row in matrix]
        
        for i, aa in enumerate(amino_acids):    
            matrix[i].insert(i, primary_lcd_vals[i])
        
        for i, res1 in enumerate(amino_acids):
            row = matrix[i]
            for j, res2 in enumerate(amino_acids):
                val = matrix[i][j]
                df[res1+res2].append(val)

    h.close()
    
    matrix = []
    for i, res1 in enumerate(amino_acids):
        row = []
        for j, res2 in enumerate(amino_acids):
            vals = df[res1+res2]
            ave_val = statistics.mean(vals)
            row.append(ave_val)
        matrix.append(row)

    return matrix
    
    
def calc_percentiles(proteome1, proteome1_df):
    
    domain_of_interest = get_domain_of_interest(proteome1)
    
    h = open('LCDfrequency_data_PerResidue_LCDoccupancy.tsv')
    header = h.readline()

    df = {aa_string:[] for aa_string in aa_strings}

    for line in h:
        proteome, taxid, domain, sci_name, common_name, total_prots, *data = line.rstrip().split('\t')

        if domain != domain_of_interest:
            continue
            
        primary_lcd_vals = data[:20]
        primary_lcd_vals = [float(x) for x in primary_lcd_vals]
        secondary_lcd_vals = data[20:]
        secondary_lcd_vals = [float(x) for x in secondary_lcd_vals]
        
        matrix = np.reshape(secondary_lcd_vals, (20, 19))
        matrix = [list(row) for row in matrix]

        for i, aa in enumerate(amino_acids):    
            matrix[i].insert(i, primary_lcd_vals[i])
        
        for i, res1 in enumerate(amino_acids):
            row = matrix[i]
            for j, res2 in enumerate(amino_acids):
                val = matrix[i][j]
                df[res1+res2].append(val)

    h.close()
    
    matrix = []
    for i, res1 in enumerate(amino_acids):
        row = []
        for j, res2 in enumerate(amino_acids):
            vals = df[res1+res2]
            proteome_val = proteome1_df[res1+res2]
            percentile_score = stats.percentileofscore(vals, proteome_val)
            row.append(percentile_score)
        matrix.append(row)

    return matrix
        
if __name__ == '__main__':
    main()
    main_percentile()
