
import matplotlib.pyplot as plt
import seaborn as sns
import statistics
import pandas as pd

domains = ['Archaea', 'Bacteria', 'Eukaryota', 'Viruses']
amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
aa_strings = []
for res1 in amino_acids:
    for res2 in amino_acids:
        if res1 == res2:
            continue
        aa_strings.append(res1+res2)
        
def main():
    
    num_lcdprots_threshold = 5
    for domain_of_interest in domains:
        h = open('All_LCDclasses_PfamClanMapping.tsv')
        header = h.readline()
        df = {}
        for line in h:
            domain, proteome, lcd_class, clan_ids, clan_set, perc_lcds_with_clan, max_perc, clans_with_max_perc, num_prots_with_max, num_lcdprots = line.rstrip().split('\t')
            if domain != domain_of_interest:
                continue
                
            if int(num_lcdprots) < num_lcdprots_threshold or max_perc == 'N/A':
                continue
                
            df[proteome] = df.get(proteome, {})
            df[proteome][lcd_class] = float(max_perc)
            
        h.close()

        matrix = make_matrix(df, domain_of_interest)
        matrix_df = pd.DataFrame(matrix)
        plotting(matrix_df, matrix, domain_of_interest, mask_diagonal=True)
        
        
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
        
        
def make_matrix(df, domain):
    
    matrix = []
    for res1 in amino_acids:
        row = []
        for res2 in amino_acids:
            lcd_class = res1 + res2
            if res1 == res2:
                row.append(-1000)
                continue
            vals = [df[proteome][lcd_class] for proteome in df if lcd_class in df[proteome]]
            if len(vals) == 0:
                ave = -1000
            else:
                ave = statistics.mean(vals)
            row.append(ave)
        matrix.append(row)

    return matrix
    
    
def plotting(df, matrix, domain, mask_diagonal=False):

    mask = []
    offdiag_vals = []
    for i in range(20):
        row = []
        for j in range(20):
            if i == j or matrix[i][j] == -1000:
                row.append(1)
            else:
                row.append(0)
                
        mask.append(row)
        
    mask = pd.DataFrame(mask)
    df = df.transpose()
    
    if mask_diagonal:
        cg = sns.heatmap(df, mask=mask, vmin=0, vmax=100, cmap='coolwarm', center=50)
    else:
        cg = sns.heatmap(df, vmin=0, vmax=100, cmap='coolwarm', center=50)
    ax = plt.gca()
    ax.set_facecolor("grey")
    plt.xticks([x+0.5 for x in range(20)], labels=list(amino_acids), rotation=0, fontname='Arial', fontsize=10)
    plt.yticks([x+0.5 for x in range(20)], labels=list(amino_acids), rotation=0, fontname='Arial', fontsize=10)
    plt.xlabel('Primary Amino Acid', fontname='Arial', fontsize=12)
    plt.ylabel('Secondary Amino Acid', fontname='Arial', fontsize=12)

    plt.ylim(20, 0)
    fig = plt.gcf()
    plt.savefig('FigS14B_' + domain + '_MaxPercentLCDproteins_with_SinglePfamClan.tif', bbox_inches='tight', dpi=600)
    plt.close()
    
        
if __name__ == '__main__':
    main()
