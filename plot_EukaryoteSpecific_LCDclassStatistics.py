
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

def main():

    minimum_prots_threshold = '10'
    h = open('ProteinLevelStatistics_EukaryoteSpecificLCDclasses.tsv')
    header = h.readline()

    df = {'LCD Class':[],
            'Value':[],
            'Category':[]}
            
    label = 'Number of Organisms with\n' + r'$\geq$' + minimum_prots_threshold + ' LCD-containing Proteins'
    fig_label = 'Fig3A_Number of Organisms with geq' + minimum_prots_threshold + ' LCD-containing Proteins'
    
    for line in h:
        items = line.rstrip().split('\t')
        lcd_class, num_passing_orgs, max_num_prots, org_with_maxprots, max_perc_prots, org_with_maxperc, all_passing_orgs, all_passing_prot_counts, all_prot_counts = items
        num_passing_orgs = int(num_passing_orgs)
        df['LCD Class'].append(lcd_class)
        df['Value'].append(num_passing_orgs)
        df['Category'].append(label)
            
    h.close()

    individual_barplots(df, label, fig_label)
    
    
def individual_barplots(df, label, fig_label):

    df = pd.DataFrame.from_dict(df)

    values = df['Value']
    lcd_classes = df['LCD Class']
    values, lcd_classes = zip(*sorted(zip(values, lcd_classes), reverse=True))
    sns.barplot(x='LCD Class', y='Value', data=df, order=lcd_classes, palette=['#1f77b4'])
    plt.xticks(fontname='Arial', fontsize=10, rotation=90)
    plt.yticks(fontname='Arial', fontsize=10)
    plt.xlabel('LCD Class', fontname='Arial', fontsize=12)
    plt.ylabel(label, fontname='Arial', fontsize=12)
    fig = plt.gcf()

    fig.set_size_inches(4, 6)
    plt.savefig(fig_label + '_EukaryoteSpecific_LCDclassStatistics_Barplot.tif', bbox_inches='tight', dpi=600, pil_kwargs={'compression':'tiff_lzw'})
    plt.close()


if __name__ == '__main__':
    main()