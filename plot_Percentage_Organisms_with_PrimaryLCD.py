
import matplotlib.pyplot as plt
import seaborn as sns

def main():

    h = open('PrimaryLCDs_NumberOfOrganisms.tsv')
    header = h.readline()
    
    df = {'Domain':[],
        'LCD Class':[],
        'Percentage':[]}
    for line in h:
        domain, lcd_class, num_with_lcd, total, perc = line.rstrip().split('\t')
        df['Domain'].append(domain)
        df['LCD Class'].append(lcd_class)
        df['Percentage'].append(float(perc))
        
    h.close()
    
    plotting(df)
    
        
def plotting(df):
    
    colors = sns.color_palette()
    sns.barplot(x='LCD Class', y='Percentage', data=df, hue='Domain', palette=colors)

    plt.xticks(fontname='Arial', fontsize=12)
    plt.yticks(fontname='Arial', fontsize=12)
    plt.xlabel('LCD Class', fontname='Arial', fontsize=14)
    plt.ylabel('Percentage of Organisms\nwith' + r'$\geq$' + '1 LCD', fontname='Arial', fontsize=14)
    ax = plt.gca()
    ax.legend(bbox_to_anchor=(1,1), loc=2, handletextpad=0.3)
    fig = plt.gcf()
    fig.set_size_inches(12, 4)
    plt.savefig('Fig1A_PrimaryLCDs_PercOrganismsWithLCD.tif', bbox_inches='tight', dpi=600, pil_kwargs={'compression':'tiff_lzw'})
    plt.close()
        

if __name__ == '__main__':
    main()