
from datetime import datetime
from tqdm import tqdm
domains = ['Archaea', 'Bacteria', 'Eukaryota', 'Viruses']
amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
aa_strings = []
for res1 in amino_acids:
    for res2 in amino_acids:
        if res1 == res2:
            continue
        aa_strings.append(res1+res2)

def main():

    print('Starting:', str(datetime.now()))
    clan_df = get_clan_info()
    pfam_map = get_pfam_map()
    print('Pfam map created:', str(datetime.now()))
    
    output = open('All_LCDclasses_PfamClanMapping.tsv', 'w')
    output.write('\t'.join(['Domain', 'Proteome', 'LCD Class', 'Clan IDs', 'Clan Set', 'Percent of LCD Prots with Clan Association (corresponds to Clan Set column)', 'Maximum Percentage of Proteins Associated with a Single Clan', 'Clans with Max. Percentage', 'Number of Proteins with Highest-Frequency Clan Annotation', 'Number of Proteins in LCD Category']) + '\n')

    for lcd_class in tqdm(aa_strings):
        lcds_df, prot_to_proteome, proteome_to_domain = get_lcd_prots(lcd_class)
        
        df = {}
        for proteome in lcds_df:
            clans = []
            domain = proteome_to_domain[proteome]
            for prot in lcds_df[proteome]:
                clan_set = set()
                if prot in pfam_map:
                    for pfam_acc in pfam_map[prot]:
                        clan_id = clan_df[pfam_acc]
                        clan_set.add(clan_id)
                clans += list(clan_set)
            num_lcd_prots = len(lcds_df[proteome])
            clans = [x for x in clans if x != '']
            clan_setlist = list(set(clans))
            clan_freqs = [clans.count(clan) for clan in clan_setlist]
            percs = [freq / num_lcd_prots * 100 for freq in clan_freqs]
            if len(percs) == 0:
                percs = ['N/A']
                max_perc = 'N/A'
                clans = ['None']
                clans_with_max = ['None']
                num_prots_with_max = 'N/A'
            else:
                max_perc = max(percs)
                clans_with_max = [clan for i, clan in enumerate(clan_setlist) if percs[i] == max_perc]
                num_prots_with_max = max(clan_freqs)

            output.write('\t'.join([domain, proteome, lcd_class, ';'.join(clans), ';'.join(clan_setlist), ';'.join([str(x) for x in percs]), str(max_perc), ';'.join(clans_with_max), str(num_prots_with_max), str(num_lcd_prots)]) + '\n')
            
    output.close()


def get_pfam_map():

    df = {}
    h = open('Pfam-A_regions_LCDprotsOnly.tsv')
    header = h.readline()
    for line in h:
        prot, pfam_acc, start, end = line.rstrip().split('\t')

        df[prot] = df.get(prot, set())
        df[prot].add(pfam_acc)
    h.close()
    
    return df
    
            
def get_clan_info():

    h = open('Pfam-A.clans.tsv')
    header = h.readline()
    df = {}
    
    for line in h:
        pfam_acc, clan_acc, clan_id, pfam_id, pfam_desc = line.rstrip().split('\t')
        df[pfam_acc] = clan_id
        
    h.close()
    
    return df
    
    
def get_lcd_prots(lcd_class):
    
    h = open(lcd_class+'_LCDs_AllOrganisms.tsv')
    header = h.readline()
    df = {}

    prot_to_proteome = {}
    proteome_to_domain = {}
    
    for line in h:
        proteome, domain, prot_desc, *remainder = line.rstrip().split('\t')
        junk, prot_id, *junk = prot_desc.split('|')
        df[proteome] = df.get(proteome, set())
        df[proteome].add(prot_id)
        prot_to_proteome[prot_id] = proteome
        proteome_to_domain[proteome] = domain
        
    h.close()

    return df, prot_to_proteome, proteome_to_domain
    
        
if __name__ == '__main__':
    main()