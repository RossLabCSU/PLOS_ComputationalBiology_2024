
import os
amino_acids = 'ACDEFGHIKLMNPQRSTVWY'

def main():

    hq_df = get_HQprots()
    allQ_df = get_all_Qrich_LCDs(hq_df)
    
    contains_nonoverlap_Q_df = crosscheck_bounds(hq_df, allQ_df)
    
    output = open('HQprots_Without_QrichLCDs.tsv', 'w')
    output.write('Proteome\tHQ Proteins Without Q-rich Primary LCD\n')
    for proteome in hq_df:
        hq_prots_only = [prot for prot in hq_df[proteome] if prot not in contains_nonoverlap_Q_df[proteome]]
        output.write('\t'.join([proteome, ','.join(hq_prots_only)]) + '\n')
    output.close()
    
    
def crosscheck_bounds(hq_df, allQ_df):

    contains_nonoverlap_Q_df = {}
    for proteome in hq_df:
        contains_nonoverlap_Q_df[proteome] = contains_nonoverlap_Q_df.get(proteome, [])
        for prot in hq_df[proteome]:
            hq_positions = []
            for bounds in hq_df[proteome][prot]:
                hq_positions += [x for x in range(bounds[0], bounds[1]+1)]
                
            if prot in allQ_df[proteome]:
                q_only_positions = []
                for qbounds in allQ_df[proteome][prot]:
                    all_positions = [x for x in range(qbounds[0], qbounds[1]+1)]
                    non_hq_positions = [x for x in all_positions if x not in hq_positions]
                    q_only_positions += non_hq_positions
                if len(q_only_positions) > 20:
                    contains_nonoverlap_Q_df[proteome].append(prot)
                    
    return contains_nonoverlap_Q_df
      
                    
def get_HQprots():

    h = open(os.path.join('.', 'HQ_LCDs_AllOrganisms.tsv'))
    header = h.readline()
    
    df = {}
    for line in h:
        proteome, domain, prot_id, lcd_seq, lcd_bounds, *remainder = line.rstrip().split('\t')
        
        if domain != 'Eukaryota':
            continue
            
        junk, prot_id, *junk = prot_id.split('|')
        
        lcd_bounds = lcd_bounds[1:-1].split('-')
        start, end = int(lcd_bounds[0]), int(lcd_bounds[1])
        lcd_bounds = (start, end)
        
        df[proteome] = df.get(proteome, {})
        df[proteome][prot_id] = df[proteome].get(prot_id, [])
        df[proteome][prot_id].append(lcd_bounds)
        
    h.close()
    
    return df
    
    
def get_all_Qrich_LCDs(hq_df):
    
    df = {}
    for proteome in hq_df:
        h = open(os.path.join('.', 'Eukaryota', proteome + '_LCDcomposer_RESULTS.tsv'))
        for i in range(11):
            h.readline()
            
        for line in h:
            prot_desc, prot_id, lcd_seq, lcd_bounds, lcd_class, *remainder = line.rstrip().split('\t')
            
            if lcd_class != 'Q':
                continue
                
            lcd_bounds = lcd_bounds[1:-1].split('-')
            start, end = int(lcd_bounds[0]), int(lcd_bounds[1])
            lcd_bounds = (start, end)
            
            df[proteome] = df.get(proteome, {})
            df[proteome][prot_id] = df[proteome].get(prot_id, [])
            df[proteome][prot_id].append(lcd_bounds)
        
        h.close()
        
    return df
    

if __name__ == '__main__':
    main()