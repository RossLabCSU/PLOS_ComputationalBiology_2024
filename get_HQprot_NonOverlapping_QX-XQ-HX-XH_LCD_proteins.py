
import os
amino_acids = 'ACDEFGHIKLMNPQRSTVWY'

def main():

    hq_df = get_HQprots()
    allQ_df = get_XQ_QX_HX_XH_LCDs(hq_df)

    contains_nonoverlap_Q_df = crosscheck_bounds(hq_df, allQ_df)

    output = open('HQprots_Without_QX_XQ_HX_XH_LCDs.tsv', 'w')
    output.write('Proteome\tHQ Proteins Without QX, XQ, HX, or XH LCD\n')

    for proteome in hq_df:
        hq_prots_only = [prot for prot in hq_df[proteome] if prot not in contains_nonoverlap_Q_df[proteome]]
        output.write('\t'.join([proteome, ','.join(hq_prots_only)]) + '\n')

    output.close()
    
    
def crosscheck_bounds(hq_df, allQ_df):

    contains_nonoverlap_Q_df = {}
    for proteome in hq_df:
        contains_nonoverlap_Q_df[proteome] = contains_nonoverlap_Q_df.get(proteome, set())
        for prot in hq_df[proteome]:
            hq_positions = []
            for bounds in hq_df[proteome][prot]:
                hq_positions += [x for x in range(bounds[0], bounds[1]+1)]
                
            for lcd_class in allQ_df:
                if proteome not in allQ_df[lcd_class]:
                    continue

                if prot in allQ_df[lcd_class][proteome]:
                    q_only_positions = []
                    for qbounds in allQ_df[lcd_class][proteome][prot]:
                        all_positions = [x for x in range(qbounds[0], qbounds[1]+1)]
                        non_hq_positions = [x for x in all_positions if x not in hq_positions]
                        q_only_positions += non_hq_positions
                    if len(q_only_positions) > 20:
                        contains_nonoverlap_Q_df[proteome].add(prot)

    return contains_nonoverlap_Q_df
    
    
def get_XQ_QX_HX_XH_LCDs(hq_df):
    
    df = {}
    
    for prim_aa in ['Q', 'H']:
        for sec_aa in amino_acids:
            lcd_class = prim_aa + sec_aa
            if prim_aa == sec_aa or lcd_class in ['QH', 'HQ']:
                continue

            h = open(os.path.join('SecondaryLCDs_by_LCDcategory', prim_aa+sec_aa + '_LCDs_AllOrganisms.tsv'))
            header = h.readline()
            
            for line in h:
                proteome, domain, prot_id, lcd_seq, lcd_bounds, *remainder = line.rstrip().split('\t')
                
                if domain != 'Eukaryota':
                    continue
                    
                junk, prot_id, *junk = prot_id.split('|')
                
                lcd_bounds = lcd_bounds[1:-1].split('-')
                start, end = int(lcd_bounds[0]), int(lcd_bounds[1])
                lcd_bounds = (start, end)
                
                df[lcd_class] = df.get(lcd_class, {})
                df[lcd_class][proteome] = df[lcd_class].get(proteome, {})
                df[lcd_class][proteome][prot_id] = df[lcd_class][proteome].get(prot_id, [])
                df[lcd_class][proteome][prot_id].append(lcd_bounds)
                
            h.close()
            
    for prim_aa in amino_acids:
        for sec_aa in ['Q', 'H']:
            lcd_class = prim_aa + sec_aa
            if prim_aa == sec_aa or lcd_class in ['QH', 'HQ']:
                continue

            h = open(os.path.join('SecondaryLCDs_by_LCDcategory', prim_aa+sec_aa + '_LCDs_AllOrganisms.tsv'))
            header = h.readline()
            
            for line in h:
                proteome, domain, prot_id, lcd_seq, lcd_bounds, *remainder = line.rstrip().split('\t')
                
                if domain != 'Eukaryota':
                    continue
                    
                junk, prot_id, *junk = prot_id.split('|')
                
                lcd_bounds = lcd_bounds[1:-1].split('-')
                start, end = int(lcd_bounds[0]), int(lcd_bounds[1])
                lcd_bounds = (start, end)
                
                df[lcd_class] = df.get(lcd_class, {})
                df[lcd_class][proteome] = df[lcd_class].get(proteome, {})
                df[lcd_class][proteome][prot_id] = df[lcd_class][proteome].get(prot_id, [])
                df[lcd_class][proteome][prot_id].append(lcd_bounds)
                
            h.close()
            
    return df
                    
                    
def get_HQprots():

    h = open(os.path.join('SecondaryLCDs_by_LCDcategory', 'HQ_LCDs_AllOrganisms.tsv'))
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


if __name__ == '__main__':
    main()