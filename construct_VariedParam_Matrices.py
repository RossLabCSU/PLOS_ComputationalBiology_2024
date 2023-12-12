
domains = ['Archaea', 'Bacteria', 'Eukaryota', 'Viruses']
window_sizes = [x for x in range(20, 70, 10)]
comp_thresholds = [x for x in range(30, 70, 10)]
amino_acids = 'ACDEFGHIKLMNPQRSTVWY'

def main():

    data_types = ['VariedCompThreshold', 'VariedWindowSize']
    df = get_data()
    
    construct_varied_compthreshold_matrix(df)
    construct_varied_windowsize_matrix(df)
    
    
def construct_varied_windowsize_matrix(df):

    output = open('AllDomains_VariedWindowSize_OrganismLevel_LCDfrequencies.tsv', 'w')
    labels = []
    for prim_aa in amino_acids:
        for window_size in window_sizes:
            labels.append(prim_aa+'_'+str(window_size)+'aaWindowSize')
    output.write('\t'.join(['Secondary LCD Classes'] + labels) + '\n')
    combined_matrix = []
    
    prim_threshold = 40
    sec_threshold = 20
    
    for sec_aa in amino_acids:  # SPECIFIC ORDERING TO MAKE THE HEATMAP EASIER TO INTERPRET
        for domain in ['Viruses', 'Archaea', 'Bacteria', 'Eukaryota']:  # SPECIFIC ORDERING TO MAKE THE HEATMAP EASIER TO INTERPRET
            row = []
            for prim_aa in amino_acids:
                for window_size in window_sizes:
                    if prim_aa == sec_aa:
                        vals = df[domain][window_size][prim_aa][prim_threshold][prim_aa][prim_threshold]
                    else:
                        if sec_threshold > prim_threshold:
                            vals = df[domain][window_size][sec_aa][sec_threshold][prim_aa][prim_threshold]
                        else:
                            vals = df[domain][window_size][prim_aa][prim_threshold][sec_aa][sec_threshold]

                    perc = sum(vals) / len(vals) * 100
                    row.append(str(perc))
            combined_matrix.append(row)
            output.write('\t'.join([domain+'_'+sec_aa] + row) + '\n')
    output.close()

    
def construct_varied_compthreshold_matrix(df):

    window_size = 20
    for domain in domains:
        output = open(domain + '_VariedCompThresholds_OrganismLevel_LCDfrequencies.tsv', 'w')
        labels = []
        for aa in amino_acids:
            for prim_threshold in comp_thresholds:
                labels.append(aa+'_'+str(prim_threshold))
        output.write('\t'.join(['Secondary LCD Classes'] + labels) + '\n')
        combined_matrix = []
        for sec_aa in amino_acids:
            for sec_threshold in [20] + comp_thresholds:
                row = []
                for prim_aa in amino_acids:
                    for prim_threshold in comp_thresholds:
                        if prim_threshold + sec_threshold > 100:
                            row.append('-1000')
                            continue
                        
                        if prim_aa == sec_aa:
                            vals = df[domain][window_size][prim_aa][prim_threshold][prim_aa][prim_threshold]
                        else:
                            if sec_threshold > prim_threshold:
                                vals = df[domain][window_size][sec_aa][sec_threshold][prim_aa][prim_threshold]
                            else:
                                vals = df[domain][window_size][prim_aa][prim_threshold][sec_aa][sec_threshold]

                        perc = sum(vals) / len(vals) * 100
                        row.append(str(perc))
                combined_matrix.append(row)
                output.write('\t'.join([str(sec_aa) + '_' + str(sec_threshold)] + row) + '\n')
        output.close()


def get_data():

    h = open('LCDcomposer_VariedParameters_Results.tsv')

    header = h.readline()
    
    # DEEP DICTIONARY WITH THE FOLLOWING LAYERS:
    # 1. DOMAIN OF LIFE
    #   2. WINDOW SIZE
    #     3. PRIMARY AA
    #       4. PRIMARY AA THRESHOLD
    #         5. SECONDARY AA
    #           6. SECONDARY AA THRESHOLD
    #             7. VALUES = LIST WITH 0's AND 1's INDICATING WHETHER THE ORGANISM CONTAINED AT LEAST ONE LCD FOR THE GIVEN PARAMETER COMBINATION
    df = {domain:{window_size:{res1:{threshold1:{res2:{threshold2:[] for threshold2 in [20] + comp_thresholds} for res2 in amino_acids} for threshold1 in comp_thresholds} for res1 in amino_acids} for window_size in window_sizes} for domain in ['Archaea', 'Bacteria', 'Eukaryota', 'Viruses']}

    for line in h:
        domain, proteome, window_size, prim_aa, sec_aa, prim_threshold, sec_threshold, num_prots = line.rstrip().split('\t')

        if sec_aa == 'none':
            sec_aa = prim_aa[:]
            sec_threshold = prim_threshold[:]
        if int(num_prots) == 0:
            contains_lcd = 0
        else:
            contains_lcd = 1

        df[domain][int(window_size)][prim_aa][int(prim_threshold)][sec_aa][int(sec_threshold)].append(contains_lcd)

    h.close()

    return df

if __name__ == '__main__':
    main()