
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
import os
proteomes = ['UP000000589_10090', 'UP000002494_10116', 'UP000005640_9606', 'UP000008227_9823', 'UP000009136_9913', 'UP000002254_9615']
organism_to_filenum = {'Human':'FigS2', 'Bovine':'FigS3B', 'Mouse':'FigS2', 'Rat':'FigS3A', 'Pig':'FigS3C', 'Dog':'FigS3D'}

def main():

    keratin_ids, prot_to_proteome, prot_to_commonname, commonname_to_prot = get_keratin_IDs()
    seqs_df = get_keratin_seqs(prot_to_proteome)
    lcd_df, all_lcd_classes = get_lcd_locations(keratin_ids, prot_to_proteome)

    plotting(lcd_df, seqs_df, all_lcd_classes, prot_to_commonname)
    
    
def get_lcd_locations(keratin_ids, prot_to_proteome):

    all_lcd_classes = set()
    for prot in keratin_ids:
        lcd_classes = keratin_ids[prot]
        for c in lcd_classes:
            all_lcd_classes.add(c)
            
    df = {}
            
    for lcd_class in all_lcd_classes:
        h = open(lcd_class + '_LCDs_AllOrganisms.tsv')
        header = h.readline()
        
        for line in h:
            proteome, domain, prot, lcd_seq, bounds, prim_aa_comp, prim_aa_disp, sec_aa_comp = line.rstrip().split('\t')
            junk, prot, *junk = prot.split('|')
            if prot not in keratin_ids:
                continue
                
            df[prot] = df.get(prot, {})
            df[prot][lcd_class] = df[prot].get(lcd_class, [])
            start, end = bounds[1:-1].split('-')
            start, end = int(start), int(end)
            positions = [x for x in range(start-1, end)]
            df[prot][lcd_class].append(positions)

        h.close()
        
    return df, all_lcd_classes
            
    
def get_keratin_seqs(prot_to_proteome):

    df = {}
    for proteome in proteomes:
        h = open('./Eukaryota/' + proteome + '.fasta')
        for seq_record in SeqIO.parse(h, 'fasta'):
            prot_id = str(seq_record.id)
            junk, prot_id, *junk = prot_id.split('|')
            if prot_id not in prot_to_proteome:
                continue
                
            seq = str(seq_record.seq)
            if seq[-1] == '*':
                seq = seq[:-1]
                
            df[prot_id] = seq
        h.close()
        
    return df
    
    
def get_keratin_IDs():

    h = open('TableS6_ModelOrganisms_KeratinAssociated_CrichSubset.tsv')
    header = h.readline()
    
    df = {}
    prot_to_proteome = {}
    prot_to_commonname = {}
    commonname_to_prot = {}
    
    for line in h:
        proteome, sciname, commonname, lcd_class, proteins = line.rstrip().split('\t')
        proteins = proteins.split(', ')
        for prot in proteins:
            df[prot] = df.get(prot, [])
            df[prot].append(lcd_class)
            prot_to_proteome[prot] = proteome
            prot_to_commonname[prot] = commonname
            commonname_to_prot[commonname] = commonname_to_prot.get(commonname, [])
            commonname_to_prot[commonname].append(prot)
            
    h.close()
    
    return df, prot_to_proteome, prot_to_commonname, commonname_to_prot
    
    
def plotting(lcd_df, seqs_df, all_lcd_classes, prot_to_commonname):

    all_organisms = set([prot_to_commonname[prot] for prot in prot_to_commonname])

    colors = sns.color_palette('Set2')
    colors += sns.color_palette('pastel')
    colors += sns.color_palette('dark')
    colors += sns.color_palette('colorblind')
    
    color_map = {lcd_class:colors[i] for i, lcd_class in enumerate(sorted(list(all_lcd_classes)))}
    
    lineEnd_offset = 6
    shift_magnitude = 0.2
    overlap_shift = shift_magnitude
    
    for org in all_organisms:
        seq_lens = []
        index = 0
        texts = []
        prots = [prot for prot in lcd_df]
        prot_lens = [len(seqs_df[prot]) for prot in prots]
        prot_lens, prots = zip(*sorted(zip(prot_lens, prots)))
        lcd_bounds = []
        labels = []
        xpositions = []
        ypositions = []
        overlap_checks = []
        sequence_lengths = []
        labels_df = {}
        
        for i, prot in enumerate(prots):
            organism = prot_to_commonname[prot]
            if organism != org:
                continue
            master_lcd_positions = []
            
            labels_df[index] = {'xpositions':[],
                            'labels':[],
                            'overlap_checks':[],
                            'sequence_lengths':[]}

            for lcd_class in lcd_df[prot]:
                lcd_positions = lcd_df[prot][lcd_class]
                lcd_bounds.append(lcd_positions)
                labels.append([lcd_class]*len(lcd_positions))
                seq = seqs_df[prot]
                seq_lens.append(len(seq))
                sequence_lengths.append([len(seq)]*len(lcd_positions))

                all_lcd_positions = [pos for pos_list in lcd_positions for pos in pos_list]
                master_lcd_positions += all_lcd_positions

                xpos_row = []
                ypos_row = []
                overlap_row = []
                for j, pos_list in enumerate(lcd_positions):
                    plt.plot(pos_list, [index for x in range(len(pos_list))], color=color_map[lcd_class], linewidth=15, solid_capstyle='butt', alpha=0.5)
                    x_pos = min(pos_list) + ((max(pos_list) - min(pos_list)) / 2)
                    xpos_row.append(x_pos)
                    labels_df[index]['xpositions'].append(x_pos)
                    labels_df[index]['labels'].append(lcd_class)
                    labels_df[index]['sequence_lengths'].append(len(seq))
                    ypos_row.append(index)
                    other_lists = []
                    for c in lcd_df[prot]:
                        if c == lcd_class:
                            continue
                        for l in lcd_df[prot][c]:
                            other_lists.append(l)
                    overlap_check = [sum([1 if pos in other_list else 0 for pos in pos_list]) for other_list in other_lists]

                    if sum(overlap_check) > 0:
                        overlap_row.append(True)
                        labels_df[index]['overlap_checks'].append(True)
                    else:
                        overlap_row.append(False)
                        labels_df[index]['overlap_checks'].append(False)

                plt.text(len(seq)+(len(seq)/100), index, str(len(seq)), ha='left', va='center', fontname='Arial')
                xpositions.append(xpos_row)
                ypositions.append(ypos_row)
                overlap_checks.append(overlap_row)
                
            nonlcd_positions = [x for x in range(len(seq)) if x not in master_lcd_positions]
            nonlcd_positions, new_pos_lists = merge_windows(nonlcd_positions, 1)
            
            for pos_list in new_pos_lists:
                plt.plot([x for x in pos_list], [index for aa in range(len(pos_list))], color='0.5')

            index += 1

        ax = plt.gca()
        
        for i in range(len(labels_df)):
            labels = labels_df[i]['labels']
            xpositions = labels_df[i]['xpositions']
            ypositions = [i]*len(xpositions)
            overlap_checks = labels_df[i]['overlap_checks']
            seqlens = labels_df[i]['sequence_lengths']
            xpositions, ypositions, labels, overlap_checks, seqlens = zip(*sorted(zip(xpositions, ypositions, labels, overlap_checks, seqlens)))
            for j, xpos in enumerate(xpositions):
                label = labels[j]
                is_overlapping = overlap_checks[j]
                ypos = i
                if is_overlapping:
                    ypos += overlap_shift
                    overlap_shift *= -1
                ax.annotate(label, (xpos, ypos), ha='center', va='center')
        
        fig = plt.gcf()
        length = max(seq_lens)
        plt.ylim(-1, index)
        plt.axis('off')
        fig.set_size_inches(7, index/2)
            
        # CAN FAIL FOR SOME PROTEINS IF THEY ARE TOO LONG
        try:
            plt.savefig(organism_to_filenum[org] + '_' + org + '_AllKeratins_IsoformCartoon.tif', bbox_inches='tight', dpi=600, pil_kwargs={'compression':'tiff_lzw'})
        except:
            fig.set_size_inches(2*(length/200), 1)
            plt.savefig(organism_to_filenum[org] + '_' + org + '_AllKeratins_IsoformCartoon.tif', bbox_inches='tight', dpi=600, pil_kwargs={'compression':'tiff_lzw'})
        plt.close()
        
        
def merge_windows(hit_positions, win_size):

    j = 0
    domain_boundaries = []
    new_pos_lists = []
    
    while hit_positions:
        start = hit_positions[0]
        while j < (len(hit_positions) - 1) and ( (hit_positions[j] + win_size) >= ( hit_positions[j+1] ) ):
            j += 1
        
        #GETS THE ENDING POSITION FOR THE CURRENT WINDOW AND STORES THE 2-TUPLE
        end = hit_positions[j]
        domain_boundaries.append( (start , end+win_size) )
        
        if start == 0:
            new_list = [x for x in range(start, end+2)]
        else:
            new_list = [x for x in range(start-1, end+2)]
        new_pos_lists.append(new_list)

        #MODIFIES hit_positions TO DELETE ALL POSITIONS THAT WERE JUST MERGED, THEN RE-SETS j=0 TO START ITERATING AT THE FIRST POSITION IN THE NEW LIST
        hit_positions = hit_positions[j+1:]
        j = 0

    return domain_boundaries, new_pos_lists
        
        
if __name__ == '__main__':
    main()