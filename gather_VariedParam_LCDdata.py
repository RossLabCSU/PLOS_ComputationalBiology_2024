

import os
from tqdm import tqdm
cwd = os.getcwd()
domains = ['Archaea', 'Bacteria', 'Eukaryota', 'Viruses']

def main():

    cwd = os.getcwd()
    
    output = open('LCDcomposer_VariedParameters_Results.tsv', 'w')
    output.write('\t'.join(['Domain', 'Organism ID', 'Window Size', 'Primary Amino Acid', 'Secondary Amino Acid', 'Primary Composition Threshold', 'Secondary Composition Threshold', '# of Proteins with LCD']) + '\n')
    
    data_types = ['VariedCompThresholds', 'VariedWindowSize']
    for i, data_type in enumerate(data_types):
        for domain in domains:
            dirname = domain + '_' + data_type
            base_path = os.path.join(cwd, 'RandomlySelectedOrganisms', dirname)
        
            files = os.listdir(os.path.join(base_path))
            files = [x for x in files if x.endswith('RESULTS.tsv')]
            for file in tqdm(files):
                fileinfo = file.split('_')
                proteome = fileinfo[0] + '_' + fileinfo[1]
                lcd_class = fileinfo[3]
                if data_type == 'VariedCompThresholds':
                    window_size = '20'
                    if len(lcd_class) == 1:
                        prim_aa = lcd_class
                        sec_aa = 'none'
                        prim_threshold = fileinfo[4]
                        sec_threshold = 'none'
                    else:
                        prim_aa, sec_aa = lcd_class.split('-')
                        prim_threshold, sec_threshold = fileinfo[4].split('-')
                        if int(sec_threshold) > int(prim_threshold):
                            prim_aa, sec_aa = sec_aa, prim_aa
                            prim_threshold, sec_threshold = sec_threshold, prim_threshold
                else:
                    window_size = fileinfo[4].replace('windowSize', '')
                    if window_size == '20':
                        continue
                    if len(lcd_class) == 1:
                        prim_aa = lcd_class
                        sec_aa = 'none'
                        prim_threshold = '40'
                        sec_threshold = 'none'
                    else:
                        prim_aa, sec_aa = lcd_class.split('-')
                        prim_threshold = '40'
                        sec_threshold = '20'
                num_prots_with_lcd = calc_num_prots_with_lcd(file, base_path)
                output.write('\t'.join([domain, proteome, window_size, prim_aa, sec_aa, prim_threshold, sec_threshold, str(num_prots_with_lcd)]) + '\n')
                if prim_threshold == sec_threshold and prim_aa != sec_aa:
                    output.write('\t'.join([domain, proteome, window_size, sec_aa, prim_aa, sec_threshold, prim_threshold, str(num_prots_with_lcd)]) + '\n')
    output.close()
    
                
def calc_num_prots_with_lcd(file, base_path):
    h = open(os.path.join(base_path, file))
    for i in range(11):
        h.readline()
        prots = set()
        
    for line in h:
        header, prot, *remainder = line.rstrip().split('\t')
        prots.add(prot)
        
    h.close()
    
    return len(prots)


if __name__ == '__main__':
    main()