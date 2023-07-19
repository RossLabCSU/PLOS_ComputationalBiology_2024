
from goatools import obo_parser
from goatools.go_enrichment import GOEnrichmentStudy
from Bio import SeqIO
from scipy import stats
import os
import sys

def main():

    gaf_df, proteome_to_organism = get_gaf_df()
    hq_prots_df = get_HQprots(gaf_df)

    go_obo = 'go-basic.obo'
    go = obo_parser.GODag(os.path.join('GOAfiles', go_obo))

    win_size = 20
    
    index = 0
    num_go_completed = 0
    for proteome in gaf_df:
        index += 1
        organism = proteome_to_organism[proteome].replace(' ', '-')
        print('Running: ', organism)
        gaf_file = gaf_df[proteome]
        assoc = get_assoc(gaf_file)
        background_prots = get_background_prots(proteome, win_size)

        lcd_prots = hq_prots_df[proteome]
        if len(lcd_prots) == 0:
            continue
        methods = ["bonferroni", "sidak", "holm", "fdr"]
        organism = organism.replace('(', '').replace(')', '').replace('/', '')
        logfile = open('HQprots-WITHOUT-QX-XQ-HX-XH-LCDs_GOruns_'+organism+'_LogFile.txt', 'w')
        g = GOEnrichmentStudy(background_prots, assoc, go, propagate_counts=True, alpha=0.05, log=logfile)
        g_res = g.run_study(lcd_prots)
        filename = proteome + '_' + organism + '_HQprots-WITHOUT-QX-XQ-HX-XH-LCDs_GOterm_RESULTS.tsv' #change file type here to txt or xlsx if desired
        output_GOresults(g_res, filename)
        print('End: ', organism)
        num_go_completed += 1
        
    print('Total number of organisms analyzed by GO analysis:', num_go_completed)
    
    
def get_gaf_df():

    h = open('GOA_FileMap_HQ_LCDs.tsv')
    df = {}
    df2 = {}
    h.readline()
    
    for line in h:
        proteome, sciname, goa_file = line.rstrip().split('\t')
        df[proteome] = goa_file
        df2[proteome] = sciname
        
    h.close()
    
    return df, df2
    
    
def get_HQprots(gaf_df):

    h = open('HQprots_Without_QX_XQ_HX_XH_LCDs.tsv')
    header = h.readline()
    
    df = {}
    
    for line in h:
        items = line.rstrip().split('\t')
        proteome = items[0]
        if proteome not in gaf_df or len(items) == 1:   # WHEN len(items) == 1, THERE IS AN EMPTY PROTEIN LIST
            continue

        prots = items[1]

        prots = prots.split(',')
        df[proteome] = df.get(proteome, set(prots))
        
    h.close()
    
    return df


def output_GOresults(g_res, filename):

    output = open(filename, 'w')
    new_header = ['# GO', 'NS', 'enrichment', 'name', 'hits_in_study', 'totalProteins_in_study', 'ratio_in_study', 'hits_in_population', 'totalProteins_in_population', 'ratio_in_pop', 'odds_ratio', 'p_uncorrected', 'depth', 'study_count', 'p_bonferroni', 'p_sidak', 'p_holm', 'study_items']

    output.write('\t'.join(new_header) + '\n')

    lines = []
    sidak_pvals = []
    oddsratios = []
    for enr_record in g_res:
        items = str(enr_record).split('\t')
        go_id = items[0]

        ratio_in_study, ratio_in_pop = items[4:6]
        hits_study, total_study = ratio_in_study.split('/')
        hits_pop, total_pop = ratio_in_pop.split('/')
        sidak_pval = float(items[10])
        
        contingency = [[int(hits_study), (int(total_study)-int(hits_study))], [int(hits_pop), (int(total_pop)-int(hits_pop))]]
        
        oddsratio, pval = stats.fisher_exact(contingency, alternative='two-sided')
        
        new_line = [go_id] + items[1:4] + [hits_study, total_study] + [items[4]] + [hits_pop, total_pop] + [items[5]] + [str(oddsratio)] + items[6:]

        lines.append(new_line)
        sidak_pvals.append(sidak_pval)
        oddsratios.append(oddsratio)
        
    # PERFORM SEQUENTIAL SORTING. SORTS BY oddsratio (HIGH-TO-LOW), THEN SORTS BY sidak_pval (LOW-TO-HIGH). END RESULT IS THAT LOWEST P-VALUES SHOW UP FIRST, THEN FOR ALL LINES WITH THE SAME P-VALUE, THOSE WITH THE HIGHEST oddsratio SHOW UP FIRST.
    if len(oddsratios) == 0:
        lines = []
    else:
        oddsratios, lines, sidak_pvals = zip(*sorted(zip(oddsratios, lines, sidak_pvals), reverse=True))
        heirarchy = [x for x in range(len(sidak_pvals))]
        sidak_pvals, heirarchy, oddsratios, lines = zip(*sorted(zip(sidak_pvals, heirarchy, oddsratios, lines)))
    for line in lines:
        output.write('\t'.join(line) + '\n')

    output.close()
    
        
def get_background_prots(proteome, win_size):
    
    prots = []
    prot_lengths = {}
    h = open(os.path.join('.', 'Eukaryota', proteome + '.fasta'))
    for seq_record in SeqIO.parse(h, 'fasta'):
        id = str(seq_record.id)
        junk, uniprot_id, junk = id.split('|')
        
        # ACCOUNTS FOR WINDOW SIZE SINCE ONLY PROTEINS EXCEEDING THE WINDOW SIZE COULD POSSIBLY CONTAIN AN LCD.
        seq = str(seq_record.seq)
        if seq[-1] == '*':
            seq = seq[:-1]
        prot_length = len(seq)
        if prot_length < win_size:
            continue
            
        prots.append( uniprot_id )

    return prots
    
    
def get_assoc(gaf_file):
    """
    Reads the proper gene association file (GAF).
    Return is a dictionary with Protein IDs as keys, and a set of associated GO terms as values (this is a formal Python set, not a list)
    
    Each list represents a separate GO term associated with the gene key, and contains the following:
    0) Database ("UniProtKB")
    1) Protein ID
    2) Common gene name
    3) Qualifier (e.g. 'NOT', 'contributes_to', 'colocalizes_with', etc. This is optional and is sometimes an empty string)
    4) GO_ID (i.e. the GO term associated with the gene. NOTE: multiple GO terms associated with a single gene will be in separate entries in the dictionary)
    5) Literature reference that the annotation was derived from
    6) Evidence code
    7) Mystery - don't really know what this column represents and doesn't seem to match README file description on GO website
    8) Aspect (i.e. molecular function, biological process, or cellular component...indicated with abbreviation)
    9) Object symbol (e.g. 'Mitochondrial 21s RNA)
    10) Object synonym
    11) Object type (i.e. gene, protein, etc.)
    12) Taxon ID (559292 for s cerevisiae)
    13) Date of annotation
    14) Assigned by (a database)
    """

    h = open(os.path.join('GOAfiles', 'HQ_GOanalysis', gaf_file))
    
    all_gos = []
    assoc = {}
    
    #creates 'pop' list and 'assoc' dictionary to pass to GOEnrichmentStudy() module in goatools
    for line in h:
        if line.startswith('!'):
            continue
        items = line.rstrip().split('\t')
        temps = items[10].split('|')
        uniprot_id = items[1]
        assoc[uniprot_id] = assoc.get(uniprot_id, [])
        assoc[uniprot_id].append(items[4])
        all_gos.append(items[4])
        
    set_gos = set(all_gos)

    for key in assoc:
        assoc[key] = set(assoc[key])
    
    h.close()

    return assoc
    
    
if __name__ == '__main__':
    main()