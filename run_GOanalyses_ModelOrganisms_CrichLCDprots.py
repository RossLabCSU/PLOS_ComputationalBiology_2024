
from goatools import obo_parser
from goatools.go_enrichment import GOEnrichmentStudy
from Bio import SeqIO
from scipy import stats

all_proteomes = ['UP000006548_3702', 'UP000009136_9913', 'UP000001940_6239', 'UP000002254_9615', 'UP000000437_7955', 'UP000002195_44689', 'UP000000803_7227', 'UP000000539_9031', 'UP000005640_9606', 'UP000000589_10090', 'UP000002311_559292', 'UP000008227_9823', 'UP000002494_10116']
all_abbrevs = ['Athaliana', 'Btaurus', 'Celegans', 'Clupusfamiliaris', 'Drerio', 'Ddiscoideum', 'Dmelanogaster', 'Ggallusdomesticus', 'Hsapiens', 'Mmusculus', 'Scerevisiae', 'Sscrofadomesticus', 'Rnorvegicus']
proteomes_to_abbrevs = {all_proteomes[i]:all_abbrevs[i] for i in range(len(all_proteomes))}
gaf_df = {'UP000006548_3702':'217288.A_thaliana.goa', 'UP000009136_9913':'273325.B_taurus_Hereford_updated.goa', 'UP000001940_6239':'9.C_elegans.goa', 'UP000002254_9615':'61890.C_familiaris.goa', 'UP000000437_7955':'4832498.D_rerio_1.goa', 'UP000002195_44689':'21395.D_discoideum.goa', 'UP000000803_7227':'17.D_melanogaster.goa', 'UP000000539_9031':'4825697.G_gallus_2.goa', 'UP000005640_9606':'25.H_sapiens.goa', 'UP000000589_10090':'59.M_musculus.goa', 'UP000002311_559292':'71242.S_cerevisiae_ATCC_204508.goa', 'UP000008227_9823':'35497.S_scrofa.goa', 'UP000002494_10116':'122.R_norvegicus.goa'}

amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
aa_strings = []
for res1 in amino_acids:
    for res2 in amino_acids:
        if res1 == res2:
            continue
        aa_strings.append(res1+res2)

rare_lcd_classes = ['C' + aa for aa in amino_acids if aa != 'C']
rare_lcd_classes += [aa + 'C' for aa in amino_acids if aa != 'C']
rare_class_indices = [aa_strings.index(lcd_class) for lcd_class in rare_lcd_classes]
        
def main():

    go_obo = 'go-basic.obo'
    go = obo_parser.GODag(go_obo)

    num_prots_threshold = 10
    win_size = 20

    h = open('ModelOrganisms_CrichLCDproteins.tsv')
    header = h.readline()
    for line in h:
        proteome, taxid, domain, sci_name, common_name, total_prots, *data = line.rstrip().split('\t')
        class_freqs = data[:len(rare_lcd_classes)]
        class_prots = data[len(rare_lcd_classes):]
        
        index = all_proteomes.index(proteome)
        abbrev = all_abbrevs[index]
        gaf_file = gaf_df[proteome]
        assoc = get_assoc(gaf_file)
        background_prots = get_background_prots(proteome, win_size)
        
        for i, lcd_class in enumerate(rare_lcd_classes):
            freq = int(class_freqs[i])
            if freq >= num_prots_threshold:
                lcd_prots = class_prots[i].split('_')
                methods = ["bonferroni", "sidak", "holm", "fdr"]
                g = GOEnrichmentStudy(background_prots, assoc, go, propagate_counts=True, alpha=0.05)
                g_res = g.run_study(lcd_prots)
                filename = proteome + '_' + abbrev + '_' + lcd_class + '_GOterm_RESULTS' #change file type here to txt or xlsx if desired
                output_GOresults(g_res, filename)
                
    h.close()
    
    
def output_GOresults(g_res, filename):

    output = open(filename + '.tsv', 'w')
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
        new_line[-1] = new_line[-1].strip()

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
    h = open('./Eukaryota/' + proteome + '.fasta')
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

    h = open(gaf_file)
    
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