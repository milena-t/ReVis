import parse_gff as gff
import genes_windows as gen_windows
import parse_repeats as repeats
from tqdm import tqdm
from Bio import SeqIO



def get_repeat_abundance(repeats_list:list, window_interval:list[int], filter_overlap_previous__, gene_number = False, gene_density=False, calc_gene_density = False) -> dict[str, int]:
    """
    get a dictionary of repeat abundances for the specified interval in the repeats_lists dictionary.
    repeats_lists is the list of repeats in a single contig! it cannot be a dictionary that contains multiple contigs. 
    --> !!the number used for the abundance is basepairs that are covered by a particular repeat category, not the number of times a repeat category is annotated!!\n
    \n
    The overlap can be filtered by cutting the beginning of the current repeat based on the end of the previous one:\n

    --------------------        ----------------                                     -------------              (completely inserted repeats vanish this way)\n
                *************************************               *************************                   (but the overall bases covered by repeats stays similar)\n
                                                                                                                (repeats completely inserted into other repeats that span a wider region are filtered out completely)\n
    --------------------*****************************               *************************-----           \n


    """
    # print(f" \t\t\t >>> filter setting {filter_overlap_previous__}")
    repeat_abundances:dict[str, int] = {}
    window_start, window_stop = window_interval
    # print(f"\t{window_interval} : {len(repeats_list)}")
    # print(f"\t{repeats_list[0]}")
    if len(repeats_list)>0:
        
        # when you have gene_number=True, then you don't care about the nuber of basepairs covered in repeats 
        # but about how many genes fall into the window interval
        # so here I will count the number of features that are placed in the window interval
        no_placed = 0
        overlap_filtered_bp = 0
        prev_end = window_start
        # repeat_cat_name = ""
        if gene_number:
            repeat_cat_name = "Gene"
        if gene_density:
            gene_number=False
            repeat_cat_name = "cds"

        for next_repeat_index, repeat in enumerate(repeats_list, start = 1):


            no_placement = True # if the feature is not placed within the window --> change to True if it is
            covered_bp = 0

            # is the repeat in the window?
            if repeat.end < window_start :
                # repeat is before the window interval continue to the next one
                continue
            elif repeat.start > window_stop:
                # if the repeat starts after the window start then exit
                # because the gff file is sorted by position, so if this occurs there will not be another repeat that is inside the window
                break
            
            # repeat is completely inside the window
            if repeat.start > window_start and repeat.end < window_stop:
            
                # if the repeat does not overlap with the previous one
                if filter_overlap_previous__ and repeat.start > prev_end:
                    covered_bp = repeat.end - repeat.start
                    prev_end = repeat.end
                # if the repeat is completely inside the previous one (the current repeat always starts after the previous one so no need to test that)
                elif filter_overlap_previous__ and repeat.end < prev_end:
                    overlap_filtered_bp += repeat.end - repeat.start
                # if the repeat ends after the previous one but the start is still in the previous repeat
                elif filter_overlap_previous__ and repeat.start < prev_end:
                    covered_bp = repeat.end - prev_end
                    prev_end = repeat.end
                    overlap_filtered_bp += prev_end - repeat.start

                elif filter_overlap_previous__ == False:
                    covered_bp = repeat.end - repeat.start

                no_placement = False
                no_placed +=1
            elif repeat.start < window_start and repeat.end > window_start:
                # repeat overlaps with the start of the interval
                covered_bp = repeat.end - window_start
                # if more than half of the gene is inside the window
                if covered_bp > 0.5*(repeat.end-repeat.start):
                    no_placement = False
                    no_placed +=1
            elif repeat.start > window_start and repeat.end > window_stop:
                # repeat overlaps with the stop of an interval
                covered_bp = window_stop - repeat.start
                # if more than half of the repeat is inside the window
                if covered_bp > 0.5*(repeat.end-repeat.start):
                    no_placement = False
                    no_placed +=1

            if gene_number == False and gene_density ==False:
                repeat_cat_name = repeat.repeat_category
            elif gene_number == False and gene_density==True:
                repeat_cat_name = "cds"
            elif gene_number == True and gene_density == False: 
                repeat_cat_name = "Gene"
                if calc_gene_density:
                    ## this does work with entire gene coordinates, not cds!
                    covered_bp_ratio = float(covered_bp)/float(window_stop - window_start)

            if repeat_cat_name not in repeat_abundances:
                repeat_abundances[repeat_cat_name] = 0

            if no_placement:
                # print(f" --> {no_placed} out of {no_repeats} placed so far, not placed: {repeat}") 
                # break
                ## (maybe) TODO what was the problem with the no_placement stuff? why was i debugging that?
                pass
            
            if gene_number:
                # don't do anything per repeat, only add the number of all placed genes at the end
                continue
            
            # if not gene_number:
            repeat_abundances[repeat_cat_name] += covered_bp
            if covered_bp<0:
                print(f" \t\t !!!! covered_bp = {covered_bp} < 0 -> repeat from {repeat.start} to {repeat.end}")

        if gene_number and not gene_density:
            repeat_abundances[repeat_cat_name] = no_placed

    else:
        # repeat_abundances["Gene"] = 0 
        pass
    
    #print(f"\t\t ----> {window_interval} : overlap filtered bp: {overlap_filtered_bp}, {no_placed} repeat features inside window")
    return(repeat_abundances, overlap_filtered_bp)
        # break


def get_contig_abundances(contig_repeats:list[gff.FeatureCategory.Repeat], windows:list[list[int]], filter_overlap_previous:bool,  verbose=True, gene_number = False, gene_density=False, calc_gene_density = False ) -> dict[int, dict[str, int]]:
    """
    returns a dictionary with { window_number : {repeat abundances} } for contig_repeats
    the repeat_abundances are basepairs covered by repeats, not the number of annotated repeats! (Except if you specify gene_number)
    """
    repeat_abundances:dict[int,dict] = {} # include the integer for the window number to make sure everything stays sorted
    repeat_categories:list[str] = []
    overlap_removed_bp = 0

    for window_num, window_list in enumerate(windows):
        repeats_dict, overlap_filtered_bp = get_repeat_abundance(contig_repeats, window_list, summarize_categories=True, gene_number = gene_number, gene_density=gene_density, filter_overlap_previous__ = filter_overlap_previous, calc_gene_density = calc_gene_density)
        overlap_removed_bp += overlap_filtered_bp
        # print(f" --> test list: {window_list}, {len(repeats_dict)}") # TODO related to the no_placement stuff above
        repeat_abundances[window_num] = repeats_dict
        repeat_categories.extend(list(repeats_dict.keys()))
        # if verbose:
        #     print(f"Window {window_num}: {window_list} ;  repeat abundances: {repeats_dict}")

    unique_repeat_categories = list(set(repeat_categories))
    if verbose and not gene_number:
        print(f"{len(unique_repeat_categories)} unique repeat categories: {unique_repeat_categories}")
    return(repeat_abundances, unique_repeat_categories, overlap_removed_bp)


def get_assembly_repeat_abundances(out_filepath:str, window_length:int, filter_overlap_previous_, statistics_outfile_name, gff_filepath = "", verbose = False, statistics = False, ) -> dict[str:dict[int, dict[str,list[int]]]]:
    """ 
    filter_overlap_previous should only be true if the result is plotted! If you want the table to do statistics then all repeat features regardless of overlap should be kept in their entirety

    get a data structure that contains all the repeat abundances separated by contig for an assembly
    { contig : { 
        window_1 : { 
            repeat_category1 : abundance,
            repeat_category2 : abundance
            },
        window_2 : { 
            repeat_category1 : abundance,
            repeat_category2 : abundance
            }
        } 
    }

    statistics will show individual statistics for each contig instead of a progress bar
    """ 

    min_contig_length = window_length
    
    if verbose:
        print(f"\nset window length: {int(window_length)} bp")
        assembly_length_incl = 0
        repeat_length = 0
        overlap_removed_bp_total = 0
        
    if len(gff_filepath)== 0:
        gff_filepath = out_filepath+".gff"

    
    contig_lengths = gen_windows.get_contig_lengths(gff_filepath)
    num_contigs = len(contig_lengths)
    assembly_length = sum([lengths[1] for lengths in contig_lengths.values()])

    contig_lengths = {contig : length for contig, length in contig_lengths.items() if length[1]>min_contig_length}
    incl_contigs = len(contig_lengths)
    incl_length = sum([lengths[1] for lengths in contig_lengths.values()])
    

    if len(contig_lengths) == 0:
        contig_lengths = gen_windows.get_contig_lengths(gff_filepath)
        longest_contig = max([lengths[1] for lengths in contig_lengths.values()])
        raise RuntimeError(f"All contigs are shorter than the window length of {int(window_length)}. The longest contig is {longest_contig} bp long, you have to pick a threshold shorter than this.")


    if verbose:
        excl_contigs = num_contigs-incl_contigs
        incl_perc = (incl_length/assembly_length) * 100
        print(f"{excl_contigs} contigs are shorter than the window length and therefore excluded. ({incl_perc:.2f} % of assembly included)")
        print(f"\n \t -  parse gff: ")
    
    
    repeat_gff = repeats.parse_repeats_repeatmasker_outfile(out_filepath, verbose = verbose, contig_list = list(contig_lengths.keys()))
    
    if verbose: 
        print(f"\n")

    all_categories:list[str] = []
    repeat_abundances = {} 

    num_contigs = len(contig_lengths)
    if statistics:
        with open(statistics_outfile_name, "a") as stats_out:
            stats_out.write(f"\n{out_filepath}\n")

            for num_curr, contig in enumerate(contig_lengths.keys()):
                assembly_length_incl += contig_lengths[contig][1]

                windows_contig = gen_windows.get_window_intervals(contig_lengths[contig], window_length)
                contig_repeats = repeat_gff[contig]
                
                contig_abundances, contig_categories, overlap_removed_bp = get_contig_abundances(contig_repeats, windows_contig, verbose=False, filter_overlap_previous = filter_overlap_previous_)
                
                # print(contig_abundances)
                if verbose:
                    contig_repeat_bp = int(sum([sum(list(window_abundances.values())) for window_abundances in contig_abundances.values()]))
                    repeat_length += contig_repeat_bp
                
                all_categories.extend(contig_categories)
                repeat_abundances[contig] = contig_abundances
                if verbose:
                    overlap_removed_bp_total += overlap_removed_bp
                    print(f"\n =============>  ( contig {num_curr+1}/{num_contigs} )  {contig} :   {contig_lengths[contig][1]} bp ")
                    stats_out.write(f"( contig {num_curr+1}/{num_contigs} )  {contig} :   {contig_lengths[contig][1]} bp\n")

                    print(f"\tthere are {len(contig_repeats)} annotated repeat sequences (the number of separately annotated features, not bp) of {len(contig_categories)} unique categories in this contig")
                    stats_out.write(f"\tthere are {len(contig_repeats)} annotated repeat sequences (the number of separately annotated features, not bp) of {len(contig_categories)} unique categories in this contig\n")
                    contig_repeat_percentage = (contig_repeat_bp/contig_lengths[contig][1]) *100
                    print(f"\t{contig_repeat_percentage:.2f}% repeats ( {contig_repeat_bp} / {contig_lengths[contig][1]} ) bp, {overlap_removed_bp} repeat-annotated bp removed because they were covered by overlapping features (filter setting: {filter_overlap_previous_})")
                    stats_out.write(f"\t{contig_repeat_percentage:.2f}% repeats ( {contig_repeat_bp} / {contig_lengths[contig][1]} ) bp, {overlap_removed_bp} repeat-annotated bp removed because they were covered by overlapping features (filter setting: {filter_overlap_previous_})\n")
                    print(f"\t{contig_categories}")
                    stats_out.write(f"\t{contig_categories}\n")
            stats_out.write(f"\n")
    
    else: # show progress bar if not verbose
        print(f" \t -  creating windows for parsed contigs\n")
        # print(contig_lengths)

        for contig in tqdm(contig_lengths.keys()):
            assembly_length_incl += contig_lengths[contig][1]

            windows_contig = gen_windows.get_window_intervals(contig_lengths[contig], window_length)
            contig_repeats = repeat_gff[contig]
            
            contig_abundances, contig_categories, overlap_removed_bp = get_contig_abundances(contig_repeats, windows_contig, verbose=False, filter_overlap_previous = filter_overlap_previous_)
            
            # print(contig_abundances)
            if verbose:
                overlap_removed_bp_total += overlap_removed_bp
                contig_repeat_bp = int(sum([sum(list(window_abundances.values())) for window_abundances in contig_abundances.values()]))
                repeat_length += contig_repeat_bp
            
            all_categories.extend(contig_categories)
            repeat_abundances[contig] = contig_abundances

    
    all_categories = list(set(all_categories))

    if verbose:
        print(f"\n\tIn total, there are {len(all_categories)} unique repeat categories: {all_categories}")
        repeat_percentage = (repeat_length/assembly_length_incl) * 100
        print(f"\t{assembly_length_incl} bp of sequence were parsed, of which {repeat_length} bp were annotated as repeats ({repeat_percentage:.2f}%)")
        overlap_percentage = (overlap_removed_bp_total/assembly_length_incl) * 100
        print(f"\t{overlap_removed_bp_total} bp were removed from the total repeat coverage because they were covered by more than one repeat feature ( {overlap_percentage:.2f} % of the assembly) (filter setting: {filter_overlap_previous_})")
    
    return repeat_abundances, all_categories

def get_masked_bp_in_window(assembly_path:str, window_length:int, verbose = False) -> dict[str,list[list[int,float]]]:
    """
    Get the number of unmasked bases in windows across the entire assembly.
    I am counting unmasked bases, uppercase A G C T, so both lowercase softmasked bases and hardmasked bases (N) are counted as masked (not-unmasked)

    {
        contig1 : {
            window_index1 : unmasked_proportion,
            window_index2 : unmasked_proportion,
            ...
        } ,
        contig2 : {
            window_index1 : unmasked_proportion,
            window_index2 : unmasked_proportion,
            ...
        }
    }
    """

    assembly = SeqIO.parse(assembly_path, "fasta")
    unmasked_bp:dict = {}
    unmasked_bases = ["A", "C", "G", "T"]

    ## get contig lengths dict { contig1 : [ 1 , cont_length] , contig2 : [ 1 , cont_length] ,  ... }
    # contig_lengths = { record.id : [1 , len(record.seq)] for record in assembly if len(record.seq) > window_length}
    # print(contig_lengths)
    if verbose:
        print(f"\n * comute unmasked bp per contig from assembly")
    for record in tqdm(assembly):
        if len(record.seq)>window_length:
            contig_length_list = [1,len(record.seq)]
        else:
            continue
        window_intervals = gen_windows.get_window_intervals(contig_coords=contig_length_list, window_length=window_length)

        unmasked_bp[record.id] = {}

        # print()
        # print(record.id)
        contig_seq = str(record.seq)
        for window_ind , window_interval in enumerate(window_intervals):
            # window_seq = Seq(contig_seq[int(window_interval[0])-1:int(window_interval[1])])
            window_seq = contig_seq[int(window_interval[0])-1:int(window_interval[1])]

            unmasked_bases_count = 0
            for base in unmasked_bases:
                unmasked_bases_count += window_seq.count(base)

            unmasked_bp[record.id][window_ind] = unmasked_bases_count

                #print(record.id, window_interval,  len(window_seq), unmasked_bases_count, masked_bases)

    return unmasked_bp
