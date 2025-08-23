#####
# Analysis and plotting of repeatmasker output, see documentation below or repeatmasker_window_analysis.py -h
# Milena R. Trabert, 2025
#####

import parse_gff as gff
import parse_repeats as repeats
from tqdm import tqdm
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.ticker import FuncFormatter
import pandas as pd
from Bio import SeqIO
import random
import time
import re
from statistics import mean

####!!
# on uppmax load biopython/1.80-py3.10.8 to use argparse!
import argparse


def parse_args():
    # Create the parser
    program_description = """\n
A script to calculate repeat abundance in non-overlapping windows over an entire assembly from repeatmasker output.
Results can be plotted in a stacked histogram (plot mode) or returned as tsv files (table mode). 


------------------------- quick start

python3 repeatmasker_window_analysis.py --masker_outfile your_assembly.fna.ori.out --masker_out_gff your_assembly.fna.out.gff --species_name Your_species --window_length 1e6 --verbose --plot

------------------------- dependencies

 *  On uppmax, load biopython/1.80-py3.10.8 to use argparse (base python doesn't include it).
 *  The parse_repeats.py file that contains the dataclass for the repeat files, the parse_gff.py file contains other utilities
        (keep in the same directory as this script or in $PATH to make the import work!). 
 *  Libraries imported in this script and in parse_repeats.py and parse_gff.py:
        - sys, re, os, subprocess, argparse
        - tqdm, random, time
        - dataclasses, enum, pandas, tempfile, SeqIO
        - matplotlib
        (they can all be installed through pip. Sorry if I forgot any, I'm sure the compiler will tell you)

------------------------- documentation

This script analyzes repeatmasker output. You can run in two modes, plot mode or table mode. Plot mode returns 
a plot with a stacked histogram of all the repeat categories (overlap filtered! Repeat annotations can be nested within each other, 
where two or more repeats cover the same stretch of sequence, which can result in a >100% repeat coverage in some windows. 
I am filtering out repeat categories that overlap with others so that each base is only covered by one repeat annotation.
This filtering likely warps the category proportion in some windows, because it is likely that there are more bases of some categories 
removed than others. See get_repeat_abundance function documentation for details). If you include a genome annotation, 
a line for the number of annotated genes in the same windows is added. You can add up to two annotations to compare them.
Table mode returns the same by-window information as a tsv file with the proportion of basepairs in each window 
covered by each repeat category (NOT overlap filtered, so there can be more bp covered by all repeats in a window than the number
of masked bp or even the window length) and also the number of bp and ratio covered by coding regions (exons). If the masked
assembly is given it will also include the number of unmasked bp in each window.

It takes as input two of the repeatmasker output files, and since I didn't exactly test this rigorously, 
it would probably be good if they were in the same directory and didn't have their names changed from what 
repeatmasker named them by default (in case i have a hardcoded string in a filename somewhere i forgot to remove).

The runtime depends on the overall repeat content and on how fragmented the assembly is. I have tried my best to optimize, 
but if your assembly is long and fragmented, it takes long to loop through many small contigs, and if there are many repeats, 
it takes long to sum all of them up per window. Short windows increase the total number of windows computed and plotted, 
which also increases the runtime. In any case, the longest runtime i managed to achieve with my data was 3:30 min.

------------------------- good luck! 

"""
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,description=program_description)

    # Add the arguments
    parser.add_argument('--masker_outfile', type=str, required=True, help="""repeatmasker output file ending in .out 
    (I really recommend .ori.out, but both work, the other one is just slower)""")
    parser.add_argument('--masker_out_gff', type=str, help="""repeatmasker output file ending in .out.gff
    If not given it will be assumed to have the same basename as .out and inferred automatically.
    (which is true if you didn't change anything about the repeatmasker outfile names.) 
    (This technically contains the same information as the .out file, but it also has "sequence region" annotations
    which have all the contig lengths. These would otherwise have to be extracted from the assembly and this is just faster.)""")

    # output format: mutually exclusive group of arguments, either table or plot
    output = parser.add_mutually_exclusive_group(required=True)
    output.add_argument("--table", action="store_true", help = """The output of the program is the tsv file with all values""")
    output.add_argument("--plot", action="store_true", help = """The output of the program is the plot""")
    
    parser.add_argument('--annotation_gff', type=str, help="""If you want to include a line that shows the gene density in each window, 
    you can give an annotation based on the same assembly (with matching contig names!)""")
    parser.add_argument('--annotation_gff2', type=str, help="""You can include a second annotation as well if you want to and show to different lines.""")
    parser.add_argument('--gene_density', action="store_true", help="""If you give an annotation, by default you will get the number of genes per window, 
    but you may also calculate the ratio: number of genes / window length with this option""")
    parser.add_argument('--assembly_path', type=str, help="""If you want to compute the table with all the per-window values you need to include the
    masked assembly, so that the number of masked/unmasked bases can be computed properly. It should be the assembly that was masked in this repeatmasker
    run and also the assembly that is the basis for the annotation""")

    parser.add_argument('--species_name', type=str, help="""species identifier string, like 'C_maculatus' 
    (Include this! If not included it will try to parse it automatically from the start of filenames, 
    which will probably only work for how i named my files)""")
    parser.add_argument('--window_length', type=float, required = True, help="window length (scientific notation like 1e6 is ok)")
    parser.add_argument('--merge_gene_windows', type=int, help="""If an annotation is included to show gene density, 
    the gene density is likely difficult to interpret visually if you show it for each repeat-window.
    Here you can choose how many windows you would like to average over, 
    default = 5
    put 1 if you don't want any averaging""")

    parser.add_argument('--verbose', action="store_true", help="print progress in the command line (recommended, on by default)")
    parser.add_argument('--statistics', action="store_true", help="""print contig-specific repeat statistics (optional, increases output length quite a bit, off by default)
    Prints repeat information (and gene numbers if applicable) for each contig, mostly for debugging purposes""")




    # Parse the arguments
    args = parser.parse_args()

    # make defaults
    if not args.verbose:
        args.verbose = True
    if not args.statistics:
        args.statistics = False
    if not args.masker_out_gff:
        args.masker_out_gff = args.masker_outfile + ".gff"
    if not args.species_name:
        args.species_name = ""
    if args.species_name:
        pattern = r'\s+'
        args.species_name = re.sub(pattern, '_', args.species_name)
    if not args.gene_density:
        args.gene_density = False
    if not args.merge_gene_windows:
        args.merge_gene_windows = 5

    return args


def random_hex_color():
    """
    custom function to format random hex color code
    """
    return "#{:06x}".format(random.randint(0, 0xFFFFFF))


def get_contig_lengths(gff_path:str) -> dict[str, list[int]]:
    """
    the gff generated by repeatmasker contains ##sequence-region comments with the contig name and length
    extract them into a dictionary with { name:str : [start:int , end:int]}
    (start:int is 1 and not 0! the gff file does not 0-index the sequence region features)
    """
    out_dict: dict[str, list[int]] = {}
    with open(gff_path, "r") as gff_file:
        seq_region_lines = list(line for line in gff_file if "sequence-region" in line)
    for seq_region in seq_region_lines:
        prefix, contig, start, end = seq_region.strip().split()
        out_dict[contig] =  [int(start), int(end)]
    return(out_dict)

def get_window_intervals(contig_coords:list[int], window_length:int) -> list[list[int]]:
    """
    get the intervals of the non-overlapping windows in a list of lists like this:
    [ [window1_start, window1_end] ,  [window2_start, window2_end] ,  ... ]
    """
    contig_end = contig_coords[1]
    contig_start = contig_coords[0]

    windows_list: list[list[int]] = []
    curr_window_start = contig_start
    curr_window_end = window_length

    while curr_window_end < contig_end:
        windows_list.append([curr_window_start, curr_window_end])
        curr_window_start = curr_window_end+1
        curr_window_end = curr_window_start+window_length-1
    
    windows_list.append([curr_window_start, contig_end])

    return(windows_list)

def get_repeat_abundance(repeats_list:list, window_interval:list[int], filter_overlap_previous__, gene_number = False, summarize_categories = True, gene_density=False, filter_overlap_following = False, calc_gene_density = False) -> dict[str, int]:
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


def get_assembly_repeat_abundances(out_filepath:str, window_length:int, filter_overlap_previous_, gff_filepath = "", verbose = False, statistics = False, ) -> dict[str:dict[int, dict[str,list[int]]]]:
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

    
    contig_lengths = get_contig_lengths(gff_filepath)
    num_contigs = len(contig_lengths)
    assembly_length = sum([lengths[1] for lengths in contig_lengths.values()])

    contig_lengths = {contig : length for contig, length in contig_lengths.items() if length[1]>min_contig_length}
    incl_contigs = len(contig_lengths)
    incl_length = sum([lengths[1] for lengths in contig_lengths.values()])
    

    if len(contig_lengths) == 0:
        contig_lengths = get_contig_lengths(gff_filepath)
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

                windows_contig = get_window_intervals(contig_lengths[contig], window_length)
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

            windows_contig = get_window_intervals(contig_lengths[contig], window_length)
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


def plot_repeat_abundance(species_abundances, species_categories, gff_filepath, window_length, gene_numbers = {}, gene_numbers2 = {}, species_name = "", y_label = "repeat abundance", filename = ""):

    print("\n* start plotting...")
    
    if len(species_name)==0:
        species_name = gff.split_at_second_occurrence(gff_filepath.split("/")[-1])
        #species_name = species_name.replace("_", ". ")

    if len(filename)==0:
        filename = f"repeat_abundance_in_{species_name}.png"
        if len(gene_numbers)>0:
            filename = f"repeat_abundance_with_gene_numbers_in_{species_name}.png"

    # Width of the bars
    # width = 0.35
    fs = 25 # fontsize is scaled with the dpi somehow which i have to do extra because i change the aspect ratio manually below

    # set figure aspect ratio
    aspect_ratio = 20 / 12
    height_pixels = 1200  # Height in pixels
    width_pixels = int(height_pixels * aspect_ratio)  # Width in pixels

    fig, ax = plt.subplots(figsize=(width_pixels / 100, height_pixels / 100), dpi=100)
    # fig = plt.figure(figsize=(width_pixels / 100, height_pixels / 100), dpi=100)
    # ax = fig.add_subplot(111)
    
    # because everyting is dictionaries, it is unordered and the barplot stacks and colors are different in every iteration
    # so i order a list to sort it
    categories_inferred = sorted(species_categories)
    categories_default = ['Unknown', 
                            'DNA', 
                            'tRNA', 
                            'rRNA',
                            #'snRNA', 
                            'LTR', 
                            'RC', 
                            'LINE', 
                            'SINE', 
                            'Low_complexity',
                            'Satellite', 
                            'Simple_repeat']

    if set(categories_default) == set(categories_inferred):
        sorted_categories = categories_default
    elif "snRNA" in categories_inferred:
        cat_removed = [cat for cat in categories_inferred if cat != "snRNA"]
        if set(cat_removed) == set(categories_default):
            sorted_categories = ['Unknown', 
                            'DNA', 
                            'snRNA', 
                            'tRNA', 
                            'rRNA',
                            'LTR', 
                            'RC', 
                            'LINE', 
                            'SINE', 
                            'Low_complexity',
                            'Satellite', 
                            'Simple_repeat']
        else:
            sorted_categories = sorted(categories_inferred)
            if 'Unknown' in sorted_categories:
                sorted_categories.remove('Unknown')
                sorted_categories.insert(0,'Unknown')
    else:
        sorted_categories = sorted(categories_inferred)
        if 'Unknown' in sorted_categories:
            sorted_categories.remove('Unknown')
            sorted_categories.insert(0,'Unknown')

        
    sorted_categories.reverse()
    
    contig_lengths_all = get_contig_lengths(gff_filepath)
    num_contigs = len(contig_lengths_all)
    assembly_length = sum([lengths[1] for lengths in contig_lengths_all.values()])

    contig_lengths = {contig : length for contig, length in contig_lengths_all.items() if length[1]>window_length}
    incl_contigs = len(contig_lengths)
    incl_length = sum([lengths[1] for lengths in contig_lengths.values()])

    incl_perc = (incl_length/assembly_length) * 100


    colors = {
        'Unknown' : "#C1C1C1" , # light grey
        # orange
        'DNA' : "#FF9000" , # Princeton orange
        # green
        'LTR' : "#6E8448" , # reseda green
        'RC' : "#8EA861" , # asparagus 
        # red
        'tRNA' : "#C14953" , # bittersweet shimmer
        'rRNA' : "#D0767E" , # old rose
        'snRNA' : "#7A2A30" , # wine
        # blue 
        'LINE' : "#3476AD" , # UCLA blue
        'SINE': "#72A8D5" , # ruddy blue
        # '' : "#2A618D" , #lapis lazuli
        # dark red-brown
        'Low_complexity' : "#3A3335" , # Jet 
        'Satellite' : "#564D4F" , #Wenge 
        'Simple_repeat' : "#827376" , #Taupe gray
    }

    genes_line_color = "#890000"
    genes_line_color2 = "#262CA6" 
    ## if there are repeat categories that are not representd in the colors, identify them and assign random colors. 

    not_colored_categories = [rep_cat for rep_cat in categories_inferred if rep_cat not in categories_default]

    # print(not_colored_categories)
    if len(not_colored_categories)>0:
        for not_colored_category in not_colored_categories:
            colors[not_colored_category] = random_hex_color()



    ##########################
    ######## PLOTTING ########
    ##########################

    curr_contig_start = 0
    x_contig_coords = []
    x_contig_labels = []

    ax2 = ax.twinx() # axis that plots the contig labels
    ax2.spines['right'].set_visible(False)

    ax2.yaxis.set_ticks_position('none')  # Hide the default ticks
    ax2.yaxis.set_ticklabels([])  # Hide the default tick labels 
    
    include_genes_line = False
    include_genes_line2 = False
    if len(gene_numbers)>0:
        include_genes_line = True

        max_genes = get_max_genes_in_window(gene_numbers)
        y_max_genes = max_genes * 2

        ax3 = ax.twinx()
        ax3.set_ylim(bottom=0, top=y_max_genes)
        ax3.set_ylabel('Number of genes per window', color = genes_line_color, fontsize = fs)
        ax3.tick_params(axis ='y', labelcolor = genes_line_color, labelsize = fs) 
    
        if len(gene_numbers2)>0:
            include_genes_line2 = True

            max_genes2 = get_max_genes_in_window(gene_numbers)
            y_max_genes2 = max_genes2 * 2

            if y_max_genes<y_max_genes2:
                ax3.set_ylim(bottom=0, top=y_max_genes2)


    for contig in contig_lengths.keys():
        window_intervals = get_window_intervals(contig_lengths[contig], window_length)
        window_centers = [window[0]+ 0.5*window_length for window in window_intervals]
        window_abundances = species_abundances[contig]

        for w in range(len(window_abundances)):
            bottom = 0
            curr_window_length = window_intervals[w][1] - window_intervals[w][0]
            curr_window_pos = window_centers[w] - (window_length - curr_window_length) * 0.5
            for category in sorted_categories:
                try: # if there's a key error because there are no repeats of a category set that category to 0
                    rep_proportion = window_abundances[w][category]/curr_window_length
                except:
                    rep_proportion = 0
                rep_percentage = rep_proportion*100

                ax.bar(curr_window_pos+curr_contig_start, rep_percentage, width=curr_window_length, label=category, bottom=bottom, color = colors[category])

                bottom += rep_percentage
            

        if include_genes_line:

            window_intervals_genes = get_window_intervals(contig_lengths[contig], window_length)
            window_centers_genes = [window[0]+ 0.5*window_length for window in window_intervals]

            genes_window = gene_numbers[contig]
            if len(genes_window)<1:
                print(gene_numbers)
                print(contig)
                raise RuntimeError("AAAAAAAA") # TODO fix message

            gene_positions:list[int] = []
            gene_numbers_in_window:list[int] = []

            for w in range(len(genes_window)):
                curr_window_length = window_intervals_genes[w][1] - window_intervals_genes[w][0]
                try:
                    gene_positions.append(window_centers_genes[w] - (window_length - curr_window_length) * 0.5 + curr_contig_start)
                except:
                    #pass
                    print(w)
                    print(contig)
                    print(genes_window)
                    print(contig_lengths_all[w])
                    gene_positions.append(contig_lengths_all[w][1] * 0.5 + curr_contig_start)   ########
                try:
                    gene_numbers_in_window.append(int(genes_window[w]["Gene"]))
                except: 
                    # print(genes_window[w])
                    # print(w)
                    # print(contig)
                    # --> Some last windows are shorter, and then there's no genes any more
                    # TODO maybe NaN? 
                    gene_numbers_in_window.append(0)
                    # raise RuntimeError("AAAAAAAA")

            if include_genes_line2:

                gene_positions2:list[int] = []
                gene_numbers_in_window2:list[int] = []
                genes_window2 = gene_numbers2[contig]

                # print(contig)

                for w in range(len(genes_window2)):
                    curr_window_length = window_intervals_genes[w][1] - window_intervals_genes[w][0]
                    try:
                        gene_positions2.append(window_centers_genes[w] - (window_length - curr_window_length) * 0.5 + curr_contig_start)
                    except:
                        #pass
                        print(w)
                        print(contig)
                        print(genes_window2)
                        print(contig_lengths_all[w])
                        gene_positions2.append(contig_lengths_all[w][1] * 0.5 + curr_contig_start)   ########
                    try:
                        gene_numbers_in_window2.append(int(genes_window2[w]["Gene"]))
                    except: 
                        # print(genes_window2[w])
                        # print(w)
                        # print(contig)
                        # --> Some last windows are shorter, and then there's no genes any more so the ["Gene"] value doesn't exist
                        gene_numbers_in_window2.append(0)
                        # raise RuntimeError("AAAAAAAA")

                    
                    # gene_orthoDB = genes_window2[w]["Gene"]
                    # gene_native = genes_window[w]["Gene"] 
                    # print(f"\tnative: {gene_native} -- {gene_orthoDB} :orthoDB")
            
            ax3.set_ylabel('Number of genes per window', color = genes_line_color, fontsize = fs)
            ax3.tick_params(axis ='y', labelcolor = genes_line_color, labelsize = fs)  
            
            gene_nos_plot, = ax3.plot(gene_positions, gene_numbers_in_window, color = genes_line_color, label = "native")
            # gene_nos_plot, = ax3.plot(gene_positions, gene_numbers_in_window, color = genes_line_color2, label = "orthoDB", linestyle = ":")

            if include_genes_line2:
                gene_nos_plot2, = ax3.plot(gene_positions2, gene_numbers_in_window2, color = genes_line_color2, linestyle = ":", label = "orthoDB")

        # x_contig_coords.append(num_curr*contig_lengths[contig][1] + contig_lengths[contig][1]*0.5)
        x_contig_coords.append(curr_contig_start + contig_lengths[contig][1]*0.5)
        x_contig_labels.append(contig)
        curr_contig_start += window_length+contig_lengths[contig][1]
        # break

    # convert x axis to Mb
    # Set the primary x-axis formatter to convert bp to Mb
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: '' if x < 0 else f'{x / 1e6:.0f} Mb'))

    # add contig names
    # Set the secondary x-axis ticks and labels
    ax2 = ax.secondary_xaxis('bottom')
    # print(x_contig_coords)
    # print(x_contig_labels)
    ax2.set_xticks(x_contig_coords)
    rot = 90
    if len(x_contig_labels[0])<2:
        rot = 0
    ax2.set_xticklabels(x_contig_labels, rotation=rot, fontsize=fs*0.7)

    # Adjust the position of the secondary x-axis
    ax2.spines['bottom'].set_position(('outward', 40))    
    ax2.xaxis.set_ticks_position('none')
    ax2.spines['bottom'].set_visible(False)

    # change tick fontsizes
    ax.tick_params(axis='x', labelsize=fs) 
    ax.tick_params(axis='y', labelsize=fs)
    ax2.tick_params(axis='x', labelsize=fs*0.8, rotation = rot)
    
    # ax.set_xlabel('Species', fontsize=fs+4)
    ax.set_ylabel(y_label, fontsize=fs)
    ax.yaxis.set_major_formatter(FuncFormatter(lambda x, pos: '' if x > 100 and x<5 else f'{int(x)}%'))
    species_name_nice = species_name.replace("_", ". ")

    if window_length>=1e6:
        rounded = int(window_length / 1e6)
        window_length_ =f'{rounded} Mb'
    elif window_length>=1e5:
        rounded = int(window_length / 1e5) * 0.1
        window_length_ =f'{rounded:.1} Mb'
    else:
        rounded = int(window_length / 1e3)
        window_length_ =f'{rounded} kb'
    
    ax.set_title(f'Repeat abundance in {species_name_nice}, (window length {window_length_}, {incl_perc:.2f}% of assembly shown)', fontsize=fs)

    if include_genes_line and include_genes_line2:
        # make legends for repeat lines (hard coded labes in plot command above!)
        # repeat_lines = [gene_nos_plot, gene_nos_plot2]
        repeat_lines = ax3.get_lines()
        legend_lines = plt.legend(repeat_lines, ["native", "uniform"], loc = "upper left", fontsize = fs*0.75, title = "annotation", title_fontsize = fs*0.85)
        plt.gca().add_artist(legend_lines)
    elif include_genes_line and not include_genes_line2:
        repeat_lines = ax3.get_lines()
        legend_lines = plt.legend(repeat_lines, ["native"], loc = "upper left", fontsize = fs*0.75, title = "annotation", title_fontsize = fs*0.85)
        # legend_lines = plt.legend(repeat_lines, ["orthoDB"], loc = "upper left", fontsize = fs, title = "annotation", title_fontsize = fs)
        plt.gca().add_artist(legend_lines)

    # for Repeat barplots: make legend patches and labels
    handles = []
    labels = [] 
    for category in  sorted_categories[::-1]:# species_abundances[list(contig_lengths.keys())[0]][1]:
        handles.append(mpatches.Patch(color=colors[category]))
        labels.append(category)

    ax.legend(handles, labels, fontsize = fs*0.75, loc='lower left', title = "repeat categories", title_fontsize = fs*0.85)

    # make space for legend on the left
    x_limits = plt.xlim()
    new_lower_limit = x_limits[0] - 0.15 * (x_limits[1] - x_limits[0])
    if len(gene_numbers)>0:
        new_lower_limit = x_limits[0] - 0.17 * (x_limits[1] - x_limits[0])
    ax.set_xlim(new_lower_limit, x_limits[1])

    ax.set_ylim(0, 105)

    plt.tight_layout()

    plt.savefig(filename, dpi = 300, transparent = True)
    print("\tFigure saved in the current working directory directory as: "+filename)

    # plt.show()


def get_assembly_gene_numbers(annotation_filepath, gff_filepath, window_length, verbose = True, statistics = False, calc_gene_density = False):

    """ 
    get a data structure that contains all the gene numbers separated by contig for an assembly
    { contig : { 
        window_1 : { 
            "Gene" : number of anntoated genes (!! total number of annotated genes, not bp covered)
            },
        window_2 : { 
            "Gene" : number of anntoated genes
            }
        } 
    }

    Use merge_gene_windows to get the average number of anntoated genes over several windows to increase readability of the graph later.
    This then means that e.g. window_1 to window_5 have the same value, which is the average over all the previously different values of these windows

    Use get_assembly_gene_abundances for getting actual bp covered by coding regions in each window
    """ 

    
    if verbose:
        print(f"\ncalculate gene densities across the windows")
        print(f"set window length for genes: {int(window_length)} bp")
        assembly_length_incl = 0


    # contig_lengths = get_contig_lengths(gff_filepath)
    ## no need to check for len(contig_lengths) > 0 here that already happens for the repeat abundances
    
    contig_lengths = get_contig_lengths(gff_filepath)
    num_contigs = len(contig_lengths)
    assembly_length = sum([lengths[1] for lengths in contig_lengths.values()])

    contig_lengths = {contig : length for contig, length in contig_lengths.items() if length[1]>window_length}
    incl_length = sum([lengths[1] for lengths in contig_lengths.values()])
    
    if verbose:
        incl_perc = (incl_length/assembly_length) * 100
        print(f"{incl_length} bp are in are shorter than the window length and therefore excluded. ({assembly_length} bp in total, {incl_perc:.2f} % of assembly included)")
        print(f"\n \t -  parse the species annotation: ")
    
    annotation = gff.parse_gff3_by_contig(annotation_filepath, verbose=verbose)
    
    if verbose: 
        print(f"\n")

    all_categories:list[str] = [] # should only be ["Gene"] but double check
    gene_abundances = {} 

    num_contigs = len(contig_lengths)
    gene_length = 0

    if statistics:
        with open(statistics_outfile_name, "a") as stats_out:
            stats_out.write(f"\n{annotation_filepath}\n")
            for num_curr, contig in enumerate(contig_lengths.keys()):
                assembly_length_incl += contig_lengths[contig][1]

                windows_contig = get_window_intervals(contig_lengths[contig], window_length)
                contig_genes = annotation[contig]
                
                contig_abundances, contig_categories, overlap_removed_bp = get_contig_abundances(contig_genes, windows_contig, verbose=statistics, gene_number = True, calc_gene_density = calc_gene_density, filter_overlap_previous_=False)
                
                # print(contig_abundances)
                
                contig_gene_bp = int(sum([sum(list(window_abundances.values())) for window_abundances in contig_abundances.values()]))
                gene_length += contig_gene_bp

                gene_numbers = [list(gene_no.values())[0] if len(list(gene_no.values()))>0 else 0 for gene_no in list(contig_abundances.values()) ]
                # print(gene_numbers)
                mean_genes_in_window = mean(gene_numbers)
                num_windows = len(contig_abundances)
                est_mean_genes_in_window = len(contig_genes)/num_windows
                
                all_categories.extend(contig_categories)
                gene_abundances[contig] = contig_abundances
            
                print(f"\t -->  ( contig {num_curr+1}/{num_contigs} ) {contig} :  {len(contig_genes)} annotated genes (the number of separately annotated features, not bp) in {contig_lengths[contig][1]} bp ")
                stats_out.write(f"( contig {num_curr+1}/{num_contigs} ) {contig} :  {len(contig_genes)} annotated genes (the number of separately annotated features, not bp) in {contig_lengths[contig][1]} bp\n")
                print(f"\t      \t( {mean_genes_in_window:.2f} genes on avg. per window are computed. for evenly distributed genes the estimate would be {est_mean_genes_in_window:.2f}")
                stats_out.write(f"\t( {mean_genes_in_window:.2f} genes on avg. per window are computed. for evenly distributed genes the estimate would be {est_mean_genes_in_window:.2f}\n")
                # contig_repeat_percentage = (contig_gene_bp/contig_lengths[contig][1]) *100
                # print(f"\t{contig_repeat_percentage:.2f}% genes ( {contig_gene_bp} / {contig_lengths[contig][1]} ) bp")
                # print(f"\t{contig_categories}")
    
    else: # show progress bar if not verbose
        print(f" \t -  creating windows for parsed contigs\n")

        for contig in tqdm(contig_lengths.keys()):
            assembly_length_incl += contig_lengths[contig][1]

            windows_contig = get_window_intervals(contig_lengths[contig], window_length)
            contig_genes = annotation[contig]
            
            contig_abundances, contig_categories, overlap_removed_bp = get_contig_abundances(contig_genes, windows_contig, verbose=False, gene_number = True, calc_gene_density = calc_gene_density, filter_overlap_previous=False)
            
            # print(contig_abundances)
            if verbose:
                contig_gene_bp = int(sum([sum(list(window_abundances.values())) for window_abundances in contig_abundances.values()]))
                gene_length += contig_gene_bp
            
            all_categories.extend(contig_categories)
            gene_abundances[contig] = contig_abundances
    
    all_categories = list(set(all_categories))

    if verbose:
        print(f"\n\tIn total, there are {len(all_categories)} unique feature categories: {all_categories} (should be only one for the genes)")
        cds_percentage = (gene_length/assembly_length_incl) * 100
        # print(f"\t{assembly_length_incl} bp of sequence were parsed, of which {gene_length} bp were annotated as genes ({cds_percentage:.2f}%)")
    
    return gene_abundances, all_categories






def get_assembly_gene_abundances(annotation_filepath, gff_filepath, window_length, verbose = True, statistics = False, calc_gene_density = False):

    """ 
    get a data structure that contains all the gene numbers separated by contig for an assembly
    { contig : { 
        window_1 : { 
            "Gene" : number of coding basepairs (!! not the number of annotated genes)
            },
        window_2 : { 
            "Gene" : number of coding basepairs
            }
        } 
    }

    Use merge_gene_windows to get the average number of anntoated genes over several windows to increase readability of the graph later.
    This then means that e.g. window_1 to window_5 have the same value, which is the average over all the previously different values of these windows
    """ 

    
    if verbose:
        print(f"\ncalculate gene densities across the windows")
        print(f"set window length for genes: {int(window_length)} bp")
        assembly_length_incl = 0


    # contig_lengths = get_contig_lengths(gff_filepath)
    ## no need to check for len(contig_lengths) > 0 here that already happens for the repeat abundances
    
    contig_lengths = get_contig_lengths(gff_filepath)
    num_contigs = len(contig_lengths)
    assembly_length = sum([lengths[1] for lengths in contig_lengths.values()])

    contig_lengths = {contig : length for contig, length in contig_lengths.items() if length[1]>window_length}
    incl_length = sum([lengths[1] for lengths in contig_lengths.values()])
    
    if verbose:
        incl_perc = (incl_length/assembly_length) * 100
        print(f"{incl_length} bp are in are shorter than the window length and therefore excluded. ({assembly_length} bp in total, {incl_perc:.2f} % of assembly included)")
        print(f"\n \t -  parse the species annotation: ")
    
    annotation = gff.parse_gff3_by_contig(annotation_filepath, verbose=verbose, featurecategory = gff.FeatureCategory.Exon)
    
    if verbose: 
        print(f"\n")

    all_categories:list[str] = [] # should only be ["Gene"] but double check
    gene_abundances = {} 

    num_contigs = len(contig_lengths)
    cds_length = 0

    all_windows:dict = {} # { contig : [[start1,end1] , [start2,end2] , ...]}

    if statistics:
        with open(statistics_outfile_name, "a") as stats_out:
            stats_out.write(f"\n{annotation_filepath}\n")
            for num_curr, contig in enumerate(contig_lengths.keys()):
                assembly_length_incl += contig_lengths[contig][1]

                windows_contig = get_window_intervals(contig_lengths[contig], window_length)
                all_windows[contig] = windows_contig

                contig_genes = annotation[contig]
                
                contig_abundances, contig_categories, overlap_removed_bp = get_contig_abundances(contig_genes, windows_contig, verbose=statistics, gene_number = True, gene_density = True, calc_gene_density = calc_gene_density, filter_overlap_previous=False)
                
                contig_gene_bp = int(sum([sum(list(window_abundances.values())) for window_abundances in contig_abundances.values()]))
                cds_length += contig_gene_bp

                gene_numbers = [list(gene_no.values())[0] if len(list(gene_no.values()))>0 else 0 for gene_no in list(contig_abundances.values()) ]
                # print(gene_numbers)
                mean_genes_in_window = mean(gene_numbers)
                num_windows = len(contig_abundances)
                est_mean_genes_in_window = len(contig_genes)/num_windows
                
                all_categories.extend(contig_categories)
                gene_abundances[contig] = contig_abundances
            
                print(f"\t -->  ( contig {num_curr+1}/{num_contigs} ) {contig} :  {len(contig_genes)} annotated cds (the number of separately annotated features, not bp) in {contig_lengths[contig][1]} bp ")
                stats_out.write(f"( contig {num_curr+1}/{num_contigs} ) {contig} :  {len(contig_genes)} annotated cds (the number of separately annotated features, not bp) in {contig_lengths[contig][1]} bp\n")
                cds_percentage = cds_length / contig_lengths[contig][1]
                print(f"\t      \t( total length of cds: {cds_length} ,  {cds_percentage:.2f} % ) ")
                stats_out.write(f"\t( total length of cds: {cds_length} ,  {cds_percentage:.2f} % )\n")
                # contig_repeat_percentage = (contig_gene_bp/contig_lengths[contig][1]) *100
                # print(f"\t{contig_repeat_percentage:.2f}% genes ( {contig_gene_bp} / {contig_lengths[contig][1]} ) bp")
                # print(f"\t{contig_categories}")
    
    else: # show progress bar if not verbose
        print(f" \t -  creating windows for parsed contigs\n")

        for contig in tqdm(contig_lengths.keys()):
            assembly_length_incl += contig_lengths[contig][1]

            windows_contig = get_window_intervals(contig_lengths[contig], window_length)
            all_windows[contig] = windows_contig

            contig_genes = annotation[contig]
            
            contig_abundances, contig_categories, overlap_removed_bp = get_contig_abundances(contig_genes, windows_contig, verbose=False, gene_number = True, gene_density=True, calc_gene_density = calc_gene_density, filter_overlap_previous=False)
            
            # print(contig_abundances)
            if verbose:
                contig_gene_bp = int(sum([sum(list(window_abundances.values())) for window_abundances in contig_abundances.values()]))
                cds_length += contig_gene_bp
            
            all_categories.extend(contig_categories)
            gene_abundances[contig] = contig_abundances
    
    all_categories = list(set(all_categories))

    if verbose:
        print(f"\n\tIn total, there are {len(all_categories)} unique feature categories: {all_categories} (should be only one for the genes)")
        cds_percentage = (cds_length/assembly_length_incl) * 100
        # print(f"\t{assembly_length_incl} bp of sequence were parsed, of which {cds_length} bp were annotated as genes ({cds_percentage:.2f}%)")
    
    return gene_abundances, all_categories, all_windows





def average_over_gene_window_abundances(gene_abundances, merge_gene_windows):
    """
    uses the output from get_assembly_gene_numbers and returns a dictionary of the same dimensions 
    but where every [merge_gene_windows] are grouped and averaged so that they each contain the mean number of genes between all of them
    """

    genes_mean = gene_abundances

    for contig in gene_abundances.keys():
        num_windows = len(gene_abundances[contig].keys())
        windows_in_last_section = num_windows%merge_gene_windows
        
        # print(f"{contig} :  {num_windows} windows")
        for window_ind in range(0, num_windows, merge_gene_windows):
            # print(window_ind)
            window_interval_genes:list[float] = []
            
            windows_in_interval = merge_gene_windows
            if window_ind+windows_in_interval > num_windows:
                windows_in_interval =   num_windows - window_ind

            # get mean from original gene abundances
            for window in range(window_ind, window_ind+windows_in_interval):
                try:
                    window_interval_genes.append(list(gene_abundances[contig][window].values())[0])
                    # print(f"\t{window} :  {list(gene_abundances[contig][window].values())[0]}")
                except:
                    window_interval_genes.append(0)
                    # print(f"\t{window} :  0")
            mean_interval_genes = int(mean(window_interval_genes))
            # print(mean_interval_genes)

            # fill new dictionary with mean abundances
            for window in range(window_ind, window_ind+windows_in_interval):
                genes_mean[contig][window]["Gene"] = mean_interval_genes
                # print(f"\t{window} :  {mean_interval_genes}")
        
    return genes_mean

    

def get_max_genes_in_window(gene_abundances):
    """
    see how many genes are in the window with the most genes to easily scale the Y axis for the gene line when plotting
    """
    max_genes = 0
    for contig in gene_abundances.keys():
        contig_abundances = gene_abundances[contig]
        for window in contig_abundances.keys():
            try:
                if contig_abundances[window]["Gene"] > max_genes:
                    max_genes = contig_abundances[window]["Gene"]
            except:
                pass
    return max_genes
        

def make_dict_to_list_for_pandas(abundances_dict, abundances_categories, windows_dict, gene_abundances = {}, gene_abundances2 = {}, unmasked_bp_dict = {}, verbose = False):
    """ 
    take the dictionary with repeats
    {  
    contig : { 
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

    and potentially include other information, 

    and make it into a list of lists for each window:
    [   
        [ contig, window_ind, window_length, cds_basepairs, cds_basepairs2, unmasked bp, repeat_category1, repeat_category2, ... ],
        [ contig, window_ind, window_length, cds_basepairs, cds_basepairs2, unmasked bp, repeat_category1, repeat_category2, ... ],
        ...
    ]
    """
    abundances_list = []
    if verbose:
        print(" *  parse dictionary for export")
        if len(gene_abundances)>0:
            # print(f"\tinclude one gene annotation :  {gene_abundances}")
            print(f"\tinclude one gene annotation")
        if len(gene_abundances2)>0:
            # print(f"\tinclude second gene annotation :  {gene_abundances2}")
            print(f"\tinclude second gene annotation")
        if len(unmasked_bp_dict)>0:
            # print(f"\tinclude unmasked basepairs :  {unmasked_bp_dict}")
            print(f"\tinclude unmasked basepairs")

    for contig in tqdm(abundances_dict.keys()):
        abundances_contig = abundances_dict[contig]

        for window_ind in abundances_contig.keys():
            abundances_window = abundances_contig[window_ind]
            window_length = windows_dict[contig][window_ind][1] - windows_dict[contig][window_ind][0] +1

            # get gene numbers
            if len(gene_abundances)>0:
                try:
                    cds_bp = int(gene_abundances[contig][window_ind]["cds"])
                except:
                    cds_bp = 0
            else:
                cds_bp = 0 
            if len(gene_abundances2)>0:
                try:
                    cds_bp2 = int(gene_abundances2[contig][window_ind]["cds"])
                except:
                    cds_bp2 = 0
            else:
                cds_bp2 = 0

            # get unmasked bp 
            if len(unmasked_bp_dict)>0:
                try:
                    unmasked_bp = int(unmasked_bp_dict[contig][window_ind])
                except:
                    print(unmasked_bp_dict[contig][window_ind])
                    raise RuntimeError
            else:
                unmasked_bp = 0

            # start to make the line in the output file
            abundances_line = [contig, window_ind, window_length, cds_bp, cds_bp/unmasked_bp, cds_bp2, cds_bp2/unmasked_bp, unmasked_bp]

            for rep_cat in abundances_categories:
                try:
                    abundances_line.append(int(abundances_window[rep_cat])/window_length)
                except:
                    abundances_line.append(0)
            abundances_list.append(abundances_line)
    
    df_header = ["contig", "window_index", "window_length", "cds_bp", "cds_ratio", "cds_bp2", "cds_ratio2", "unmasked_bp"]
    df_header.extend(abundances_categories)

    abundances_df = pd.DataFrame(abundances_list, columns=df_header)

    return(abundances_df)


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
        window_intervals = get_window_intervals(contig_coords=contig_length_list, window_length=window_length)

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


if __name__ == "__main__":

    start_time_complete = time.perf_counter()

    args = parse_args()

    include_annotation = False
    if args.annotation_gff:
        include_annotation = True
        include_second_annotation = False
        if args.annotation_gff2:
            include_second_annotation = True

    out_filepath = args.masker_outfile
    gff_filepath = args.masker_out_gff

    window_length = args.window_length
    statistics = args.statistics

        

    if len(args.species_name)==0:
        species_name = gff.split_at_second_occurrence(gff_filepath.split("/")[-1])
    else:
        species_name = args.species_name
    
    if statistics:
        if args.plot:
            mode="plot_mode"
        elif args.table:
            mode="table_mode"
        statistics_outfile_name = f"by_contig_statistics_{mode}_{species_name}.txt"
        with open(statistics_outfile_name, "w") as stats_out:
            stats_out.write(f"{species_name} -- by-contig statistics from running mape_repeat_windows.py with --statistics\n")
    
    if args.plot:

        if args.verbose:
            print(f"\n\tPLOT mode --> the output will be a png file\n")

        # get repeat abundances in windows
        species_abundances, species_categories = get_assembly_repeat_abundances(out_filepath, window_length, gff_filepath=gff_filepath, verbose=args.verbose, statistics=args.statistics, filter_overlap_previous_=True)
        
        # get gene density
        species_gene_abundances ={}
        species_gene_abundances2 ={}

        if include_annotation and not include_second_annotation:
            species_gene_abundances, speices_gene_categories = get_assembly_gene_numbers(args.annotation_gff, gff_filepath, window_length, verbose = True, calc_gene_density = args.gene_density)
            avg_species_gene_abundances = average_over_gene_window_abundances(species_gene_abundances, args.merge_gene_windows)
            # print(avg_species_gene_abundances)
            plot_repeat_abundance(species_abundances, species_categories, gff_filepath, window_length, gene_numbers = species_gene_abundances, species_name=args.species_name)

        elif include_annotation and include_second_annotation:
            # first annotation
            species_gene_abundances, speices_gene_categories = get_assembly_gene_numbers(args.annotation_gff, gff_filepath, window_length, verbose = True, statistics=args.statistics, calc_gene_density = args.gene_density)
            avg_species_gene_abundances = average_over_gene_window_abundances(species_gene_abundances, args.merge_gene_windows)
            # second annotation
            species_gene_abundances2, speices_gene_categories2 = get_assembly_gene_numbers(args.annotation_gff2, gff_filepath, window_length, verbose = True, statistics=args.statistics, calc_gene_density = args.gene_density)
            avg_species_gene_abundances2 = average_over_gene_window_abundances(species_gene_abundances2, args.merge_gene_windows)
            # print(avg_species_gene_abundances)
            plot_repeat_abundance(species_abundances, species_categories, gff_filepath, window_length, gene_numbers = species_gene_abundances, gene_numbers2 = species_gene_abundances2, species_name=args.species_name)


        else:
            plot_repeat_abundance(species_abundances, species_categories, gff_filepath, window_length, species_name=args.species_name)
            pass



    elif args.table:

        if args.verbose:
            print(f"\n\tTABLE mode --> the output will be a tsv file\n")

        #turn off the overlap filtering for making the table! this is necessary when you want to do reliable statistics on 
        species_abundances, species_categories = get_assembly_repeat_abundances(out_filepath, window_length, gff_filepath=gff_filepath, verbose=args.verbose, statistics=args.statistics, filter_overlap_previous_=False)

        # get gene density
        species_gene_abundances ={}
        species_gene_abundances2 ={}
        unmasked_dict = {}

        if include_annotation and not include_second_annotation:
            species_gene_abundances, speices_gene_categories, all_windows = get_assembly_gene_abundances(args.annotation_gff, gff_filepath, window_length, verbose = True, calc_gene_density = args.gene_density)
                        
        elif include_annotation and include_second_annotation:
            # first annotation
            species_gene_abundances, speices_gene_categories, all_windows = get_assembly_gene_abundances(args.annotation_gff, gff_filepath, window_length, verbose = True, statistics=args.statistics, calc_gene_density = args.gene_density)
            # second annotation
            species_gene_abundances2, speices_gene_categories2, all_windows = get_assembly_gene_abundances(args.annotation_gff2, gff_filepath, window_length, verbose = True, statistics=args.statistics, calc_gene_density = args.gene_density)

        else:
            pass

        if args.assembly_path:
            unmasked_dict = get_masked_bp_in_window(args.assembly_path, window_length, verbose=args.verbose)
            

        output_df = make_dict_to_list_for_pandas(species_abundances, species_categories, gene_abundances = species_gene_abundances, gene_abundances2 = species_gene_abundances2, windows_dict = all_windows, unmasked_bp_dict = unmasked_dict, verbose = args.verbose) ## --> this function does not work properly yet
        # print(output_df)

        ## save to tsv
  
        outfile_name = f"repeat_statistics_{species_name}_table.tsv"
        
        output_df.to_csv(outfile_name, sep="\t", header=True, index=False)
        print(f"\table saved in the current working directory directory as: {outfile_name}")


    if args.verbose:
        end_time_complete = time.perf_counter()
        execution_time = end_time_complete - start_time_complete
        print(f"* complete runtime: {execution_time:.2f} seconds\n")
