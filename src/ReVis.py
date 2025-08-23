#####
# Analysis and plotting of repeatmasker output, see documentation below or repeatmasker_window_analysis.py -h
# Milena R. Trabert, 2025
#####

import parse_gff as gff
import repeats_windows as rep_windows
import genes_windows as gen_windows
from tqdm import tqdm
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.ticker import FuncFormatter
import pandas as pd
import random
import time
import re

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
where two or more repeats cover the same stretch of sequence, which can result in a >100% repeat coverage in some rep_windows. 
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

    parser.add_argument('--plot_overlap_filtered', action="store_true", help="use the overlap-filtered repeats for plotting")

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

        max_genes = gen_windows.get_max_genes_in_window(gene_numbers)
        y_max_genes = max_genes * 2

        ax3 = ax.twinx()
        ax3.set_ylim(bottom=0, top=y_max_genes)
        ax3.set_ylabel('Number of genes per window', color = genes_line_color, fontsize = fs)
        ax3.tick_params(axis ='y', labelcolor = genes_line_color, labelsize = fs) 
    
        if len(gene_numbers2)>0:
            include_genes_line2 = True

            max_genes2 = gen_windows.get_max_genes_in_window(gene_numbers)
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
            stats_out.write(f"{species_name} -- by-contig statistics from running mape_repeat_rep_windows.py with --statistics\n")
    
    if args.plot:

        if args.verbose:
            print(f"\n\tPLOT mode --> the output will be a png file\n")

        # get repeat abundances in windows
        species_abundances, species_categories = rep_windows.get_assembly_repeat_abundances(out_filepath, window_length, gff_filepath=gff_filepath, verbose=args.verbose, statistics=args.statistics, filter_overlap_previous_=args.plot_overlap_filtered, statistics_outfile_name=statistics_outfile_name, statistics_outfile_name=statistics_outfile_name)
        
        # get gene density
        species_gene_abundances ={}
        species_gene_abundances2 ={}

        if include_annotation and not include_second_annotation:
            species_gene_abundances, speices_gene_categories = gen_windows.get_assembly_gene_numbers(args.annotation_gff, gff_filepath, window_length, verbose = True, calc_gene_density = args.gene_density, statistics_outfile_name=statistics_outfile_name)
            avg_species_gene_abundances = gen_windows.average_over_gene_window_abundances(species_gene_abundances, args.merge_gene_windows)
            # print(avg_species_gene_abundances)
            plot_repeat_abundance(species_abundances, species_categories, gff_filepath, window_length, gene_numbers = species_gene_abundances, species_name=args.species_name)

        elif include_annotation and include_second_annotation:
            # first annotation
            species_gene_abundances, speices_gene_categories = gen_windows.get_assembly_gene_numbers(args.annotation_gff, gff_filepath, window_length, verbose = True, statistics=args.statistics, calc_gene_density = args.gene_density, statistics_outfile_name=statistics_outfile_name)
            avg_species_gene_abundances = gen_windows.average_over_gene_window_abundances(species_gene_abundances, args.merge_gene_windows)
            # second annotation
            species_gene_abundances2, speices_gene_categories2 = gen_windows.get_assembly_gene_numbers(args.annotation_gff2, gff_filepath, window_length, verbose = True, statistics=args.statistics, calc_gene_density = args.gene_density, statistics_outfile_name=statistics_outfile_name)
            avg_species_gene_abundances2 = gen_windows.average_over_gene_window_abundances(species_gene_abundances2, args.merge_gene_windows)
            # print(avg_species_gene_abundances)
            plot_repeat_abundance(species_abundances, species_categories, gff_filepath, window_length, gene_numbers = species_gene_abundances, gene_numbers2 = species_gene_abundances2, species_name=args.species_name)


        else:
            plot_repeat_abundance(species_abundances, species_categories, gff_filepath, window_length, species_name=args.species_name)
            pass



    elif args.table:

        if args.verbose:
            print(f"\n\tTABLE mode --> the output will be a tsv file\n")

        #turn off the overlap filtering for making the table! this is necessary when you want to do reliable statistics on 
        species_abundances, species_categories = rep_windows.get_assembly_repeat_abundances(out_filepath, window_length, gff_filepath=gff_filepath, verbose=args.verbose, statistics=args.statistics, filter_overlap_previous_=False)

        # get gene density
        species_gene_abundances ={}
        species_gene_abundances2 ={}
        unmasked_dict = {}

        if include_annotation and not include_second_annotation:
            species_gene_abundances, speices_gene_categories, all_windows = gen_windows.get_assembly_gene_abundances(args.annotation_gff, gff_filepath, window_length, verbose = True, calc_gene_density = args.gene_density, statistics_outfile_name=statistics_outfile_name)
                        
        elif include_annotation and include_second_annotation:
            # first annotation
            species_gene_abundances, speices_gene_categories, all_windows = gen_windows.get_assembly_gene_abundances(args.annotation_gff, gff_filepath, window_length, verbose = True, statistics=args.statistics, calc_gene_density = args.gene_density, statistics_outfile_name=statistics_outfile_name)
            # second annotation
            species_gene_abundances2, speices_gene_categories2, all_windows = gen_windows.get_assembly_gene_abundances(args.annotation_gff2, gff_filepath, window_length, verbose = True, statistics=args.statistics, calc_gene_density = args.gene_density, statistics_outfile_name=statistics_outfile_name)

        else:
            pass

        if args.assembly_path:
            unmasked_dict = gen_windows.get_masked_bp_in_window(args.assembly_path, window_length, verbose=args.verbose)
            

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
