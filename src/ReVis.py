#####
# Analysis and plotting of repeatmasker output, see documentation below or repeatmasker_window_analysis.py -h
# Milena R. Trabert, 2025
#####

import parse_gff as gff
import repeats_windows as rep_windows
import genes_windows as gen_windows
import plot_stacked_hist as plot_windows

from tqdm import tqdm
import pandas as pd
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

python3 ReVis.py --masker_outfile your_assembly.fna.ori.out --masker_out_gff your_assembly.fna.out.gff --species_name Your_species --window_length 1e6 --plot_overlap_filtered --verbose --plot

------------------------- dependencies

 *  On uppmax, load biopython/1.80-py3.10.8 to use argparse (base python doesn't include it).
 *  Libraries imported in this script and in parse_repeats.py and parse_gff.py:
        - sys, re, os, subprocess, argparse
        - tqdm, random, time
        - dataclasses, enum, pandas, tempfile, SeqIO
        - matplotlib
        (they can all be installed through pip. Sorry if I forgot any, I'm sure the compiler will tell you)

------------------------- documentation

This script analyzes repeatmasker output. You can run in two modes, plot mode or table mode. Plot mode returns 
a plot with a stacked histogram of all the repeat categories (overlap filtered! Repeat annotations can be on top of each other, 
where two or more repeats cover the same stretch of sequence, which can result in a >100% repeat coverage in some rep_windows. 
You can filter out repeat categories that overlap with others so that each base is only covered by one repeat annotation.
This filtering warps the category proportion in some windows, because it is likely that there are more bases of some categories 
removed than others. See get_repeat_abundance function documentation for details). If you include a genome annotation, 
a line for the number of annotated genes in the same windows is added. You can add up to two annotations to compare them.
Table mode returns the same by-window information as a tsv file with the proportion of basepairs in each window 
covered by each repeat category (NOT overlap filtered, so there can be more bp covered by all repeats in a window than the number
of masked bp or even the window length) and also the number of bp and ratio covered by coding regions (exons). If the masked
assembly is given it will also include the number of unmasked bp in each window.

It takes as input two of the repeatmasker output files, *.ori.out (but just *.out also works, only slower), and *.out.gff

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
    parser.add_argument('--out_dir', type=float, help="path to output directory. If not given all files will be saved in current working directory")
    parser.add_argument('--merge_gene_windows', type=int, help="""If an annotation is included to show gene density, 
    the gene density is likely difficult to interpret visually if you show it for each repeat-window.
    Here you can choose how many windows you would like to average over, 
    default = 5
    put 1 if you don't want any averaging""")

    parser.add_argument('--plot_overlap_filtered', action="store_true", help="use the overlap-filtered repeats for plotting")
    parser.add_argument('--plot_white_background', action="store_true", help="the plot does NOT have a transparent background, but white instead")

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
    plot_transparent_backrgound = True
    if args.plot_white_background:
        plot_transparent_backrgound = False

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
    else:
        statistics_outfile_name = ""

    if args.plot:

        if args.verbose:
            print(f"\n\tPLOT mode --> the output will be a png file\n")

        # get repeat abundances in windows
        species_abundances, species_categories = rep_windows.get_assembly_repeat_abundances(out_filepath, window_length, gff_filepath=gff_filepath, verbose=args.verbose, statistics=args.statistics, filter_overlap_previous_=args.plot_overlap_filtered, statistics_outfile_name=statistics_outfile_name)
        
        # get gene density
        species_gene_abundances ={}
        species_gene_abundances2 ={}

        if include_annotation and not include_second_annotation:
            species_gene_abundances, speices_gene_categories = gen_windows.get_assembly_gene_numbers(args.annotation_gff, gff_filepath, window_length, verbose = True, calc_gene_density = args.gene_density, statistics_outfile_name=statistics_outfile_name)
            avg_species_gene_abundances = gen_windows.average_over_gene_window_abundances(species_gene_abundances, args.merge_gene_windows)
            # print(avg_species_gene_abundances)
            plot_windows.plot_repeat_abundance(species_abundances, species_categories, gff_filepath, window_length, transparent_bg = plot_transparent_backrgound, gene_numbers = species_gene_abundances, species_name=args.species_name)

        elif include_annotation and include_second_annotation:
            # first annotation
            species_gene_abundances, speices_gene_categories = gen_windows.get_assembly_gene_numbers(args.annotation_gff, gff_filepath, window_length, verbose = True, statistics=args.statistics, calc_gene_density = args.gene_density, statistics_outfile_name=statistics_outfile_name)
            avg_species_gene_abundances = gen_windows.average_over_gene_window_abundances(species_gene_abundances, args.merge_gene_windows)
            # second annotation
            species_gene_abundances2, speices_gene_categories2 = gen_windows.get_assembly_gene_numbers(args.annotation_gff2, gff_filepath, window_length, verbose = True, statistics=args.statistics, calc_gene_density = args.gene_density, statistics_outfile_name=statistics_outfile_name)
            avg_species_gene_abundances2 = gen_windows.average_over_gene_window_abundances(species_gene_abundances2, args.merge_gene_windows)
            # print(avg_species_gene_abundances)
            plot_windows.plot_repeat_abundance(species_abundances, species_categories, gff_filepath, window_length, transparent_bg = plot_transparent_backrgound, gene_numbers = species_gene_abundances, gene_numbers2 = species_gene_abundances2, species_name=args.species_name)


        else:
            plot_windows.plot_repeat_abundance(species_abundances, species_categories, gff_filepath, window_length, transparent_bg = plot_transparent_backrgound, species_name=args.species_name)
            pass

        print(f"plot white background: {args.plot_white_background}")


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
