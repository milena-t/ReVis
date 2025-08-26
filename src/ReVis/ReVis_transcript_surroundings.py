"""
Associate the surroundings of genes that are part of rapidly evolving gene families
for TE presence via a stacked barplot. Each bar is a base, and the coverage is how many TEs of each class are
annotated on that base
"""

import parse_gff as gff
import parse_orthogroups as OGs
import make_transcript_surrounds_table as tr_surrounds

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from matplotlib.lines import Line2D

import argparse

def parse_args():
    # Create the parser
    program_description = """
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

    parser.add_argument('--species_name', type=str, help="""species identifier string, like 'C_maculatus', '_' will be replaced with '. ' 
    (Include this! If not included it will try to parse it automatically from the start of filenames, 
    which will probably only work for how i named my files)""")
    parser.add_argument('--window_length', type=float, required = True, help="window length (scientific notation like 1e6 is ok)")
    parser.add_argument('--out_dir', type=str, help="path to output directory. If not given all files will be saved in current working directory")
    parser.add_argument('--merge_gene_windows', type=int, help="""If an annotation is included to show gene density, 
    the gene density is likely difficult to interpret visually if you show it for each repeat-window.
    Here you can choose how many windows you would like to average over, 
    default = 5
    put 1 if you don't want any averaging""")

    parser.add_argument('--plot_white_background', action="store_true", help="the plot does NOT have a transparent background, but white instead")

    parser.add_argument('--verbose', action="store_true", help="print progress in the command line (recommended, on by default)")


    # Parse the arguments
    args = parser.parse_args()

    return args






def filter_sig_OGs_by_size(orthoDB_orthogroups:dict, species:str, q:int, verbose=False):
    """
    return a list of orthogroup IDs, where the GF size in the species is above the q'th percentile
    """
    GF_sizes_species = {}
    sizes = []
    for OG_id, transcripts_list in orthoDB_orthogroups.items():
        GF_sizes_species[OG_id] = len(transcripts_list)
        sizes.append(len(transcripts_list))

    ## calculate size threshold
    OGs_filtered = []
    sizes = np.array(sizes)
    percentile_size = np.percentile(sizes, q = q)
    if verbose:
        print(f" ---> {species}: \n\tmax gene family: {max(sizes)} \n\tmin gene family: {min(sizes)}\n--> {q}th percentile : {percentile_size}")

    for OG_id, size in GF_sizes_species.items():
        if size > percentile_size:
            OGs_filtered.append(OG_id)
    if verbose:
        print(f"before filtering: {len(orthoDB_orthogroups)} --> after filtering: {len(OGs_filtered)}")

    return OGs_filtered




def plot_TE_abundance(before_filepath:str, after_filepath:str, sig_transcripts:int, all_before_filepath:str = "", all_after_filepath:str = "", all_transcripts:int = 0, filename = "cumulative_repeat_presence_around_transcripts.png", legend = True):
    """
    plot the cumulative repeat presence per base before and after a transcript (before and after infile paths)
    infiles generated from make_cumulative_TE_table and saved to text file
    """
    
    before_dict = gff.read_dict_from_file(before_filepath)
    before_dict = { key : [int(v)/sig_transcripts*100 for v in value] for key, value in before_dict.items()}
    after_dict = gff.read_dict_from_file(after_filepath)
    after_dict = { key : [int(v)/sig_transcripts*100 for v in value] for key, value in after_dict.items()}

    all_before_dict = {}
    all_after_dict = {}
    if all_before_filepath !="" and all_after_filepath !="" :
        
        all_before_dict = gff.read_dict_from_file(all_before_filepath)
        all_before_dict = { key : [int(v)/all_transcripts*100 for v in value] for key, value in all_before_dict.items()}
        all_after_dict = gff.read_dict_from_file(all_after_filepath)
        all_after_dict = { key : [int(v)/all_transcripts*100 for v in value] for key, value in all_after_dict.items()}
        
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
        'SINE?': "#72A8D5" , # ruddy blue
        # '' : "#2A618D" , #lapis lazuli
        # dark red-brown
        'Low_complexity' : "#3A3335" , # Jet 
        'Satellite' : "#564D4F" , #Wenge 
        'Simple_repeat' : "#827376" , #Taupe gray
    }

    fs = 25 # set font size
    if legend:
        fig, ax = plt.subplots(1, 1, figsize=(20, 10))
    else:
        fig, ax = plt.subplots(1, 1, figsize=(20, 12))
    rep_classes = list(before_dict.keys())
    num_bp = len(before_dict[rep_classes[0]])
    x_before = range(-num_bp, 0)
    x_after = range(num_bp)

    max_percentage = 0
    for rep_class in rep_classes:

        max_before = max(before_dict[rep_class])
        max_after = max(after_dict[rep_class])
        if max_before>max_percentage:
            max_percentage=max_before
        if max_after>max_percentage:
            max_percentage=max_after

        ax.plot(x_before, before_dict[rep_class], label = rep_class, color = colors[rep_class])
        ax.plot(x_after, after_dict[rep_class], color = colors[rep_class])
        
        if all_before_dict !={} and all_after_dict !={} and all_transcripts!=0:
            max_before = max(all_before_dict[rep_class])
            max_after = max(all_after_dict[rep_class])
            if max_before>max_percentage:
                max_percentage=max_before
            if max_after>max_percentage:
                max_percentage=max_after

            ax.plot(x_before, all_before_dict[rep_class], color = colors[rep_class], linestyle = (0, (1, 10)))                    
            ax.plot(x_after, all_after_dict[rep_class], color = colors[rep_class], linestyle = (0, (1, 10)))                
    
    max_percentage=int(max_percentage*1.4)
    if max_percentage == 0 or max_percentage>100:
        max_percentage= 100

    plt.vlines(x= 0, ymin=0, ymax=max_percentage, colors="#000000", linestyles="dashed", label="transcript border")
    plt.xticks(range(-num_bp, num_bp+1, int(num_bp/5)), fontsize = fs)
    plt.yticks(range(0, max_percentage+1, 10), fontsize = fs)

    species = gff.split_at_second_occurrence(before_filepath.split("/")[-1])
    species = species.replace("_", ". ")

    if legend:
        ax.set_xlim([-num_bp, num_bp*1.35])
        plt.legend(loc = "upper right", fontsize = fs)
        plt.title(f"{species} transcript surroundings {num_bp} bp up and downstream \n({num_sig_transcripts} significant transcripts of {all_transcripts} in CAFE analysis)", fontsize = fs*1.25)
    else:
        
        solid = Line2D([0], [0], color='black', linestyle='-', linewidth=2)
        dotted = Line2D([0], [0], color='black', linestyle=':', linewidth=2)
        handles = [solid, dotted]
        labels = []
        # handles.append(mpatches.Patch(fill=False, linestyle=None))
        # handles.append(mpatches.Patch(fill=False, linestyle=None))
        labels.append(f"significant transcripts ({num_sig_transcripts})")
        labels.append(f"all CAFE transcripts ({all_transcripts})")
        plt.legend(handles, labels, loc = "upper center", fontsize = fs)
        plt.title(f"{species} transcript surroundings {num_bp} bp up and downstream", fontsize = fs*1.25)
    
    plt.xlabel(f"basepairs upstream and downstream from transcript", fontsize = fs)

    plt.ylabel(f"percent of transcripts in which this base is a repeat", fontsize = fs)
    ax.yaxis.set_major_formatter(FuncFormatter(lambda x, pos: '' if x > 99 and x<1 else f'{int(x)}%'))

    plt.tight_layout()
    plt.savefig(filename, dpi = 300, transparent = True)
    print("Figure saved in the current working directory directory as: "+filename)



def filepaths():

    repeats_dir = "/Users/milena/work/repeatmasking/repeat_gffs/"
    repeats_out = {
        "A_obtectus" : f"{repeats_dir}A_obtectus_assembly_genomic.fna.out",
        "A_verrucosus" : f"{repeats_dir}A_verrucosus_assembly_genomic.fna.out",
        "B_siliquastri" : f"{repeats_dir}B_siliquastri_assembly_genomic.fna.out",
        # "C_analis" : f"{repeats_dir}C_analis_assembly_genomic.fna.out",
        "C_chinensis" : f"{repeats_dir}C_chinensis_assembly_genomic.fna.out",
        "C_maculatus" : f"{repeats_dir}C_maculatus_assembly_genomic.fna.out",
        "C_septempunctata" : f"{repeats_dir}C_septempunctata_assembly_genomic.fna.out",
        "D_melanogaster" : f"{repeats_dir}D_melanogaster_assembly_genomic.fna.out",
        "D_ponderosae" : f"{repeats_dir}D_ponderosae_assembly_genomic.fna.out",
        "I_luminosus" : f"{repeats_dir}I_luminosus_assembly_genomic.fna.out",
        "P_pyralis" : f"{repeats_dir}P_pyralis_assembly_genomic.fna.out",
        "R_ferrugineus" : f"{repeats_dir}R_ferrugineus_assembly_genomic.fna.out",
        "T_castaneum" : f"{repeats_dir}T_castaneum_assembly_genomic.fna.out",
        "T_molitor" : f"{repeats_dir}T_molitor_assembly_genomic.fna.out",
        "Z_morio" : f"{repeats_dir}Z_morio_assembly_genomic.fna.out",
    }

    repeats_dir_work = "/Users/miltr339/work/repeatmasking/repeat_gffs/"
    repeats_out_work = {
        "A_obtectus" : f"{repeats_dir_work}A_obtectus_masking.ori.out",
        "A_verrucosus" : f"{repeats_dir_work}A_verrucosus_masking.ori.out",
        "B_siliquastri" : f"{repeats_dir_work}B_siliquastri_masking.ori.out",
        "C_analis" : f"{repeats_dir_work}C_analis_masking.ori.out",
        "C_chinensis" : f"{repeats_dir_work}C_chinensis_masking.ori.out",
        "C_maculatus" : f"{repeats_dir_work}C_maculatus_superscaffolded_masking.ori.out",
        "C_septempunctata" : f"{repeats_dir_work}C_septempunctata_masking.ori.out",
        "D_melanogaster" : f"{repeats_dir_work}D_melanogaster_masking.ori.out",
        "D_ponderosae" : f"{repeats_dir_work}D_ponderosae_masking.ori.out",
        "I_luminosus" : f"{repeats_dir_work}I_luminosus_masking.ori.out",
        "P_pyralis" : f"{repeats_dir_work}P_pyralis_masking.ori.out",
        "R_ferrugineus" : f"{repeats_dir_work}R_ferrugineus_masking.ori.out",
        "T_castaneum" : f"{repeats_dir_work}T_castaneum_masking.ori.out",
        "T_molitor" : f"{repeats_dir_work}T_molitor_masking.ori.out",
        "Z_morio" : f"{repeats_dir_work}Z_morio_masking.ori.out",
    }

    orthoDB_annot_dir = "/Users/milena/work/orthoDB_annotations/"
    orthoDB_annotations = {
        "A_obtectus" : f"{orthoDB_annot_dir}A_obtectus_orthoDB_filtered.gff",
        "A_verrucosus" : f"{orthoDB_annot_dir}A_verrucosus_orthoDB_filtered.gff",
        "B_siliquastri" : f"{orthoDB_annot_dir}B_siliquastri_orthoDB_filtered.gff",
        "C_analis" : f"{orthoDB_annot_dir}C_analis_orthoDB_filtered.gff",
        "C_chinensis" : f"{orthoDB_annot_dir}C_chinensis_orthoDB_filtered.gff",
        "C_maculatus" : f"{orthoDB_annot_dir}C_maculatus_orthoDB_filtered.gff",
        "C_septempunctata" : f"{orthoDB_annot_dir}C_septempunctata_s_orthoDB_filtered.gff",
        "D_melanogaster" : f"{orthoDB_annot_dir}D_melanogaster_orthoDB_filtered.gff",
        "D_ponderosae" : f"{orthoDB_annot_dir}D_ponderosae_orthoDB_filtered.gff",
        "I_luminosus" : f"{orthoDB_annot_dir}I_luminosus_orthoDB_filtered.gff",
        "P_pyralis" : f"{orthoDB_annot_dir}P_pyralis_orthoDB_filtered.gff",
        "R_ferrugineus" : f"{orthoDB_annot_dir}R_ferrugineus_orthoDB_filtered.gff",
        "T_castaneum_s" : f"{orthoDB_annot_dir}T_castaneum_s_orthoDB_filtered.gff",
        "T_molitor" : f"{orthoDB_annot_dir}T_molitor_orthoDB_filtered.gff",
        "Z_morio" : f"{orthoDB_annot_dir}Z_morio_orthoDB_filtered.gff",
    }

    orthoDB_annot_dir_work = "/Users/miltr339/work/orthoDB_annotations/"
    orthoDB_annotations_work = {
        "A_obtectus" : f"{orthoDB_annot_dir_work}A_obtectus_braker_isoform_filtered.gff",
        "A_verrucosus" : f"{orthoDB_annot_dir_work}A_verrucosus_braker_isoform_filtered.gff",
        "B_siliquastri" : f"{orthoDB_annot_dir_work}B_siliquastri_braker_isoform_filtered.gff",
        "C_analis" : f"{orthoDB_annot_dir_work}C_analis_braker_isoform_filtered.gff",
        "C_chinensis" : f"{orthoDB_annot_dir_work}C_chinensis_braker_isoform_filtered.gff",
        "C_maculatus" : f"{orthoDB_annot_dir_work}C_maculatus_superscaffolded_annotation_isoform_filtered.gff",
        "C_septempunctata" : f"{orthoDB_annot_dir_work}C_septempunctata_braker_isoform_filtered.gff",
        "D_melanogaster" : f"{orthoDB_annot_dir_work}D_melanogaster_braker_isoform_filtered.gff",
        "D_ponderosae" : f"{orthoDB_annot_dir_work}D_ponderosae_braker_isoform_filtered.gff",
        "I_luminosus" : f"{orthoDB_annot_dir_work}I_luminosus_braker_isoform_filtered.gff",
        "P_pyralis" : f"{orthoDB_annot_dir_work}P_pyralis_braker_isoform_filtered.gff",
        "R_ferrugineus" : f"{orthoDB_annot_dir_work}R_ferrugineus_braker_isoform_filtered.gff",
        "T_castaneum" : f"{orthoDB_annot_dir_work}T_castaneum_braker_isoform_filtered.gff",
        "T_molitor" : f"{orthoDB_annot_dir_work}T_molitor_braker_isoform_filtered.gff",
        "Z_morio" : f"{orthoDB_annot_dir_work}Z_morio_braker_isoform_filtered.gff",
    }

    orthogroups_native = "/Users/miltr339/work/orthofinder/native_orthogroups/N0.tsv"
    orthogroups_orthoDB = "/Users/miltr339/work/orthofinder/orthoDB_uniform_masking_orthogroups/N0.tsv"
    # orthogroups_native = "/Users/milena/work/orthofinder/native_orthogroups/N0.tsv"
    # orthogroups_orthoDB = "/Users/milena/work/orthofinder/orthoDB_uniform_masking_orthogroups/N0.tsv"
    sig_native = "/Users/miltr339/Box Sync/code/CAFE/native_from_N0_Base_family_results.txt"
    sig_orthoDB = "/Users/miltr339/Box Sync/code/CAFE/orthoDB_TE_filtered_Base_family_results.txt"

    return repeats_out, repeats_out_work, orthoDB_annotations, orthoDB_annotations_work, orthogroups_native, orthogroups_orthoDB, sig_native, sig_orthoDB



def tables_filepaths():
    work_out_dir = "/Users/miltr339/work/PhD_code/"
    sig_before_transcript = {
        "A_obtectus" : f"{work_out_dir}A_obtectus_cumulative_repeats_before_sig_transcripts.txt",
        "A_verrucosus" : f"{work_out_dir}A_verrucosus_cumulative_repeats_before_sig_transcripts.txt",
        "B_siliquastri" : f"{work_out_dir}B_siliquastri_cumulative_repeats_before_sig_transcripts.txt",
        "C_analis" : f"{work_out_dir}C_analis_cumulative_repeats_before_sig_transcripts.txt",
        "C_chinensis" : f"{work_out_dir}C_chinensis_cumulative_repeats_before_sig_transcripts.txt",
        "C_maculatus" : f"{work_out_dir}C_maculatus_cumulative_repeats_before_sig_transcripts.txt",
        "C_septempunctata" : f"{work_out_dir}C_septempunctata_cumulative_repeats_before_sig_transcripts.txt",
        "D_melanogaster" : f"{work_out_dir}D_melanogaster_cumulative_repeats_before_sig_transcripts.txt",
        "D_ponderosae" : f"{work_out_dir}D_ponderosae_cumulative_repeats_before_sig_transcripts.txt",
        "I_luminosus" : f"{work_out_dir}I_luminosus_cumulative_repeats_before_sig_transcripts.txt",
        "P_pyralis" : f"{work_out_dir}P_pyralis_cumulative_repeats_before_sig_transcripts.txt",
        "R_ferrugineus" : f"{work_out_dir}R_ferrugineus_cumulative_repeats_before_sig_transcripts.txt",
        "T_castaneum" : f"{work_out_dir}T_castaneum_cumulative_repeats_before_sig_transcripts.txt",
        "T_molitor" : f"{work_out_dir}T_molitor_cumulative_repeats_before_sig_transcripts.txt",
        "Z_morio" : f"{work_out_dir}Z_morio_cumulative_repeats_before_sig_transcripts.txt",
    }
    sig_after_transcript = {
        "A_obtectus" : f"{work_out_dir}A_obtectus_cumulative_repeats_after_sig_transcripts.txt",
        "A_verrucosus" : f"{work_out_dir}A_verrucosus_cumulative_repeats_after_sig_transcripts.txt",
        "B_siliquastri" : f"{work_out_dir}B_siliquastri_cumulative_repeats_after_sig_transcripts.txt",
        "C_analis" : f"{work_out_dir}C_analis_cumulative_repeats_after_sig_transcripts.txt",
        "C_chinensis" : f"{work_out_dir}C_chinensis_cumulative_repeats_after_sig_transcripts.txt",
        "C_maculatus" : f"{work_out_dir}C_maculatus_cumulative_repeats_after_sig_transcripts.txt",
        "C_septempunctata" : f"{work_out_dir}C_septempunctata_cumulative_repeats_after_sig_transcripts.txt",
        "D_melanogaster" : f"{work_out_dir}D_melanogaster_cumulative_repeats_after_sig_transcripts.txt",
        "D_ponderosae" : f"{work_out_dir}D_ponderosae_cumulative_repeats_after_sig_transcripts.txt",
        "I_luminosus" : f"{work_out_dir}I_luminosus_cumulative_repeats_after_sig_transcripts.txt",
        "P_pyralis" : f"{work_out_dir}P_pyralis_cumulative_repeats_after_sig_transcripts.txt",
        "R_ferrugineus" : f"{work_out_dir}R_ferrugineus_cumulative_repeats_after_sig_transcripts.txt",
        "T_castaneum" : f"{work_out_dir}T_castaneum_cumulative_repeats_after_sig_transcripts.txt",
        "T_molitor" : f"{work_out_dir}T_molitor_cumulative_repeats_after_sig_transcripts.txt",
        "Z_morio" : f"{work_out_dir}Z_morio_cumulative_repeats_after_sig_transcripts.txt",
    }
    
    all_before_transcript = {
        "A_obtectus" : f"{work_out_dir}A_obtectus_cumulative_repeats_before_all_transcripts.txt",
        "A_verrucosus" : f"{work_out_dir}A_verrucosus_cumulative_repeats_before_all_transcripts.txt",
        "B_siliquastri" : f"{work_out_dir}B_siliquastri_cumulative_repeats_before_all_transcripts.txt",
        "C_analis" : f"{work_out_dir}C_analis_cumulative_repeats_before_all_transcripts.txt",
        "C_chinensis" : f"{work_out_dir}C_chinensis_cumulative_repeats_before_all_transcripts.txt",
        "C_maculatus" : f"{work_out_dir}C_maculatus_cumulative_repeats_before_all_transcripts.txt",
        "C_septempunctata" : f"{work_out_dir}C_septempunctata_cumulative_repeats_before_all_transcripts.txt",
        "D_melanogaster" : f"{work_out_dir}D_melanogaster_cumulative_repeats_before_all_transcripts.txt",
        "D_ponderosae" : f"{work_out_dir}D_ponderosae_cumulative_repeats_before_all_transcripts.txt",
        "I_luminosus" : f"{work_out_dir}I_luminosus_cumulative_repeats_before_all_transcripts.txt",
        "P_pyralis" : f"{work_out_dir}P_pyralis_cumulative_repeats_before_all_transcripts.txt",
        "R_ferrugineus" : f"{work_out_dir}R_ferrugineus_cumulative_repeats_before_all_transcripts.txt",
        "T_castaneum" : f"{work_out_dir}T_castaneum_cumulative_repeats_before_all_transcripts.txt",
        "T_molitor" : f"{work_out_dir}T_molitor_cumulative_repeats_before_all_transcripts.txt",
        "Z_morio" : f"{work_out_dir}Z_morio_cumulative_repeats_before_all_transcripts.txt",
    }
    all_after_transcript = {
        "A_obtectus" : f"{work_out_dir}A_obtectus_cumulative_repeats_after_all_transcripts.txt",
        "A_verrucosus" :f"{work_out_dir}A_verrucosus_cumulative_repeats_after_all_transcripts.txt",
        "B_siliquastri" :f"{work_out_dir}B_siliquastri_cumulative_repeats_after_all_transcripts.txt",
        "C_analis" :f"{work_out_dir}C_analis_cumulative_repeats_after_all_transcripts.txt",
        "C_chinensis" :f"{work_out_dir}C_chinensis_cumulative_repeats_after_all_transcripts.txt",
        "C_maculatus" :f"{work_out_dir}C_maculatus_cumulative_repeats_after_all_transcripts.txt",
        "C_septempunctata" :f"{work_out_dir}C_septempunctata_cumulative_repeats_after_all_transcripts.txt",
        "D_melanogaster" :f"{work_out_dir}D_melanogaster_cumulative_repeats_after_all_transcripts.txt",
        "D_ponderosae" :f"{work_out_dir}D_ponderosae_cumulative_repeats_after_all_transcripts.txt",
        "I_luminosus" :f"{work_out_dir}I_luminosus_cumulative_repeats_after_all_transcripts.txt",
        "P_pyralis" :f"{work_out_dir}P_pyralis_cumulative_repeats_after_all_transcripts.txt",
        "R_ferrugineus" :f"{work_out_dir}R_ferrugineus_cumulative_repeats_after_all_transcripts.txt",
        "T_castaneum" :f"{work_out_dir}T_castaneum_cumulative_repeats_after_all_transcripts.txt",
        "T_molitor" :f"{work_out_dir}T_molitor_cumulative_repeats_after_all_transcripts.txt",
        "Z_morio" :f"{work_out_dir}Z_morio_cumulative_repeats_after_all_transcripts.txt",
    }

    repeats_tables ="/Users/miltr339/work/PhD_code/PhD_chapter1/data/repeats_tables/"
    threshold_before_transcript = {
        'A_obtectus' : f'{repeats_tables}A_obtectus_cumulative_repeats_before_sig_transcripts_90th_GF_size_percentile.txt',
        'A_verrucosus' : f'{repeats_tables}A_verrucosus_cumulative_repeats_before_sig_transcripts_90th_GF_size_percentile.txt',
        'B_siliquastri' : f'{repeats_tables}B_siliquastri_cumulative_repeats_before_sig_transcripts_90th_GF_size_percentile.txt',
        'C_analis' : f'{repeats_tables}C_analis_cumulative_repeats_before_sig_transcripts_90th_GF_size_percentile.txt',
        'C_chinensis' : f'{repeats_tables}C_chinensis_cumulative_repeats_before_sig_transcripts_90th_GF_size_percentile.txt',
        'C_maculatus' : f'{repeats_tables}C_maculatus_cumulative_repeats_before_sig_transcripts_90th_GF_size_percentile.txt',
        'C_septempunctata' : f'{repeats_tables}C_septempunctata_cumulative_repeats_before_sig_transcripts_90th_GF_size_percentile.txt',
        'D_melanogaster' : f'{repeats_tables}D_melanogaster_cumulative_repeats_before_sig_transcripts_90th_GF_size_percentile.txt',
        'D_ponderosae' : f'{repeats_tables}D_ponderosae_cumulative_repeats_before_sig_transcripts_90th_GF_size_percentile.txt',
        'I_luminosus' : f'{repeats_tables}I_luminosus_cumulative_repeats_before_sig_transcripts_90th_GF_size_percentile.txt',
        'P_pyralis' : f'{repeats_tables}P_pyralis_cumulative_repeats_before_sig_transcripts_90th_GF_size_percentile.txt',
        'R_ferrugineus' : f'{repeats_tables}R_ferrugineus_cumulative_repeats_before_sig_transcripts_90th_GF_size_percentile.txt',
        'T_castaneum' : f'{repeats_tables}T_castaneum_cumulative_repeats_before_sig_transcripts_90th_GF_size_percentile.txt',
        'T_molitor' : f'{repeats_tables}T_molitor_cumulative_repeats_before_sig_transcripts_90th_GF_size_percentile.txt',
        'Z_morio' : f'{repeats_tables}Z_morio_cumulative_repeats_before_sig_transcripts_90th_GF_size_percentile.txt',
    }
    threshold_after_transcript = {
        'A_obtectus' : f'{repeats_tables}A_obtectus_cumulative_repeats_after_sig_transcripts_90th_GF_size_percentile.txt',
        'A_verrucosus' : f'{repeats_tables}A_verrucosus_cumulative_repeats_after_sig_transcripts_90th_GF_size_percentile.txt',
        'B_siliquastri' : f'{repeats_tables}B_siliquastri_cumulative_repeats_after_sig_transcripts_90th_GF_size_percentile.txt',
        'C_analis' : f'{repeats_tables}C_analis_cumulative_repeats_after_sig_transcripts_90th_GF_size_percentile.txt',
        'C_chinensis' : f'{repeats_tables}C_chinensis_cumulative_repeats_after_sig_transcripts_90th_GF_size_percentile.txt',
        'C_maculatus' : f'{repeats_tables}C_maculatus_cumulative_repeats_after_sig_transcripts_90th_GF_size_percentile.txt',
        'C_septempunctata' : f'{repeats_tables}C_septempunctata_cumulative_repeats_after_sig_transcripts_90th_GF_size_percentile.txt',
        'D_melanogaster' : f'{repeats_tables}D_melanogaster_cumulative_repeats_after_sig_transcripts_90th_GF_size_percentile.txt',
        'D_ponderosae' : f'{repeats_tables}D_ponderosae_cumulative_repeats_after_sig_transcripts_90th_GF_size_percentile.txt',
        'I_luminosus' : f'{repeats_tables}I_luminosus_cumulative_repeats_after_sig_transcripts_90th_GF_size_percentile.txt',
        'P_pyralis' : f'{repeats_tables}P_pyralis_cumulative_repeats_after_sig_transcripts_90th_GF_size_percentile.txt',
        'R_ferrugineus' : f'{repeats_tables}R_ferrugineus_cumulative_repeats_after_sig_transcripts_90th_GF_size_percentile.txt',
        'T_castaneum' : f'{repeats_tables}T_castaneum_cumulative_repeats_after_sig_transcripts_90th_GF_size_percentile.txt',
        'T_molitor' : f'{repeats_tables}T_molitor_cumulative_repeats_after_sig_transcripts_90th_GF_size_percentile.txt',
        'Z_morio' : f'{repeats_tables}Z_morio_cumulative_repeats_after_sig_transcripts_90th_GF_size_percentile.txt',
    }

    return sig_before_transcript, sig_after_transcript, all_before_transcript, all_after_transcript, threshold_before_transcript, threshold_after_transcript


if __name__ == "__main__":

    # read infile paths 
    repeats_out, repeats_out_work, orthoDB_annotations, orthoDB_annotations_work, orthogroups_native, orthogroups_orthoDB, sig_native, sig_orthoDB = filepaths()

    ##########################################################
    ######### make TE abundance tables for plotting ##########
    ##########################################################
    
    ## Tables are written to csv files since computing them takes a bit
    repeats_tables = "/Users/miltr339/work/PhD_code/PhD_chapter1/data/repeats_tables/"
    ## uv run python3 for quicker runtimes

    

    ## test the function about the GF size filtering
    if False:
        all_species = list(repeats_out.keys())
        for species in all_species:
            sig_orthoDB_list, all_orthogroups_list = OGs.get_sig_orthogroups(sig_orthoDB)
            orthoDB_orthogroups = OGs.parse_orthogroups_dict(orthogroups_orthoDB, sig_orthoDB_list, species=species)
            OG_id_list = filter_sig_OGs_by_size(orthoDB_orthogroups, species, q=90)

            print(f"{species}: \n\tunfiltered list:{sig_orthoDB_list[0:6]} ({len(sig_orthoDB_list)})\n\tfiltered list: {OG_id_list[0:6]} ({len(OG_id_list)})")

    
    ########
    ## tables for significant transcripts
    ########
    size_percentile_threshold = 90
    
    if False:
        all_species = list(repeats_out.keys())
        # failed = ['A_verrucosus', 'C_chinensis', 'D_ponderosae', 'I_luminosus', 'R_ferrugineus', 'T_molitor', 'Z_morio']
        failed = []
        all_species = ['C_maculatus']
        for species in all_species:
            print(f"orthoDB {species}: ")
            
            sig_orthoDB_list, all_orthogroups_list = OGs.get_sig_orthogroups(sig_orthoDB)
            orthoDB_orthogroups = OGs.parse_orthogroups_dict(orthogroups_orthoDB, sig_orthoDB_list, species=species)
            sig_OGs_size_filtered = filter_sig_OGs_by_size(orthoDB_orthogroups=orthoDB_orthogroups, species=species, q=size_percentile_threshold)
            if True:
            #try:
                before_transcript, after_transcript = make_cumulative_TE_table(orthogroups_orthoDB, n=10000, species=species, repeats_annot_path=repeats_out_work[species], genome_annot_path=orthoDB_annotations_work[species], sig_orthogroups=sig_OGs_size_filtered)
                gff.write_dict_to_file(before_transcript, f"{repeats_tables}{species}_cumulative_repeats_before_sig_transcripts_{size_percentile_threshold}th_GF_size_percentile.txt")
                gff.write_dict_to_file(after_transcript, f"{repeats_tables}{species}_cumulative_repeats_after_sig_transcripts_{size_percentile_threshold}th_GF_size_percentile.txt")
            # except: 
            #     failed.append(species)
        print(f"failed species: {failed}")


    #########
    ## tables for all CAFE transcripts
    #########
    if False:
        all_species = list(repeats_out.keys())
        # failed = ['A_verrucosus', 'C_chinensis', 'D_ponderosae', 'I_luminosus', 'R_ferrugineus', 'T_molitor', 'Z_morio']
        failed = []
        for species in all_species:
            print(f"orthoDB {species}: ")
            ## uv run python3 for quicker runtimes
            sig_orthoDB_list, all_orthogroups_list = OGs.get_sig_orthogroups(sig_orthoDB)
            try:
                # before_transcript, after_transcript = make_cumulative_TE_table(orthogroups_orthoDB, n=50, species=species, repeats_annot_path=repeats_out[species], genome_annot_path=orthoDB_annotations[species], sig_orthogroups=sig_orthoDB_list)
                before_transcript, after_transcript = make_cumulative_TE_table(orthogroups_orthoDB, n=10000, species=species, repeats_annot_path=repeats_out_work[species], genome_annot_path=orthoDB_annotations_work[species])
                gff.write_dict_to_file(before_transcript, f"{species}_cumulative_repeats_before_all_transcripts.txt")
                gff.write_dict_to_file(after_transcript, f"{species}_cumulative_repeats_after_all_transcripts.txt")
            except: 
                failed.append(species)
        print(f"failed species: {failed}")



    ######################################################
    ############ plot above computed tables ##############
    ######################################################

    sig_before_transcript, sig_after_transcript, all_before_transcript, all_after_transcript, threshold_before_transcript, threshold_after_transcript = tables_filepaths()

    if True:
        repeats_plots = "/Users/miltr339/work/PhD_code/PhD_chapter1/data/repeats_cumulative_around_transcripts/"
        all_species = list(repeats_out.keys())
        failed_species = []
        for species in all_species:
            # species = "B_siliquastri"
            print(f"\n\n ------------------------------------------------------------- \nplot {species}")
            # get total number of transcripts that are part of significantly rapidly evolving orthogroups in this species
            sig_orthoDB_list, all_orthogroups_list = OGs.get_sig_orthogroups(sig_orthoDB)

            orthoDB_orthogroups = OGs.parse_orthogroups_dict(orthogroups_orthoDB, sig_orthoDB_list, species=species)
            sig_OGs_size_filtered = filter_sig_OGs_by_size(orthoDB_orthogroups=orthoDB_orthogroups, species=species, q=size_percentile_threshold)
            
            all_transcript_IDs = tr_surrounds.get_sig_transcripts(orthoDB_orthogroups)
            num_sig_transcripts_ = len(all_transcript_IDs)
            # print(f"\t{num_sig_transcripts_} significant transcripts according to reading the CAFE and orthoDB output")
            
            all_transcript_IDs = tr_surrounds.make_cumulative_TE_table(orthogroups_orthoDB, n=10000, species=species, repeats_annot_path=repeats_out_work[species], genome_annot_path=orthoDB_annotations_work[species], sig_orthogroups=sig_OGs_size_filtered, count_transcripts=True)
            num_sig_transcripts = len(all_transcript_IDs)
            if num_sig_transcripts_ != num_sig_transcripts:
                print(f"\t!!! {num_sig_transcripts} significant transcrips according to making the table\n\t!!! {num_sig_transcripts_} transcripts according to CAFE and orthoDB output\n")
                failed_species.append(species)

            # plot_TE_abundance(threshold_before_transcript[species], threshold_after_transcript[species], sig_transcripts = num_sig_transcripts, filename=f"{repeats_plots}cumulative_repeat_presence_around_transcripts_sig_only_{species}_{size_percentile_threshold}th_percentile_GF_size.png")

            orthoDB_orthogroups = OGs.parse_orthogroups_dict(orthogroups_orthoDB, all_orthogroups_list, species=species)
            all_transcript_IDs = tr_surrounds.get_sig_transcripts(orthoDB_orthogroups)
            num_all_transcripts_ = len(all_transcript_IDs)
            # print(f"\t{num_all_transcripts_} CAFE transcripts according to reading the CAFE and orthoDB output")
            # write table to output file
            all_transcript_IDs = tr_surrounds.make_cumulative_TE_table(orthogroups_orthoDB, n=10000, species=species, repeats_annot_path=repeats_out_work[species], genome_annot_path=orthoDB_annotations_work[species], count_transcripts=True)
            num_all_transcripts = len(all_transcript_IDs)
            if num_all_transcripts_ != num_all_transcripts:
                print(f"\t!!! {num_all_transcripts} CAFE transcrips according to making the table\n\t!!! {num_all_transcripts_} transcripts according to CAFE and orthoDB output")
                failed_species.append(species)

            plot_TE_abundance(threshold_before_transcript[species], threshold_after_transcript[species], sig_transcripts = num_sig_transcripts, all_before_filepath=all_before_transcript[species], all_after_filepath=all_after_transcript[species], all_transcripts=num_all_transcripts, filename=f"{repeats_plots}cumulative_repeat_presence_around_transcripts_sig_and_all_{species}_{size_percentile_threshold}th_percentile_GF_size.png", legend=False)
            # break
        print(f"\n-----------------------------------------------------\nspecies where the transcript numbers don't match: {failed_species}")



