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

import time
import argparse

def parse_args():
    # Create the parser
    program_description = """
----------------- Quick start:
* compute tables
python3 ReVis_transcript_surroundings.py --compute_tables --out_dir ../../example_data --masker_outfile ../../example_data/bruchidius_siliquastri_repeats.fna.out --annotation_gff ../../example_data/bruchidius_siliquastri.gff --orthogroups ../../example_data/N0.tsv --CAFE5_results ../../example_data/CAFE5_Base_family_results.txt --species_name B_siliquastri --bp 500 --GF_size_percentile 90 --verbose
* read tables for plotting
python3 ReVis_transcript_surroundings.py --plot --out_dir ../../example_data --all_before_table ../../example_data/B_siliquastri_cumulative_repeats_before_all_transcripts.txt --sig_before_table ../../example_data/B_siliquastri_cumulative_repeats_before_sig_transcripts_90th_GF_size_percentile.txt --all_after_table ../../example_data/B_siliquastri_cumulative_repeats_after_all_transcripts.txt --sig_after_table ../../example_data/B_siliquastri_cumulative_repeats_after_sig_transcripts_90th_GF_size_percentile.txt --num_transcripts ../../example_data/B_siliquastri_transcript_numbers.txt --species_name B_siliquastri --verbose --plot_white_background

----------------- Summary:

ReVis_transcript_surroundings.py comptues four tables:
* foreground transcripts sequences:
    * before transcript start
    * after transcript start
* background transcripts sequences:
    * before transcript start
    * after transcript start

Files that need to be included to compute the tables from scratch:
* repeats_out (takes both .out and .ori.out and then acts accordingly, like base ReVis. Not out.gff!)
* annotations_gff
* orthogroups_orthofinder
* CAFE_sig_OGs

Other flags to include:
* bp (how many bp up- and downstream of a transcript to compute and plot)
* gene family size percentile (recommended! for the sig. table, only genes that are part of actually expanding gene families (not just all sig. evolving orthogroups) are included)
* species name (matching one in the orthofinder output!!!)

tables that need to be included if they were previously computed and should only be plotted:
* "background" repeat abundance
    * all_before_transcript
    * all_after_transcript
* "foreground" repeat abundance: one of two options
    * sig_before_transcript
    * sig_after_transcript
This can be one of two kinds. they are formatted the same, but can be computed one of two ways
* CAFE5-significant repeat abundances (all orthogroups)
* CAFE5-significant repeats abundances (only genes that are in actually expanding orthogroups) 
"""

    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, description=program_description)
    
    # make tables or only plot
    output = parser.add_mutually_exclusive_group(required=True)
    output.add_argument("--compute_tables_from_OG", action="store_true", help = """it will compute all the tables required for plotting and then also make the plot""")
    output.add_argument("--plot", action="store_true", help = """it will ONLY plot and you have to pass the table filepaths as input (if you go more than 1kb up/downstream, computing the tables will take a long time)""")
    output.add_argument("--compute_tables_from_list", action="store_true", help = """it will compute all the tables from two lists of transcript IDs, a foreground and a background list (or a 'significant' and 'all' list respectively)""")
    
    parser.add_argument('--out_dir', type=str, help="path to output directory. If not given all files will be saved in current working directory")
    
    # arguments for computing the tables
    parser.add_argument('--masker_outfile', type=str, help="""repeatmasker output file ending in .out (.ori.out, is recommended, but both work, the other one is just slower)""")
    parser.add_argument('--annotation_gff', type=str, help="""genome annotation based on the same assembly as the repeatmasker output""")
    
    input_bg = parser.add_mutually_exclusive_group(required=False)
    input_fg = parser.add_mutually_exclusive_group(required=False)
    input_bg.add_argument('--orthogroups', type=str, help="""hierarchical orthogroups file computed by orthofinder (N0.tsv) matching the results from CAFE5""")
    input_fg.add_argument('--CAFE5_results', type=str, help="""family_results.txt file from CAFE5 with orthogroup names matching the orthofinder output""")
    
    input_bg.add_argument('--all_list', type=str, help="""path to csv file containing a list of gene or transcript IDs from annotation_gff to be used in the all table (background table).""")
    input_fg.add_argument('--sig_list', type=str, help="""path to csv file containing a list of gene or transcript IDs from annotation_gff to be used in the sig table (foreground table).""")

    # arguments for given tables, only do the plotting
    parser.add_argument('--all_before_table', type=str, help="""when only plotting: table with background repeat abundance before genes ('all' genes)""")
    parser.add_argument('--all_after_table', type=str, help="""when only plotting: table with background repeat abundance after genes ('all' genes)""")
    parser.add_argument('--sig_before_table', type=str, help="""when only plotting: table with foreground repeat abundance before genes ('significant' genes)""")
    parser.add_argument('--sig_after_table', type=str, help="""when only plotting: table with foreground repeat abundance after genes ('significant' genes)""")
    parser.add_argument('--num_transcripts', type=str, help="""when only plotting: file containing the number of genes/transcripts used to compute the foreground and background tables (to be able to calculate proportion)""")
    

    parser.add_argument('--species_name', type=str, help="""species identifier string, MUST match one species name in the orthofinder output""")
    parser.add_argument('--bp', type=int, help="how many bp up and downstream of the transcript borders should be included. Default is 500 for testing purposes, but to see patterns you should use at least 5kbp")
    parser.add_argument('--GF_size_percentile', type=int, help="""only gene families in the upper nth percentile of gene family size areincluded in the significant gene families. This helps ensure that only gene families that are really expandingin species_name specifically are included, and not the ones that are significant because they are expanding in other species.The default is 90 percent, which includes almost all gene families except the ones with only 1 or 0 members in most cases""")

    parser.add_argument('--plot_white_background', action="store_true", help="the plot does NOT have a transparent background, but white instead")
    parser.add_argument('--plot_no_legend', action="store_true", help="the plot does NOT include a legend with the colors for all the repeat categories")
    parser.add_argument('--verbose', action="store_true", help="print progress in the command line (recommended, on by default)")


    # Parse the arguments
    args = parser.parse_args()

    ## make defaults
    if not args.GF_size_percentile:
        args.GF_size_percentile = 90
    if not args.bp:
        args.bp = 500
    if args.out_dir[-1] != "/":
        args.out_dir = f"{args.out_dir}/"
    
    ## enforce one type of input files for plotting or tables computation
    if args.compute_tables_from_OG:
        for required in ["masker_outfile", "annotation_gff", "orthogroups", "CAFE5_results", "species_name", "bp"]:
            if getattr(args, required) is None:
                parser.error(f"--{required} is required when using --compute_tables_from_OG")

    elif args.compute_tables_from_list:
        for required in ["masker_outfile", "annotation_gff", "all_list", "sig_list", "bp"]:
            if getattr(args, required) is None:
                parser.error(f"--{required} is required when using --compute_tables_from_list")
    elif args.plot:
        for required in ["all_before_table", "all_after_table","sig_before_table", "sig_after_table"]:
            if getattr(args, required) is None:
                parser.error(f"--{required} is required when using --plot")


    return args


def filter_sig_OGs_by_size(orthoDB_orthogroups:dict, species:str, q:int, verbose=False):
    """
    return a list dict orthogroups, where the GF size in the species is above the q'th percentile
    """
    GF_sizes_species = {}
    sizes = []
    for OG_id, transcripts_list in orthoDB_orthogroups.items():
        GF_sizes_species[OG_id] = len(transcripts_list)
        sizes.append(len(transcripts_list))

    ## calculate size threshold
    OGs_filtered = {}
    sizes = np.array(sizes)
    percentile_size = np.percentile(sizes, q = q)
    if verbose:
        print(f" ---> {species}: \n\tmax gene family: {max(sizes)} \n\tmin gene family: {min(sizes)}\n--> {q}th percentile : {percentile_size}")

    for OG_id, size in GF_sizes_species.items():
        if size > percentile_size:
            OGs_filtered[OG_id]= orthoDB_orthogroups[OG_id]
    if verbose:
        print(f"before filtering: {len(orthoDB_orthogroups)} --> after filtering: {len(OGs_filtered)}")

    return OGs_filtered


def read_transcripts_num_file(transcripts_num_filepath:str):
    num_sig_tr = 0
    num_all_r = 0
    with open(transcripts_num_filepath, "r") as num_file:
        for line_ in num_file.readlines():
            line = line_.strip()
            desc,num = line.split(": ")
            if "number of sig. transcripts: " in line:
                num_sig_tr = int(num)
            elif "number of all transcripts: " in line:
                numm_all_tr = int(num)
            else:
                raise RuntimeError(f"wrong formatting in transcript numbers file!\n\t {line_} does not match the expectation")
    return num_sig_tr, numm_all_tr


def num_transcripts_from_OG_dict(OG_dict:dict):
    """ 
    get the number of transcripts in a dictionary of orthogroups like { OG_id : [transcripts,list] }
    """
    num_transcripts = 0
    for OG_id, transcripts_list in OG_dict.items():
        if transcripts_list == ['']:
            continue
        num_transcripts+= len(transcripts_list)
    return num_transcripts

def plot_TE_abundance(before_filepath:str, after_filepath:str, sig_transcripts:int, general_legend_names = False, all_before_filepath:str = "", all_after_filepath:str = "", all_transcripts:int = 0, filename = "cumulative_repeat_presence_around_transcripts.png", legend = True, plot_white_bg = False):
    """
    plot the cumulative repeat presence per base before and after a transcript (before and after infile paths)
    infiles generated from make_cumulative_TE_table and saved to text file
    """
    plot_transparent_bg = True
    if plot_white_bg:
        plot_transparent_bg = False
    
    before_dict = gff.read_dict_from_file(before_filepath)
    before_dict = { key : [int(v)/sig_transcripts*100 for v in value] for key, value in before_dict.items()}
    # show percentages of they are too high
    for key, value in before_dict.items():
        for v in value:
            v_ = int(v)
            perc = v_/sig_transcripts*100
            if perc>150:
                print(f"{key} : \t {v_}/{sig_transcripts} * 100 = {perc:.2f}%")

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
        fig, ax = plt.subplots(1, 1, figsize=(23, 10))
    else:
        fig, ax = plt.subplots(1, 1, figsize=(21, 12))

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

        rep_label = rep_class.replace("_", " ")
        ax.plot(x_before, before_dict[rep_class], label = rep_label, color = colors[rep_class], linewidth=2)
        ax.plot(x_after, after_dict[rep_class], color = colors[rep_class], linewidth=2)
        
        if all_before_dict !={} and all_after_dict !={} and all_transcripts!=0:
            max_before = max(all_before_dict[rep_class])
            max_after = max(all_after_dict[rep_class])
            if max_before>max_percentage:
                max_percentage=max_before
            if max_after>max_percentage:
                max_percentage=max_after

            ax.plot(x_before, all_before_dict[rep_class], color = colors[rep_class], linestyle = (0, (1, 10)))                    
            ax.plot(x_after, all_after_dict[rep_class], color = colors[rep_class], linestyle = (0, (1, 10)))                
    
    max_percentage=int(max_percentage*1.3)
    if max_percentage == 0 or max_percentage>100:
        max_percentage= 100

    plt.vlines(x= 0, ymin=0, ymax=max_percentage, colors="#000000", linestyles="dashed", label="transcript border", linewidth=3)
    plt.xticks(range(-num_bp, num_bp+1, int(num_bp/5)), fontsize = fs)
    plt.yticks(range(0, max_percentage+1, 10), fontsize = fs)

    species = gff.split_at_second_occurrence(before_filepath.split("/")[-1])
    species = species.replace("_", ". ")

    if legend:
        # plot color legend
        ax.set_xlim([-num_bp, num_bp*1.55])
        legend_colors = ax.legend(loc = "center right", fontsize = fs)
        plt.gca().add_artist(legend_colors)
        # plt.title(f"{species} transcript surroundings {num_bp} bp up and downstream \n({num_sig_transcripts} significant transcripts of {all_transcripts} in CAFE analysis)", fontsize = fs*1.25)
        
    # plot dotted/bold legend
    solid = Line2D([0], [0], color='black', linestyle='-', linewidth=2)
    dotted = Line2D([0], [0], color='black', linestyle=':', linewidth=2)
    handles = [solid, dotted]
    labels = []
    if general_legend_names:
        labels.append(f"foreground transcripts ({num_sig_transcripts})")
        labels.append(f"background transcripts ({all_transcripts})")
    else:
        labels.append(f"significant transcripts ({num_sig_transcripts})")
        labels.append(f"all CAFE transcripts ({all_transcripts})")
    plt.legend(handles, labels, loc = "upper left", fontsize = fs)

    plt.title(f"{species} transcript surroundings {num_bp} bp up and downstream", fontsize = fs*1.25)
    plt.xlabel(f"basepairs upstream and downstream from transcript", fontsize = fs)

    plt.ylabel(f"percent of transcripts in which this base is a repeat", fontsize = fs)
    ax.yaxis.set_major_formatter(FuncFormatter(lambda x, pos: '' if x > 99 or x<1 else f'{int(x)}%'))
    
    plt.tight_layout()
    plt.savefig(filename, dpi = 300, transparent = plot_transparent_bg)
    print("Figure saved in the current working directory directory as: "+filename)

def read_transcript_IDs_csv(filepath:str)->list:
    """ 
    read the csv path with transcript ID lists into an actual list
    """
    tr_list = []
    with open(filepath, "r") as file:
        lines = file.readlines()
        assert len(lines) == 1
        line = lines[0].strip()
        tr_list = line.split(",")
    tr_list = list(set(tr_list))## get rid of duplicates
    return tr_list


if __name__ == "__main__":

    start_time_complete = time.perf_counter()
    args = parse_args()

    repeats_out = args.masker_outfile
    orthoDB_annotation = args.annotation_gff
    orthogroups_orthoDB = args.orthogroups
    sig_orthoDB = args.CAFE5_results

    size_percentile_threshold = args.GF_size_percentile
    species = args.species_name
    verbose = args.verbose

    num_bp = args.bp

    ##########################################################
    ######### make TE abundance tables for plotting ##########
    ##########################################################
    
    ########
    ## tables for significant transcripts.
    ########

    transcript_nums_file = f"{args.out_dir}{species}_transcript_numbers.txt"
    num_sig_transcripts = 0
    num_all_transcripts = 0

    if args.compute_tables_from_list:    
        if verbose:
            print(f"  * making foreground tables from input list for {species}: ")
        all_transcripts_list = read_transcript_IDs_csv(args.all_list)
        num_all_transcripts = len(all_transcripts_list)
        sig_list = read_transcript_IDs_csv(args.sig_list)
        num_sig_transcripts = len(sig_list)
        before_transcript, after_transcript = tr_surrounds.make_cumulative_TE_table(sig_list, n=num_bp, species=species, repeats_annot_path=repeats_out, genome_annot_path=orthoDB_annotation, sig_orthogroups=sig_list, verbose = verbose)

    if args.compute_tables_from_OG:
        if verbose:
            print(f"  * making foreground tables from orthofinder/CAFE output for {species}: ")
        sig_orthoDB_list, all_orthogroups_list = OGs.get_sig_orthogroups(sig_orthoDB)
        orthoDB_orthogroups_dict = OGs.parse_orthogroups_dict(orthogroups_orthoDB, sig_orthoDB_list, species=species)
        
        sig_OGs_size_filtered_dict = filter_sig_OGs_by_size(orthoDB_orthogroups=orthoDB_orthogroups_dict, species=species, q=size_percentile_threshold)
        sig_OGs_size_filtered = list(sig_OGs_size_filtered_dict.keys())
        list_sig_transcripts = tr_surrounds.get_sig_transcripts(sig_OGs_size_filtered_dict)
        num_sig_transcripts = len(list_sig_transcripts)

        ### compute table
        before_transcript, after_transcript = tr_surrounds.make_cumulative_TE_table(orthogroups_orthoDB, n=num_bp, species=species, repeats_annot_path=repeats_out, genome_annot_path=orthoDB_annotation, sig_orthogroups=sig_OGs_size_filtered, verbose = verbose)
    
    if not args.plot:
        ## return values are filepaths
        sig_before_transcript = gff.write_dict_to_file(before_transcript, f"{args.out_dir}{species}_cumulative_repeats_before_sig_transcripts_{size_percentile_threshold}th_GF_size_percentile.txt")
        sig_after_transcript = gff.write_dict_to_file(after_transcript, f"{args.out_dir}{species}_cumulative_repeats_after_sig_transcripts_{size_percentile_threshold}th_GF_size_percentile.txt")

    #########
    ## tables for all background transcripts
    #########
    if args.compute_tables_from_list:    
        if verbose:
            print(f"  * making foreground tables from input list for {species}: ")
        all_transcripts_list = read_transcript_IDs_csv(args.all_list)
        num_all_transcripts = len(all_transcripts_list)
        sig_list = read_transcript_IDs_csv(args.sig_list)
        num_sig_transcripts = len(sig_list)
        before_transcript, after_transcript = tr_surrounds.make_cumulative_TE_table(all_transcripts_list, n=num_bp, species=species, repeats_annot_path=repeats_out, genome_annot_path=orthoDB_annotation, verbose = verbose)

    if args.compute_tables_from_OG:
        if verbose:
            print(f"  * making background tables from orthofinder/CAFE output for {species}: ")
        sig_orthoDB_list, all_orthogroups_list = OGs.get_sig_orthogroups(sig_orthoDB)
        # get number of transcripts to compute percentiles for "background"
        orthoDB_orthogroups_dict = OGs.parse_orthogroups_dict(orthogroups_orthoDB, species=species)
        list_all_transcripts = tr_surrounds.get_sig_transcripts(orthoDB_orthogroups_dict)
        num_all_transcripts = len(list_all_transcripts)

        ## compute table
        before_transcript, after_transcript = tr_surrounds.make_cumulative_TE_table(orthogroups_orthoDB, n=num_bp, species=species, repeats_annot_path=repeats_out, genome_annot_path=orthoDB_annotation)
    
    if not args.plot:
        all_before_transcript = gff.write_dict_to_file(before_transcript, f"{args.out_dir}{species}_cumulative_repeats_before_all_transcripts.txt")
        all_after_transcript = gff.write_dict_to_file(after_transcript, f"{args.out_dir}{species}_cumulative_repeats_after_all_transcripts.txt")
        with open(transcript_nums_file, "w") as tr_num_file:
            tr_num_file.write(f"number of sig. transcripts: {num_sig_transcripts}\n")
            tr_num_file.write(f"number of all transcripts: {num_all_transcripts}")


    elif args.plot:
        ## if only plotting is specified, use the given table filepaths to get all tables and values required
        all_before_transcript = args.all_before_table
        all_after_transcript = args.all_after_table
        sig_before_transcript = args.sig_before_table
        sig_after_transcript = args.sig_after_table

        transcript_nums_file = args.num_transcripts
        num_sig_transcripts, num_all_transcripts = read_transcripts_num_file(transcript_nums_file)

    ######################################################
    ############ plot above computed tables ##############
    ######################################################


    if verbose:
        print(f"\n  * plot {species}")
        print(f"all transcripts: {num_all_transcripts}, sig transcripts: {num_sig_transcripts}")
    # plot_TE_abundance(threshold_before_transcript[species], threshold_after_transcript[species], sig_transcripts = num_sig_transcripts, filename=f"{repeats_plots}cumulative_repeat_presence_around_transcripts_sig_only_{species}_{size_percentile_threshold}th_percentile_GF_size.png")
    plot_legend = True
    if args.plot_no_legend:
        plot_legend=False

    
    plot_TE_abundance(before_filepath = sig_before_transcript, after_filepath=sig_after_transcript, sig_transcripts = num_sig_transcripts, general_legend_names=args.compute_tables_from_list, all_before_filepath=all_before_transcript, all_after_filepath=all_after_transcript, all_transcripts=num_all_transcripts, filename=f"{args.out_dir}{species}_cumulative_repeat_presence_around_transcripts.png", legend=plot_legend, plot_white_bg=args.plot_white_background)
    # print(f"{sig_before_transcript}")
    # break



