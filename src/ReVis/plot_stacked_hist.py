import parse_gff as gff
import genes_windows as gen_windows

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.ticker import FuncFormatter

import random


def random_hex_color():
    """
    custom function to format random hex color code
    """
    return "#{:06x}".format(random.randint(0, 0xFFFFFF))



def plot_repeat_abundance(species_abundances, species_categories, gff_filepath, window_length, transparent_bg, gene_numbers = {}, gene_numbers2 = {}, species_name = "", y_label = "repeat abundance", output_dir  = "", plot_overlap_filtered = True):

    print("\n* start plotting...")

    if len(species_name)==0:
        species_name = gff.split_at_second_occurrence(gff_filepath.split("/")[-1])
        #species_name = species_name.replace("_", ". ")
    if output_dir == "./":
        filename = f"repeat_abundance_in_{species_name}.png"
    else:    
        filename = f"{output_dir}repeat_abundance_in_{species_name}.png"
    if len(gene_numbers)>0:
        filename = f"{output_dir}repeat_abundance_with_gene_numbers_in_{species_name}.png"

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
    
    contig_lengths_all = gen_windows.get_contig_lengths(gff_filepath)
    num_contigs = len(contig_lengths_all)
    assembly_length = sum([lengths[1] for lengths in contig_lengths_all.values()])

    contig_lengths = {contig : length for contig, length in contig_lengths_all.items() if length[1]>window_length}
    contig_lengths_keys = {length[1] : contig for contig, length in contig_lengths_all.items() if length[1]>window_length}
    lengths = [lengths[1] for lengths in contig_lengths.values()]
    lengths.sort(reverse=True)
    contig_names_sorted_by_length = [contig_lengths_keys[length] for length in lengths]

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

    for contig in contig_names_sorted_by_length:
#    for contig in contig_lengths.keys():
        window_intervals = gen_windows.get_window_intervals(contig_lengths[contig], window_length)
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

            window_intervals_genes = gen_windows.get_window_intervals(contig_lengths[contig], window_length)
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
            ax3.yaxis.set_major_formatter(FuncFormatter(lambda x, pos: '' if x == 0 else int(x)))
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
        legend_lines = plt.legend(repeat_lines, ["1", "2"], loc = "upper left", fontsize = fs*0.75, title = "genome annotation", title_fontsize = fs*0.85)
        plt.gca().add_artist(legend_lines)
    elif include_genes_line and not include_genes_line2:
        repeat_lines = ax3.get_lines()
        legend_lines = plt.legend(repeat_lines, ["genome annotation"], loc = "upper left", fontsize = fs*0.75, title_fontsize = fs*0.85)
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
    new_lower_limit = x_limits[0] - 0.17 * (x_limits[1] - x_limits[0])
    if len(gene_numbers)>0:
        new_lower_limit = x_limits[0] - 0.2 * (x_limits[1] - x_limits[0])
    ax.set_xlim(new_lower_limit, x_limits[1])
    
    if plot_overlap_filtered:
        ax.set_ylim(0, 105)

    plt.tight_layout()

    plt.savefig(filename, dpi = 300, transparent = transparent_bg)
    print("\tFigure saved in the current working directory directory as: "+filename)

    # plt.show()
