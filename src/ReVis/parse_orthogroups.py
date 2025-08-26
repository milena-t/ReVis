# check the position of the significant orthogroups in each annotation:
#  - do they look like tandem duplications or are they dispersed?
#  - Are the n basepairs around the gene enriched for a certain kind of TE compared to the rest of the genome?
#    (pick n like 1e5? pretty short probably)
#  - Are they intron-less?


import parse_gff as gff
import numpy as np
from statistics import mean



def get_orthogroup_sizes(orthogroup_dict, q = 0):
    """
    returns a dict of all orthogroups sizes. Do NOT confuse it with gene family sizes! (get_GF_sizes function instead)
    if a percentile q is specified other than 0, it returns only orthogroups whose size is > than the pth percentile of the size distribution
    """
    
    OG_sizes_dict = {}
    sizes = []
    for OG_id, species_dict in orthogroup_dict.items():
        size = 0
        for transcripts_list in species_dict.values():
            size += len(transcripts_list)
        
        OG_sizes_dict[OG_id] = size
        sizes.append(size)
    
    if q == 0:
        return OG_sizes_dict
    else:
        OG_sizes_filtered = {}
        sizes = np.array(sizes)
        percentile_size = np.percentile(sizes, q = q)
        for OG_id, size in OG_sizes_dict.items():
            if size > percentile_size:
                OG_sizes_filtered[OG_id] = size

        return OG_sizes_filtered


def parse_CAFE_output(filepath):
    """
    parse the CAFE family_results.tsv output to a dictionary with OG_id:str keys and p-value:float values
    """
    out_dict = {}
    with open(filepath, "r") as file:
        next(file) # skip first line with file headers
        for line in file.readlines():
            line = line.strip().split("\t")
            out_dict[line[0]] = float(line[1])
    return out_dict


def get_sig_orthogroups(filepath, p_sig = 0.05):
    """ 
    get a list of significant orthogroup IDs from CAFE output
    """
    sig_list = []
    all_list = []
    with open(filepath, "r") as file:
        next(file) # skip first line
        for line in file:
            orthogroup, p_value, sig_bool = line.strip().split("\t")
            if float(p_value)<p_sig:
                sig_list.append(orthogroup)
            all_list.append(orthogroup)

    return(sig_list, all_list)


def parse_orthogroups_dict(filepath, sig_list:list[str] = [], species = ""):
    """
    Get a dictionary of all the orthogroups like below. if a list of significant orthogroups from CAFE is included then only parse those ones.
    
    { 
    orthogroup1_ID : 
        { species : [ transcript1 , transcript2 , ...],
          species : [ transcript1 , transcript2 , ...],
          ...
        }
      ...
    }

    if a species name is specified it will only read the orthogroups of one species, resulting in just { orthogroup_ID : [ transcript1, ...]}
    """
    out_dict = {}
    # og_df = pd.read_csv(filepath, sep="\t")
    # print(og_df)

    ## There's some orthogroups with duplicate IDs, and I don't know why. count how many there are
    old_og_number = 0
    duplicates_count = 0

    with open(filepath, "r") as N0_infile:
        N0_infile = N0_infile.readlines()
        headers = N0_infile[0].strip().split("\t")
        headers =headers[3:] # remove non-species headers for the orthogroup identifiers
        headers_clean = [gff.split_at_second_occurrence(header) for header in headers]

        # if a weird species name is entered
        if len(species)>0 and (species not in headers_clean or species not in headers):
            try:
                species = [head for head in headers_clean if species in head][0]

            except:
                raise RuntimeError(f"""the species name {species} is not represented in the headers,
please pick one of the header names instead: \n\t{headers} \nor:\n\t{headers_clean}""")

        for i, orthogroup_line in enumerate(N0_infile[1:]):
            orthogroup_line = orthogroup_line.strip().split("\t")
            if len(orthogroup_line[3:])< len(headers): # the orthogroups line cuts off after the last species with members
                len_diff = len(headers)-len(orthogroup_line[3:])
                for i in range(len_diff):
                    orthogroup_line.append('')
                
            orthogroup = orthogroup_line[0] # column 0 for the hierarchical one, otherwise 1 for the old one (which is deprecated because of duplicates)
            orthogroup_line=orthogroup_line[3:]

            if len(sig_list)>0 and orthogroup not in sig_list:
                continue # skip this orthogroup
            
            if species == "":
                out_dict[orthogroup] = {}
                for column in range(len(orthogroup_line)):
                    species_col = gff.split_at_second_occurrence(headers[column])
                    try:
                        if len(orthogroup_line[column].split(", "))>0:
                            # the transcripts contain a "species_name_" prefix, remove that here
                            out_dict[orthogroup][species_col] = [transcript.replace(f"{species_col}_", "") for transcript in orthogroup_line[column].split(", ")]
                        else:
                            out_dict[orthogroup][species_col] = []
                    except:
                        out_dict[orthogroup][species_col] = []

            else: 
                try:
                    column = headers.index(species)
                except:
                    column = headers_clean.index(species)

                try:
                    transcripts_list = orthogroup_line[column].split(", ")
                    OG_species = [transcript.replace(f"{species}_", "") for transcript in transcripts_list]
                except:
                    OG_species = ['']
                    
                out_dict[orthogroup] = OG_species
            

    # print(f"{duplicates_count} orthogroups with duplicate names" )
    return(out_dict)


def get_GF_sizes(orthogroups_dict):
    """
    modify the dictionary from the parse_orthogroups function
    
    {
        orthogroup : {
            species1 : num,
            species2 : num,
            ...
        }
        orthogroup2 : {
            species1 : num,
            species2 : num,
            ...
        }
    }
    """
    orthogroups_modified = {}
    for orthogroup_id , species_dict in orthogroups_dict.items():
        species_dict_nums = { species : len(transcripts_list) for species, transcripts_list in species_dict.items() if transcripts_list!=['']}
        orthogroups_modified[orthogroup_id] = species_dict_nums

    return orthogroups_modified

            
def get_mean_GF_size(orthogroups_dict:dict, species:str):
    GF_sizes = []
    for gene_families in orthogroups_dict.values():
        try:
            num_families = len(gene_families[species])
        except:
            num_families = 0
        GF_sizes.append(num_families)
    
    return mean(GF_sizes)



def get_OG_sizes_by_species(orthoDB_orthogroups:dict, orthoDB_sig_dict:dict, native_orthogroups:dict, native_sig_dict:dict, outfile_path = "/Users/miltr339/work/PhD_code/PhD_chapter1/data/mean_orthogroup_sizes.csv"):
    """
    Make table for the R statistical analysis, with mean orthogroup size and mean CAFE sig. orthogroup size 
    """

    assert type(orthoDB_orthogroups) == dict
    assert type(native_orthogroups) == dict
    header = f"species, uniform_all_mean_size, uniform_sig_mean_size, native_all_mean_size, native_sig_mean_size\n"

    # get complete species list:
    OGs_all_list = list(orthoDB_orthogroups.keys())
    species_list = []
    for OG_id in OGs_all_list:
        species_list.extend(list(orthoDB_orthogroups[OG_id].keys()))
    species_list = list(set(species_list))


    with open(outfile_path,"w") as outfile:
        outfile.write(header)
        for species in species_list:
            orthoDB_all_mean = get_mean_GF_size(orthoDB_orthogroups, species)
            orthoDB_sig_mean = get_mean_GF_size(orthoDB_sig_dict, species)
            native_all_mean = get_mean_GF_size(native_orthogroups, species)
            native_sig_mean = get_mean_GF_size(native_sig_dict, species)

            outfile_line = f"{species}, {orthoDB_all_mean}, {orthoDB_sig_mean}, {native_all_mean}, {native_sig_mean}\n"
            outfile.write(outfile_line)
    
    print(f"outfile written to: {outfile_path}")
    
    





if __name__ == "__main__":
    
    repeats_dir = "/Users/milena/work/repeatmasking/repeat_gffs/"
    repeats_out = {
        "A_obtectus" : f"{repeats_dir}A_obtectus_assembly_genomic.fna.out",
        "A_verrucosus" : f"{repeats_dir}A_verrucosus_assembly_genomic.fna.out",
        "B_siliquastri" : f"{repeats_dir}B_siliquastri_assembly_genomic.fna.out",
        "C_analis" : f"{repeats_dir}C_analis_assembly_genomic.fna.out",
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

    orthogroups_native_filepath = "/Users/miltr339/work/PhD_code/PhD_chapter1/data/orthofinder_native/N0.tsv"
    orthogroups_orthoDB_filepath = "/Users/miltr339/work/PhD_code/PhD_chapter1/data/orthofinder_uniform/N0.tsv"
    sig_native = "/Users/miltr339/work/PhD_code/PhD_chapter1/data/CAFE_native_Base_Family_results.txt"
    sig_orthoDB = "/Users/miltr339/work/PhD_code/PhD_chapter1/data/CAFE_uniform_Base_Family_results.txt"

    print("orthoDB : ")
    sig_orthoDB_list, all_orthogroups_list = get_sig_orthogroups(sig_orthoDB)
    orthoDB_orthogroups = parse_orthogroups_dict(orthogroups_orthoDB_filepath)
    print(f"{len(sig_orthoDB_list)} significant orthogroups in CAFE output")
    orthoDB_orthogroups_sig_only = parse_orthogroups_dict(orthogroups_orthoDB_filepath, sig_orthoDB_list)
    # print(f"{len(orthoDB_orthogroups)} orthogroups in total, {len(orthoDB_orthogroups_sig_only)} left after filtering for only CAFE-significant ones")

    print("\nnative :")
    sig_native_list, all_orthogroups_list = get_sig_orthogroups(sig_native)
    native_orthogroups = parse_orthogroups_dict(orthogroups_native_filepath)
    native_orthogroups_sig_only = parse_orthogroups_dict(orthogroups_native_filepath, sig_native_list)

    print("\ncompute_sizes")
    print(type(orthoDB_orthogroups))
    get_OG_sizes_by_species(orthoDB_orthogroups, orthoDB_orthogroups_sig_only, native_orthogroups, native_orthogroups_sig_only)