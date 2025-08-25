from enum import Enum
import time
from tqdm import tqdm
import re



### Utility functions

def split_at_second_occurrence(s, char = "_"): 
    """
    Split a string at the second occurence of "_" (by default)
    This is used when filenames are like 'Species_name_annotation' to extract the species name
    """
    if s.count(char)<2:
        return s
    else:
        first_occurrence = s.find(char)
        second_occurrence = s.find(char, first_occurrence + 1)
        species = s[:second_occurrence]
        return species
    
def write_dict_to_file(dict, filename, header = "", separator = ""):
    with open(filename, "w") as outfile:
        if header != "":
            outfile.write(f"{header}\n")
        if separator == "":
            separator_ = ","
        else:
            separator_ = separator

        for key, value_list in dict.items():
            if type(value_list) == list:
                try:
                    value_str = separator_.join(value_list)
                except:
                    try:
                        value_str = separator_.join([str(element) for element in value_list])
                    except:
                        raise RuntimeError(f"{key} : {value_list} entry could not be parsed")
            else:
                try:
                    value_str = str(value_list)
                except:
                    raise RuntimeError(f"{key} : {value_list} entry could not be parsed")
            if separator == "":
                outfile.write(f"{key}:{value_str}\n")
            else:
                outfile.write(f"{key}{separator}{value_str}\n")
        print("file saved in current working directory as: "+filename)
    
def read_dict_from_file(filepath):
    out_dict = {}
    with open(filepath, "r") as infile:
        lines =infile.readlines()
        for line in lines:
            key, value = line.strip().split(":")
            out_dict[key] = value.split(",")
    return out_dict


def make_species_order_from_tree(newick_tree_path):
    """
    Takes a newick tree and extracts a list of all leaf names in order of occurence
    Makes it possible to align x-axes in later plots with trees
    """
    # Regular expression to extract leaf names
    # This matches strings between commas, parentheses, and before colons.
    leaf_pattern = r'(?<=\(|,)([a-zA-Z0-9_]+)(?=:)'
    with open(newick_tree_path, "r") as newick_tree_file:
        newick_tree_string = newick_tree_file.readlines()[0]
        # print(newick_tree_string)
        leaf_names = re.findall(leaf_pattern, newick_tree_string)
        # leaf names are like "A_obtectus_filtered_proteinfasta" but we only care about the species names in the beginning
        species_names = [split_at_second_occurrence(leaf, "_") for leaf in leaf_names]
    return species_names


def get_single_exon_transcripts(parsed_gff_dict, verbose=False):
    transcripts_single_exon = {}
    transcripts_multi_exon = {}
    for feature_id, feature_object in parsed_gff_dict.items():
        if parsed_gff_dict[feature_id].category == FeatureCategory.Transcript:
            exons = [child_feature for child_feature in feature_object.child_ids_list if parsed_gff_dict[child_feature].category == FeatureCategory.Exon]
            if len(exons) == 1:
                transcripts_single_exon[feature_id] = feature_object
            else:
                transcripts_multi_exon[feature_id] = feature_object
        else:
            continue
    if verbose:
        print(f"{len(transcripts_single_exon)+len(transcripts_multi_exon)} total transcripts\n\t{len(transcripts_multi_exon)} of which are multi-exon, \n\t{len(transcripts_single_exon)} of which are single exon")
    return(transcripts_single_exon, transcripts_multi_exon)


def get_transcript_lengths(subset_transcripts_dict, all_parsed_dict,verbose=False):
    """
    takes the subset of the annotation by get_single_exon_transcripts which is only exons, 
    and the complete parsed annotation
    """
    transcript_lengths = {} 
    for trans_id, attributes in subset_transcripts_dict.items():
        if attributes.category == FeatureCategory.Transcript:
            exons = attributes.child_ids_list
            length = 0
            for exon in exons:
                exon = all_parsed_dict[exon]
                exon_length = int(exon.end) - int(exon.start)
                exon_length = abs(exon_length)
                length += exon_length
            transcript_lengths[trans_id] = length
    if verbose:
        mean_transcript_length = sum(list(transcript_lengths.values()))/len(transcript_lengths)
        print(f"mean transcript length: {mean_transcript_length}")
    return(transcript_lengths)

def print_single_exon_stats(filepath, include_list = True):
    try:
        gff_dict = parse_gff3_general(filepath)
    except Exception as e:
        raise RuntimeError(f"parsing gone wrong! \nError: \n{e}")
    no_transcripts = {key : value for key, value in gff_dict.items() if value.category == FeatureCategory.Transcript}
    no_genes = {key : value for key, value in gff_dict.items() if value.category == FeatureCategory.Gene}
    print(f"total number of transcripts: {len(no_transcripts)}")
    print(f"total number of top-level features (genes): {len(no_genes)}")

    aobt_single_exon, aobt_multi_exon = get_single_exon_transcripts(gff_dict)
    print(f"no. single exon transcripts: {len(aobt_single_exon)}")
    if include_list:
        single_exon_IDs = ",".join(aobt_single_exon)
        print(f"list of single-exon transcript IDs: {single_exon_IDs}")
    print(f"no. multi exon transcripts: {len(aobt_multi_exon)}")
    if include_list:
        multi_exon_IDs = ",".join(aobt_multi_exon)
        print(f"list of single-exon transcript IDs: {multi_exon_IDs }")

    single_exon_transcript_lengths = get_transcript_lengths(subset_transcripts_dict= aobt_single_exon, all_parsed_dict=gff_dict)
    mean_len = sum(single_exon_transcript_lengths.values())/len(single_exon_transcript_lengths)
    print(f"{len(single_exon_transcript_lengths)} single-exon transcripts with average length: {mean_len:.6}")

    multi_exon_transcript_lengths = get_transcript_lengths(subset_transcripts_dict= aobt_multi_exon, all_parsed_dict=gff_dict)
    mean_len = sum(multi_exon_transcript_lengths.values())/len(multi_exon_transcript_lengths)
    print(f"{len(multi_exon_transcript_lengths)} multi-exon transcripts with average length: {mean_len:.6}")





### Genome Feature Class

class FeatureCategory(str,Enum):
    """
    All possible feature categories that can be present in an annotation 
    and example strings of how they could be represented in the gff file
    """
    Gene = "gene"
    CDS = "CDS"
    Exon = "exon"
    Transcript = "transcript"
    Intron = "intron"
    Start_codon = "start_codon"
    Stop_codon = "stop_codon"
    TP_UTR = "three_prime_UTR"
    FP_UTR = "five_prime_UTR"
    RNA = "mRNA" # this and below are specific to the native annotations only
    Region = "region"
    Repeat = "dispersed_repeat"


string_to_category = { # define aliases when some strings should map to the same feature category
    "gene": FeatureCategory.Gene,
    "pseudogene": FeatureCategory.Gene,
    "sequence_feature": FeatureCategory.Region,
    "mobile_genetic_element": FeatureCategory.Region,
    "CDS": FeatureCategory.CDS,
    "cds": FeatureCategory.CDS, 
    "exon": FeatureCategory.Exon,
    "transcript": FeatureCategory.Transcript,
    "primary_transcript": FeatureCategory.Transcript,
    "mRNA": FeatureCategory.Transcript, # I will also use mRNA as "transcript" because even though braker calls transcripts "transcript", 
                                        # other annotation pipelines often say "mRNA" instead when they mean the same thing. 
                                        # Correctly identifying the transcirpt-section of a feature is important for working with e.g. orthofinder output
                                        # which works with transcript IDs only, so if I want to also get transcript stats for native and breaker annotations, 
                                        # I need to group it with the "transcript" feature category instead of the "RNA" feature category.
    "intron": FeatureCategory.Intron,
    "start_codon": FeatureCategory.Start_codon,
    "stop_codon": FeatureCategory.Stop_codon,
    "three_prime_UTR": FeatureCategory.TP_UTR,
    "five_prime_UTR": FeatureCategory.FP_UTR,
    "RNA": FeatureCategory.RNA,
    "antisense_RNA": FeatureCategory.RNA,
    "miRNA": FeatureCategory.RNA,
    "piRNA": FeatureCategory.RNA,
    "rRNA": FeatureCategory.RNA,
    "snoRNA": FeatureCategory.RNA,
    "snRNA": FeatureCategory.RNA,
    "ncRNA": FeatureCategory.RNA,
    "tRNA": FeatureCategory.RNA,
    "lnc_RNA": FeatureCategory.RNA,
    "RNase_MRP_RNA": FeatureCategory.RNA,
    "RNase_P_RNA": FeatureCategory.RNA,
    "SRP_RNA": FeatureCategory.RNA,
    "dispersed_repeat": FeatureCategory.Repeat,
    "region": FeatureCategory.Region, # since it has no parent_ID
    "repeat_region": FeatureCategory.Region, # since it has no parent_ID
    
}

def categorize_string(s):
    """
    sort the category string into one of the existing feature categories
    """
    return string_to_category.get(s)


class Feature:
    """
    any feature that could be present in the gff file
    the gff columns are:
    contig,source,category_,start,stop,score,strandedness,frame,otherproperties
    """
    def __init__(self, feature_id:str, contig:str, category:FeatureCategory, start:int, end:int, strandedness:str, frame:str, parent_id=None):
        
        self.feature_id = feature_id

        self.contig = contig
        self.category = category 
        self.start = start
        self.end = end
        self.strandedness = strandedness
        self.frame = frame

        # since this is general for all features, parent ID is by default None and child ID an empty list, 
        # since they're not always present
        self.parent_id = parent_id 
        self.child_ids_list = []

    def add_child(self, child_id:str):
        self.child_ids_list.append(child_id)

    def __repr__(self):
        return "Feature"
    
    def __str__(self):
        return(
            f"""
            Feature ID: {self.feature_id}, Feature category: {self.category}
            \tParent ID: {self.parent_id}, 
            \tchild ID: {self.child_ids_list}
            \tOn contig: {self.contig}, from {self.start} to {self.end} (on strand {self.strandedness})
            """
        )


def parse_gff3_general(filepath:str, verbose = True, only_genes = False, keep_feature_category=None):
    """
    Read a gff file specified in the filepath and parse it into a dictionary of Feature IDs and instances of the Feature class
    {
        Gene_ID :        Feature (with reference immediate child feature (transcript)) (no parent feature)
        Transcript_ID :  Feature (with parent feature id (gene) and list of child features (exons, mRNA, ...))
        Exon_ID :        Feature (with parent feature (transcript)) (no child features)
    }
    Output can be limited to only Gene features, or include every annotated feature
    """

    genome_annotation:dict = {}

    if verbose:
        start_time = time.perf_counter()
        print(f"File that is being parsed: {filepath}")
        if only_genes:
            print(f"------------------------------------------------------\n!! only_genes = True !!  ==> all non-gene features are excluded\n------------------------------------------------------")
        if keep_feature_category is not None:
            print(f"------------------------------------------------------\n!! only features of {keep_feature_category} included \n------------------------------------------------------")
    with open(filepath, "r") as file:
        linelist = file.readlines()

        # in the attributes column, the "ID" tag and the actual ID can be two things:
        #  * ID=gene1234 (gff2 i think?)
        #  * ID gene1234 (gff3)
        # determine separator (" " or "=") so that I don't have to test in every line
        # look 10 lines from the end to avoid leading or tailing comment lines in the file
        tail_line = linelist[-10].split("\t")[-1].split(";")[0]
        separator = " "
        if "=" in tail_line:
            separator = "="

        count_mRNA = 0
        count_trans = 0
        count_gene = 0
        count_exon = 0
        for line in tqdm(linelist):
            
            # Skip empty lines or lines starting with a comment character
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            # more info on file format and columns here: https://www.ensembl.org/info/website/upload/gff.html?redirect=no
            contig,source,category_,start,stop,score,strandedness,frame,attributes_=[c for c in line.split("\t") if len(c)>0]
            category = categorize_string(category_)

            if only_genes and category != FeatureCategory.Gene:
                # If there's only genes supposed to be included, skip everything that isn't a gene
                continue
            if keep_feature_category != category_ and keep_feature_category != None:
                continue
            # for the child of Gene features, some annotations use "transcript" and some use "mRNA"
            if verbose:
                if category_ == "mRNA":
                    count_mRNA+=1
                elif category_ == "transcript":
                    count_trans +=1
                elif category_ == "gene":
                    count_gene +=1
                elif category_ == "exon":
                    count_exon +=1
            
            attributes={}
            for attr in attributes_.strip().split(";"):
                attr = attr.strip()
                try:
                    key,value=attr.split(separator)[-2:]
                except:
                    continue

                if "," in attr:
                    attributes[key]=value.split(",")[0]
                else:
                    attributes[key]=value
            
            ## check that ID and Parent are detected correctly
            if "ID" not in attributes:
                    raise RuntimeError(f"no id property found for gene in line: {line}")
            if "Parent" not in attributes and not category==FeatureCategory.Gene and not category==FeatureCategory.Region:
                raise RuntimeError(f"feature is not a gene and no parentid property found for feature in line: {line}")
            
            ## add the feature to its parent if relevant
            parent_id = None
            if "Parent" in attributes and keep_feature_category==None:
                parent_id = attributes["Parent"]
                if parent_id in genome_annotation:
                    genome_annotation[parent_id].add_child(attributes["ID"])
                    test_break = True
                else:
                    raise RuntimeError(f"Feature {attributes['ID']} has a parent ID {parent_id} that does not exist prior to it.")           

            ## add Feature to the output dict
            new_feature=Feature(feature_id=attributes["ID"],contig = contig,category=category,start=int(start),end=int(stop),strandedness=strandedness, frame=frame, parent_id=parent_id)
            genome_annotation[new_feature.feature_id]=new_feature

    if verbose:
        end_time = time.perf_counter()
        execution_time = end_time - start_time
        print(f"\tparsing time: {execution_time:.2f} seconds")
        print(f"\t  * gene features: {count_gene}")
        if count_mRNA>0:
            print(f"\t  * mRNA features: {count_mRNA} (0 'transcript' features)")
        if count_trans>0:
            print(f"\t  * transcript features: {count_trans} (0 mRNA features)")
        print(f"\t  * exon features: {count_exon}")
        print(f"\t  * total number of features: {len(genome_annotation)}")

    assert len(genome_annotation)>0      
    return genome_annotation


def parse_gff3_by_contig(filepath:str, verbose = True, featurecategory = FeatureCategory.Gene):
    """
    parse the gff into a dictionary by contig { contig : [ FeatureCategory.Gene ] } 
    This includes by default only genes, not any of the child features, 
    the featurecategory variable is optional and can be set to pick other features than genes, for example exons
    """

    genome_annotation:dict = {}

    if verbose:
        start_time = time.perf_counter()
        print(f" --> parse by contig")
        print(f"File that is being parsed: {filepath}")
        print(f"feature category that was set for the output: {featurecategory}")

    with open(filepath, "r") as file:
        linelist = file.readlines()

        # in the attributes column, the "ID" tag and the actual ID can be two things:
        #  * ID=gene1234 (gff2 i think?)
        #  * ID gene1234 (gff3)
        # determine separator (" " or "=") so that I don't have to test in every line
        # look 10 lines from the end to avoid leading or tailing comment lines in the file
        tail_line = linelist[-10].split("\t")[-1].split(";")[0]
        separator = " "
        if "=" in tail_line:
            separator = "="

        count_feature = 0
        for line in tqdm(linelist):
            
            # Skip empty lines or lines starting with a comment character
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            # more info on file format and columns here: https://www.ensembl.org/info/website/upload/gff.html?redirect=no
            contig,source,category_,start,stop,score,strandedness,frame,attributes_=[c for c in line.split("\t") if len(c)>0]
            category = categorize_string(category_)

            if category != featurecategory:
                # There's only genes supposed to be included, skip everything that isn't a gene
                continue
            
            else:
                if verbose:
                    count_feature += 1

                attributes={}
                for attr in attributes_.strip().split(";"):
                    attr = attr.strip()
                    try:
                        key,value=attr.split(separator)[-2:]
                    except:
                        continue

                    if "," in attr:
                        attributes[key]=value.split(",")[0]
                    else:
                        attributes[key]=value
            
                ## check that ID and Parent are detected correctly
                if "ID" not in attributes:
                        raise RuntimeError(f"no id property found for gene in line: {line}")

                ## add Feature to the output dict
                new_feature=Feature(feature_id=attributes["ID"],contig = contig,category=category,start=int(start),end=int(stop),strandedness=strandedness, frame=frame)
                if new_feature.contig not in genome_annotation:
                    genome_annotation[new_feature.contig]=[new_feature]
                else:
                    genome_annotation[new_feature.contig].append(new_feature)

    if verbose:
        end_time = time.perf_counter()
        execution_time = end_time - start_time
        print(f"\tparsing time: {execution_time:.2f} seconds")
        print(f"\t  * number of contigs: {len(genome_annotation)}")
        print(f"\t  * {featurecategory} features: {count_feature}")
            
    return genome_annotation


def parse_gff3_for_attributes(filepath, feature_of_interest = FeatureCategory.Transcript):
    """
    Parse gff_file to only include ID and attribute information, of only one type of feature (since i usually only need transcripts for this)
    """

    genome_annotation:dict = {}

    with open(filepath, "r") as file:
        linelist = file.readlines()
        
        # determine separator for attributes line
        tail_line = linelist[-10].split("\t")[-1].split(";")[0]
        separator = " "
        if "=" in tail_line:
            separator = "="
        

        for line in linelist:    
            # Skip empty lines or lines starting with a comment character
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            contig,source,category_,start,stop,score,strandedness,frame,attributes_=[c for c in line.split("\t") if len(c)>0]
            category = categorize_string(category_)

            if category != feature_of_interest:
                continue

            # parse all attributes
            attributes={}
            for attr in attributes_.strip().split(";"):
                attr = attr.strip()
                try:
                    key,value=attr.split(separator)[-2:]
                except:
                    continue

                if "," in attr:
                    attributes[key]=value
                else:
                    attributes[key]=value
            
            ## check that ID and Parent are detected correctly
            if "ID" not in attributes:
                    raise RuntimeError(f"no id property found for gene in line: {line}")
        
            feature_id = attributes["ID"]
            genome_annotation[feature_id] = attributes
    
    return genome_annotation


# if __name__ == "__main__":
#     print("1", "contig1", FeatureCategory.Gene, 1, 100, "+", ".", parent_id=None))
#     obtectus_test_native = "/Users/milena/work/native_annot_gff/acanthoscelides_obtectus.gff.isoform_filtered"
#     parse_gff3_general(obtectus_test_native)