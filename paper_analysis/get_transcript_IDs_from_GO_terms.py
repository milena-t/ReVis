"""
get a list of transcript IDs that are assigned GO terms related to toxins or olfactory function to do the ReVis enrichment analysis
the GO terms are in Nvit_gene_ontology.gaf, and the gene IDs listed there correspond to Dbxref in the gff attributes column
"""

import parse_gff as gff
import pandas as pd

def read_GO_terms_dict(go_file_path:str):
    """
    read the GO terms table into a dictionary
    """
    with open(go_file_path, "r") as go_file:
        out_dict = {}
        for line in go_file.readlines():
            GO_term,function = line.strip().split(",")
            out_dict[GO_term] = function
    return out_dict


def get_Dbxref_IDs_of_GO_terms(gene_ontology_path:str, go_dict:dict, add_prefix = "GeneID:"):
    """
    get a list of the Dbxref gene IDs that correspond to our GO terms of interest from the gene ontology file
    in the anntoatio the id is parsed with a prefix, add it here for the string matching later
    """
    header = ["DB","GeneID","Symbol","Qualifier","GO_ID","Reference","Evidence_Code","With,From","Aspect","Gene_Name","Gene_Synonym","Type","Taxon","Date","Assigned_By","Annot_Ext","Gene_Product_Form_ID"]
    gene_ontology_df = pd.read_csv(gene_ontology_path, sep="\t", comment="!", names=header)
    
    go_gene_IDs = {go : [] for go in go_dict.keys()}
    for GO_term in go_gene_IDs.keys():
        print(f"{GO_term} --->  {go_dict[GO_term]}")
        go_genes_df = gene_ontology_df.loc[gene_ontology_df["GO_ID"] == GO_term]
        # print(go_genes_df)
        go_geneIDs_list = go_genes_df["GeneID"].values.tolist()
        len_before = len(go_geneIDs_list)
        # only unique IDs
        go_geneIDs_list = list(set(go_geneIDs_list))
        go_geneIDs_list = [f"{add_prefix}{geneID}" for geneID in go_geneIDs_list] 
        len_after = len(go_geneIDs_list)
        print(f"{len_before} occurences, {len_after} unique IDs: {go_geneIDs_list[:3]}...")
        go_gene_IDs[GO_term] = go_geneIDs_list
    return go_gene_IDs


def find_transcriptID_from_Dbxref(annotation:dict, Dbxref_ID:str):
    """
    find the transcript ID that corresponds to the Dbxref_ID in the gene ontology file that matches the GO term of interest
    """
    for tr_id, feature in annotation.items():
        if feature.other == Dbxref_ID:
            return tr_id


def get_trIDs_from_Dbxref(annotation:dict, Dbxref_list:str):
    """
    get all the transcript IDs that match the Dbxref IDs
    """
    return [find_transcriptID_from_Dbxref(annotation, Dbxref_ID=Dbxref) for Dbxref in Dbxref_list]



if __name__ == "__main__":

    Gene_ontology_file = "Nvit_gene_ontology.gaf"
    annotation_file = "Nvit_genomic.gff"
    outfile_transcripts = "Nvit_GO_smell_transcripts_list.txt"
    outfile_transcripts_all = "Nvit_foreground_transcripts_list.txt"
    outfile_transcripts_bg = "Nvit_background_transcripts_list.txt"

    GO_terms = {
        "toxins" : "GO_terms_toxins.txt",
        "smell" : "GO_terms_smell.txt"
    }

    # go_toxins = read_GO_terms_dict(GO_terms["toxins"])
    # not really any GO terms related to toxins, so I am only using the olfactory
    go_smell = read_GO_terms_dict(GO_terms["smell"])
    IDs_smell = get_Dbxref_IDs_of_GO_terms(Gene_ontology_file, go_smell)

    annotation = gff.parse_gff3_general(annotation_file, other_attribute_string="Dbxref")
    annot_transcripts = {id_str : feature for id_str, feature in annotation.items() if feature.category == gff.FeatureCategory.Transcript}
    print(f"{len(annotation)} features annotated in total")
    print(f"{len(annot_transcripts)} transcripts annotated")
    
    outfile_string_all = {GO_term : "" for GO_term in IDs_smell.keys()}
    with open(outfile_transcripts, "w") as outfile:
        # get transcript IDs for all the Dbxref IDs
        for GO_term, Dbxref_list in IDs_smell.items():
            print(f">>>>>> {GO_term} ({go_smell[GO_term]}): {len(Dbxref_list)} genes")
            tr_IDs = get_trIDs_from_Dbxref(annot_transcripts, Dbxref_list=Dbxref_list)
            tr_IDs_string = ",".join(tr_IDs)
            outfile_string_all[GO_term] = tr_IDs
            outfile.write(f"{GO_term}:{tr_IDs_string}\n")
    
    all_list_tr_ID = []
    for GO_term, list_str in outfile_string_all.items():
        all_list_tr_ID.extend(list_str)
    unique_list_tr_ID = list(set(all_list_tr_ID))
    print(f"\t{len(all_list_tr_ID)} gene IDs accessed in total")
    print(f"\t{len(unique_list_tr_ID)} of them are unique")

    with open(outfile_transcripts_all, "w") as outfile:
        outfile.write(",".join(unique_list_tr_ID))

    ## get background transcripts
    background_annot = annot_transcripts
    for transcript_id in unique_list_tr_ID:
        background_annot.pop(transcript_id)
    print(f"{len(background_annot)} transcripts after fg transcripts are removed")
    with open(outfile_transcripts_bg, "w") as outfile:
        background_IDs = list(background_annot.keys())
        outfile.write(",".join(background_IDs))

