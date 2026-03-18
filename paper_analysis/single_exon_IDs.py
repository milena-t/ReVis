"""
get all single-exon gene IDs and multi-exon gene IDs as separate lists
"""

import parse_gff as gff

def make_lists(annot_path, se_path, me_path):
    annot = gff.parse_gff3_general(annot_path)# , keep_feature_category=gff.FeatureCategory.Transcript)
    se_list = []
    me_list = []
    # print(annot["rna-XM_063365854.1"])
    for gene_ID, gene_feature in annot.items():
        if gene_feature.category == gff.FeatureCategory.Transcript:
            if len(gene_feature.child_ids_list) == 1:
                se_list.append(gene_ID)
            elif len(gene_feature.child_ids_list) > 1:
                me_list.append(gene_ID)
    with open(se_path, "w") as se, open(me_path, "w")as me:
        se.write(",".join(se_list))
        me.write(",".join(me_list))

    print(f"outfiles written:\n\t* {se_path}\n\t* {me_path}")


            
            


if __name__ == "__main__":
    workdir = "/proj/naiss2023-6-65/Milena/ReVis_paper/Brsri_analysis/annotation/"
    annot_path=f"{workdir}GCF_032445375.1_Brsri_v3_genomic_isoform_filtered.gff"
    se_file = f"{workdir}single_exons_list.txt"
    me_file = f"{workdir}multi_exons_list.txt" 

    make_lists(annot_path=annot_path, se_path=se_file, me_path=me_file)

