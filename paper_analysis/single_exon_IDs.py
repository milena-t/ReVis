"""
get all single-exon gene IDs and multi-exon gene IDs as separate lists
"""

import parse_gff as gff

def make_lists(annot_path, se_path, me_path):
    annot = gff.parse_gff3_general(annot_path, keep_feature_category=gff.FeatureCategory.Transcript)
    print(len(annot))


if __name__ == "__main__":
    workdir = "/proj/naiss2023-6-65/Milena/ReVis_paper/Brsri_analysis/annotation/"
    annot_path=f"{workdir}GCF_032445375.1_Brsri_v3_genomic_isoform_filtered.gff"
    se_file = f"{workdir}single_exons_list.txt"
    me_file = f"{workdir}multi_exons_list.txt" 
    
    make_lists(annot_path=annot_path)

