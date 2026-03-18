"""
get all single-exon gene IDs and multi-exon gene IDs as separate lists
"""

import parse_gff as gff

def make_lists(annot_path):
    annot = gff.parse_gff3_general(annot_path, keep_feature_category=gff.FeatureCategory.Transcript)
    print(len(annot))


if __name__ == "__main__":
    annot_path="/proj/naiss2023-6-65/Milena/ReVis_paper/Brsri_analysis/annotation/GCF_032445375.1_Brsri_v3_genomic_isoform_filtered.gff"


