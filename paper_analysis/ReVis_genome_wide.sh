#!/bin/bash -l
#SBATCH -A uppmax2026-1-8
#SBATCH -c 1
#SBATCH -t 1:00:00
#SBATCH -J ReVis_Broc
#SBATCH -o ReVis_Broc.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user milena.trabert@ebc.uu.se

## don't load the biopython module if using the venv!
source /proj/coleoptera-genomics-2025/snic2021-6-30/Milena/python_venvs/venv/bin/activate

REVIS_PATH=/proj/coleoptera-genomics-2025/snic2021-6-30/Milena/ReVis_paper/ReVis/src/ReVis/
ASSEMBLY_DIR=/proj/coleoptera-genomics-2025/snic2021-6-30/Milena/ReVis_paper/Brsri_analysis/assembly/
ANNOTATION_DIR=/proj/coleoptera-genomics-2025/snic2021-6-30/Milena/ReVis_paper/Brsri_analysis/annotation/

python3 ${REVIS_PATH}ReVis.py \
    --masker_outfile ${ASSEMBLY_DIR}GCF_032445375.1_Brsri_v3_genomic.fna.out \
    --masker_out_gff ${ASSEMBLY_DIR}GCF_032445375.1_Brsri_v3_genomic.fna.out.gff \
    --annotation_gff ${ANNOTATION_DIR}GCF_032445375.1_Brsri_v3_genomic_isoform_filtered_AGAT_cleaned.gff \
    --merge_gene_windows 2 \
    --out_dir //proj/coleoptera-genomics-2025/snic2021-6-30/Milena/ReVis_paper/Brsri_analysis/ReVis_plots \
    --species_name B_roccius_r \
    --window_length 5e6 \
    --plot_overlap_filtered \
    --verbose \
    --plot