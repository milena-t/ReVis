#!/bin/bash -l
#SBATCH -A uppmax2026-1-8
#SBATCH -c 1
#SBATCH -t 2:00:00
#SBATCH -J ReVis_transcript_Nvit
#SBATCH -o ReVis_transcript_Nvit.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user milena.trabert@ebc.uu.se

## don't load the biopython module if using the venv!
source /proj/naiss2023-6-65/Milena/python_venvs/venv/bin/activate

REVIS_PATH=/proj/naiss2023-6-65/Milena/ReVis_paper/ReVis/src/ReVis/
ASSEMBLY_DIR=/proj/naiss2023-6-65/Milena/ReVis_paper/Brsri_analysis/assembly/
ANNOTATION_DIR=/proj/naiss2023-6-65/Milena/ReVis_paper/Brsri_analysis/annotation/

python3 ${REVIS_PATH}ReVis_transcript_surroundings.py \
    --compute_tables_from_list \
    --out_dir /proj/naiss2023-6-65/Milena/ReVis_paper/Brsri_analysis/ReVis_plots \
    --masker_outfile ${ASSEMBLY_DIR}GCF_032445375.1_Brsri_v3_genomic.fna.out \
    --annotation_gff ${ANNOTATION_DIR}GCF_032445375.1_Brsri_v3_genomic_isoform_filtered.gff \
    --all_list /proj/naiss2023-6-65/Milena/ReVis_paper/Brsri_analysis/Brsri_Revis_enrichment_analysisBrsri_background_transcripts_list.txt \
    --sig_list /proj/naiss2023-6-65/Milena/ReVis_paper/Brsri_analysis/Brsri_Revis_enrichment_analysisBrsri_foreground_transcripts_list.txt \
    --species_name B_rossius_r \
    --polreg_fourier_denoise \
    --bp 10000 \
    --verbose