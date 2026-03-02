#!/bin/bash -l
#SBATCH -A uppmax2026-1-8
#SBATCH -c 1
#SBATCH -t 2:00:00
#SBATCH -J ReVis_orphans_Aobt
#SBATCH -o ReVis_orphans_Aobt.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user milena.trabert@ebc.uu.se

## don't load the biopython module if using the venv!
source /proj/naiss2023-6-65/Milena/python_venvs/venv/bin/activate

REVIS_PATH=/proj/naiss2023-6-65/Milena/ReVis_paper/ReVis/src/ReVis/
ASSEMBLY_DIR=/proj/naiss2023-6-65/Milena/annotation_pipeline/only_orthodb_annotation/A_obtectus/
ANNOTATION_DIR=/proj/naiss2023-6-65/Milena/annotation_pipeline/only_orthodb_annotation/A_obtectus/braker/

python3 ${REVIS_PATH}ReVis_transcript_surroundings.py \
    --compute_tables_from_list \
    --out_dir /proj/naiss2023-6-65/Milena/ReVis_paper/test_orphans \
    --masker_outfile ${ASSEMBLY_DIR}assembly_genomic.fna.out \
    --annotation_gff ${ANNOTATION_DIR}braker_isoform_filtered.gff \
    --all_list /proj/naiss2023-6-65/Milena/ReVis_paper/test_orphans/Aobt_BUSCOs.txt \
    --sig_list /proj/naiss2023-6-65/Milena/ReVis_paper/test_orphans/Aobt_orphans.txt \
    --species_name A_obtectus \
    --polreg_fourier_denoise \
    --bp 10000 \
    --verbose

    # g52254.t1