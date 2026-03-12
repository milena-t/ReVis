#!/bin/bash -l
#SBATCH -A uppmax2026-1-8
#SBATCH -c 1
#SBATCH -t 2:00:00
#SBATCH -J ReVis_orphans_Cmac
#SBATCH -o ReVis_orphans_Cmac.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user milena.trabert@ebc.uu.se

## don't load the biopython module if using the venv!
source /proj/naiss2023-6-65/Milena/python_venvs/venv/bin/activate

REVIS_PATH=/proj/naiss2023-6-65/Milena/ReVis_paper/ReVis/src/ReVis/
ASSEMBLY_DIR=/proj/naiss2023-6-65/Milena/annotation_pipeline/only_orthodb_annotation/C_maculatus/
ANNOTATION_DIR=/proj/naiss2023-6-65/Milena/annotation_pipeline/only_orthodb_annotation/C_maculatus/braker/

## make busco list with awk
# awk -F'\t' '!/^#/{v=$3; sub(/^C_maculatus_/,"",v); sub(/_1$/,"",v); printf "%s%s",(c++?",":""),v} END{print ""}' /proj/naiss2023-6-65/Milena/annotation_pipeline/busco/busco2/BUSCO_orthoDB_annotation/C_maculatus_filtered_proteinfasta.fa/run_arthropoda_odb10/full_table.tsv > /proj/naiss2023-6-65/Milena/ReVis_paper/test_orphans/Cmac_BUSCOs.txt
## make orphan list with awk
# awk '!/^#/{sub(/^C_maculatus_/,""); sub(/_1$/,""); printf "%s%s",(c++?",":""),$0} END{print ""}' /proj/naiss2023-6-65/Ashwath/mas_thesis/results/final_orphans_2.0/C_maculatus_e1e-1_cov0_bits0/C_maculatus_true_orphans.ids.txt > /proj/naiss2023-6-65/Milena/ReVis_paper/test_orphans/Cmac_orphans.txt


python3 ${REVIS_PATH}ReVis_transcript_surroundings.py \
    --compute_tables_from_list \
    --out_dir /proj/naiss2023-6-65/Milena/ReVis_paper/test_orphans \
    --masker_outfile ${ASSEMBLY_DIR}assembly_genomic.fna.out \
    --annotation_gff ${ANNOTATION_DIR}braker_isoform_filtered.gff \
    --all_list /proj/naiss2023-6-65/Milena/ReVis_paper/test_orphans/Cmac_BUSCOs.txt \
    --sig_list /proj/naiss2023-6-65/Milena/ReVis_paper/test_orphans/Cmac_orphans.txt \
    --species_name C_maculatus \
    --polreg_fourier_denoise \
    --bp 50000 \
    --verbose
