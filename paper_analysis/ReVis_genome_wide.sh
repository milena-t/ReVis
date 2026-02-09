#!/bin/bash -l
#SBATCH -A uppmax2026-1-8
#SBATCH -c 1
#SBATCH -t 1:00:00
#SBATCH -J ReVis_Nvit
#SBATCH -o ReVis_Nvit.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user milena.trabert@ebc.uu.se

source /proj/naiss2023-6-65/Milena/python_venvs/venv/bin/activate

REVIS_PATH=/proj/naiss2023-6-65/Milena/ReVis_paper/ReVis/src/ReVis/
ASSEMBLY_DIR=/proj/naiss2023-6-65/Milena/ReVis_paper/Nvit_analysis/assembly

python3 ${REVIS_PATH}ReVis.py \
    --masker_outfile GCF_009193385.2_Nvit_psr_1.1_genomic_short_headers.fna.out \
    --masker_out_gff GCF_009193385.2_Nvit_psr_1.1_genomic_short_headers.fna.out.gff \
    --out_dir /proj/naiss2023-6-65/Milena/ReVis_paper/Nvit_analysis/ReVis_plots \
    --species_name N_vitripennis \
    --window_length 1e6 \
    --plot_overlap_filtered \
    --verbose \
    --plot