#!/bin/bash -l
#SBATCH -A uppmax2026-1-8
#SBATCH -c 1
#SBATCH -t 2:00:00
#SBATCH -J filter_for_longest_isoform
#SBATCH -o filter_for_longest_isoform.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user milena.trabert@ebc.uu.se

#filter gff files for the longest isoforms

module load AGAT/1.6.1-GCCcore-13.3.0

# ANNOT_GFF=/proj/coleoptera-genomics-2025/snic2021-6-30/Milena/ReVis_paper/Nvit_analysis/annotation/GCF_009193385.2_Nvit_psr_1.1_genomic.gff
# FILTERED_GTF=/proj/coleoptera-genomics-2025/snic2021-6-30/Milena/ReVis_paper/Nvit_analysis/annotation/GCF_009193385.2_Nvit_psr_1.1_genomic_isoform_filtered.gff
# ASSEMBLY=/proj/coleoptera-genomics-2025/snic2021-6-30/Milena/ReVis_paper/Nvit_analysis/assembly/GCF_009193385.2_Nvit_psr_1.1_genomic_short_headers.fna

ANNOT_GFF=/proj/coleoptera-genomics-2025/snic2021-6-30/Milena/ReVis_paper/Brsri_analysis/annotation/GCF_032445375.1_Brsri_v3_genomic.gff
FILTERED_GTF=/proj/coleoptera-genomics-2025/snic2021-6-30/Milena/ReVis_paper/Brsri_analysis/annotation/GCF_032445375.1_Brsri_v3_genomic_isoform_filtered.gff
ASSEMBLY=/proj/coleoptera-genomics-2025/snic2021-6-30/Milena/ReVis_paper/Brsri_analysis/assembly/GCF_032445375.1_Brsri_v3_genomic.fna.masked

# interactive -A uppmax2026-1-8 -t 5:00:00
perl /proj/coleoptera-genomics-2025/snic2021-6-30/Milena/gene_family_analysis/filter_longest_isoform/agat_sp_keep_longest_isoform.pl --gff $ANNOT_GFF -o $FILTERED_GTF
ANNOT_GFF_ID=${ANNOT_GFF_RAW}_AGAT_ID.gff
agat_sp_manage_IDs.pl --gff $FILTERED_GTF -o $ANNOT_GFF_ID

# module load gffread/0.12.7-GCCcore-13.3.0
# gffread $ANNOT_GFF -g $ASSEMBLY -M -K -o $FILTERED_GTF


