#!/bin/bash -l
#SBATCH -A uppmax2026-1-8
#SBATCH -c 20
#SBATCH -t 2-00:00:00
#SBATCH -J repeatmasking_Nvit
#SBATCH -o repeatmasking_Nvit.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user milena.trabert@ebc.uu.se

module load RepeatModeler/2.0.7-foss-2024a RepeatMasker/4.2.1-foss-2024a SAMtools/1.22-GCC-13.3.0

"""
NOTES ON USEAGE:

if the LIBRARIES_DIR exists, then the script will assume that there is a working repeat library in it and skip the library creation step and go straight to masking.
LIBRARIES_DIR should contain identifying information, like a species name.

Good blog post with more detailed explanations: https://darencard.net/blog/2022-07-09-genome-repeat-annotation/

====================================
"""


# if [ $# -lt 2 ]; then
#     echo "Usage: $0 path_to_assembly repeat_library_dir species_identifier_string"
#     echo "you have $#"
#     exit 1
# fi
# 
# ASSEMBLY=$1
# LIBRARIES_DIR=$2
# SPECIES_IDENT=$3

ASSEMBLY="/proj/naiss2023-6-65/Milena/ReVis_paper/Nvit_analysis/assembly/GCF_009193385.2_Nvit_psr_1.1_genomic_short_headers.fna"
LIBRARIES_DIR="/proj/naiss2023-6-65/Milena/ReVis_paper/Nvit_analysis/Nvit_repeat_library"
SPECIES_IDENT=N_vitripennis


## make custom repeat library based on the species assembly
# RepeatModeler uses a NCBI BLASTDB as input to the repeat modeling pipeline, BuildDatabase is a wrapper to make this database for all future steps

#### make the custom repeat libraries with repeatmodeller

if [ -d $LIBRARIES_DIR ]; then
  echo "Directory '$LIBRARIES_DIR' already exists, assume it has a repeat library in it: ${SPECIES_IDENT}_repeats-families.fa"
else
  mkdir -p "$LIBRARIES_DIR"
  BuildDatabase -name "${LIBRARIES_DIR}/${SPECIES_IDENT}_repeats" $ASSEMBLY  # this takes like 15 mins for Cmac
  echo "=====================> build database done"
  RepeatModeler -database "${LIBRARIES_DIR}/${SPECIES_IDENT}_repeats" -threads 20 -LTRStruct  # this takes over a day for Cmac
  echo "=====================> repeatmodeller done"
fi


echo "REPEAT LIBRARY: ${LIBRARIES_DIR}/${SPECIES_IDENT}_repeats-families.fa"
echo "ASSEMBLY: $ASSEMBLY"

echo "  --> truncate assembly fasta headers to max 50 characters because repeatmasker wants that"
head -n 1 $ASSEMBLY

# cut to 49 characters
# echo "sed -i 's/^(>.{48}).*/\1/' $ASSEMBLY"
# sed -E -i 's/^(>.{48}).*/\1/' $ASSEMBLY

# cut after first underscore
echo "'s/^>([^_]*)_.*/>\1/' $ASSEMBLY"
sed -r -i 's/^>([^_]*)_.*/>\1/' $ASSEMBLY
head -n 1 $ASSEMBLY

echo "  --> index the assembly to match the new fasta headers"
samtools faidx $ASSEMBLY # indexed file will be in the same directory as $ASSEMBLY, not working directory
head -n 1 "${ASSEMBLY}.fai"


echo " "
echo "  --> start repeatmasking"
RepeatMasker -lib "${LIBRARIES_DIR}/${SPECIES_IDENT}_repeats-families.fa" -xsmall -s -u -engine ncbi -gff -pa 20 $ASSEMBLY
