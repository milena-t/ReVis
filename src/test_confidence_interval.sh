## do the plotting based on existing tables


U_NAME=miltr339
REP_TABLES=/Users/${U_NAME}/work/PhD_code/PhD_chapter1/data/repeats_tables
SPECIES=A_obtectus

python3 /Users/${U_NAME}/work/PhD_code/ReVis/src/ReVis/ReVis_transcript_surroundings.py \
    --plot \
    --out_dir /Users/${U_NAME}/work/PhD_code/ReVis/tests \
    --all_before_table ${REP_TABLES}/${SPECIES}_cumulative_repeats_before_all_transcripts.txt \
    --sig_before_table ${REP_TABLES}/${SPECIES}_cumulative_repeats_before_sig_transcripts_90th_GF_size_percentile.txt \
    --all_after_table ${REP_TABLES}/${SPECIES}_cumulative_repeats_after_all_transcripts.txt \
    --sig_after_table ${REP_TABLES}/${SPECIES}_cumulative_repeats_after_sig_transcripts_90th_GF_size_percentile.txt \
    --num_transcripts ${REP_TABLES}/${SPECIES}_transcript_numbers.txt \
    --species_name $SPECIES \
    --verbose
##