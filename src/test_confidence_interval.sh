## do the plotting based on existing tables


U_NAME=miltr339
REP_TABLES=/Users/${U_NAME}/work/PhD_code/PhD_chapter1/data/repeats_tables
SPECIES=B_siliquastri

for SPECIES in A_obtectus A_verrucosus C_chinensis C_maculatus C_septempunctata D_melanogaster D_ponderosae I_luminosus P_pyralis R_ferrugineus T_castaneum T_molitor Z_morio
do 
python3 /Users/${U_NAME}/work/PhD_code/ReVis/src/ReVis/ReVis_transcript_surroundings.py \
    --plot \
    --out_dir /Users/${U_NAME}/work/PhD_code/ReVis/tests \
    --all_before_table ${REP_TABLES}/${SPECIES}_cumulative_repeats_before_all_transcripts.txt \
    --sig_before_table ${REP_TABLES}/${SPECIES}_cumulative_repeats_before_sig_transcripts_90th_GF_size_percentile.txt \
    --all_after_table ${REP_TABLES}/${SPECIES}_cumulative_repeats_after_all_transcripts.txt \
    --sig_after_table ${REP_TABLES}/${SPECIES}_cumulative_repeats_after_sig_transcripts_90th_GF_size_percentile.txt \
    --num_transcripts ${REP_TABLES}/${SPECIES}_transcript_numbers.txt \
    --overlapping_windows \
    --polreg_win_smooth 1 \
    --species_name $SPECIES \
    --verbose
done