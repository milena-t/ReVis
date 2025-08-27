import parse_gff as gff
import parse_repeats as repeats
import parse_orthogroups as OGs

from tqdm import tqdm


def get_sig_transcripts(orthoDB_orthogroups):
    """
    get a list of all transcripts in the orthoDB dictionary (presumeably significant transcripts)
    """
    all_transcript_IDs = []
    for transcripts_list in orthoDB_orthogroups.values():
        for transcript_id in transcripts_list:
            transcript_id = transcript_id[:-2] # remove the "_1" suffix
            if len(transcript_id)>0:
                all_transcript_IDs.append(transcript_id)
    return all_transcript_IDs


def get_all_repeat_categories(repeats_dict):
    all_contigs = []
    for contig in repeats_dict.keys():
        all_contigs.append(list(set([repeat_instance.repeat_category for repeat_instance in repeats_dict[contig]])))
    all_types = [rep_type for contig_list in all_contigs for rep_type in contig_list]
    all_types = list(set(all_types))
    return(all_types)



def make_cumulative_TE_table(orthogroups_path, n:int, species:str, repeats_annot_path:str, genome_annot_path:str, sig_orthogroups = [], count_transcripts = False, verbose=True):
    """
    make a table for the surrounding n bases upstream and downstream of each gene where each base is a row 
    and each column is the sum of how often this base is annotated as the TE-category across all transcripts

    optionally: provide a sig_orthogroups list to not include all transcripts that are in orthogroups_path
    orthogroups_path can also be a list of transcript IDs to be considered, then reading the orthogroups is skipped entirely
    
    if count_transcripts=True it only returns a list that includes all the transcripts that were included in the computtion
    """
    if type(orthogroups_path) == str:
        orthoDB_orthogroups = OGs.parse_orthogroups_dict(orthogroups_path, sig_orthogroups, species=species)
        all_transcript_IDs = get_sig_transcripts(orthoDB_orthogroups)
        if verbose:
            print(f"{species} --> {len(orthoDB_orthogroups)} orthogroups with {len(all_transcript_IDs)} transcripts (includes orthogroups with no members in {species}")
    elif type(orthogroups_path) == list:
        all_transcript_IDs = orthogroups_path

    print(f"\t*  parse gene-annotation {genome_annot_path}")
    genes_dict = gff.parse_gff3_general(genome_annot_path, verbose=False, keep_feature_category=gff.FeatureCategory.Transcript)
    print(f"\t*  parse repeat-annotation {repeats_annot_path}")
    repeats_dict = repeats.parse_repeats_repeatmasker_outfile(repeats_annot_path, verbose=False)

    repeats_categories = get_all_repeat_categories(repeats_dict=repeats_dict)
    repeats_categories.sort()
    print(f"\t*  all repeat categories: {repeats_categories}")

    before_transcript = { cat : [0]*n for cat in repeats_categories} ## dict with lists for each category from 0 to n, each sub-list covering the category at base n
    after_transcript = { cat : [0]*n for cat in repeats_categories} 
    """
    number of transcript where the base i positions before transcript start is covered by each TE category. 
    Order given by the list in repeats_categories
    {                 0                 n
        category1 : [ 0 , 5, 43, 0, 0, ...],
                      n                2*n
        category2 : [ 0 , 7, 37, 4, 0, ...], 
        ...
    }
    """
    missing_in_annot_transcripts = []
    contigs_with_no_repeats = [] 

    all_transcripts_list = []
    if verbose:
        print(f"\t*  calculate surroundings for {len(all_transcript_IDs)} genes")
    for transcript_id in tqdm(all_transcript_IDs):
        try:
            transcript = genes_dict[transcript_id]
            if count_transcripts:
                all_transcripts_list.append(transcript_id)
                continue
        except:
            raise RuntimeError(f"{transcript_id} can not be found in {genome_annot_path}. Are you using the correct annotation?")
            missing_in_annot_transcripts.append(transcript_id)
            print(f"\ttranscript: {transcript_id}")
            continue
        # start and end of the interval surrounding this transcript
        int_start = transcript.start - n
        int_stop = transcript.end + n
        
        # print(f"\t\t -  {transcript_id} interval from {int_start} to {int_stop}")
        
        ################
        ### start filling first half before the coding region
        
            # collect all repeats that are in the pre-transcript interval
            # this fails if there's no repeats on the transcript.contig, which can happen in very fragmented assemblies
        try:
            repeat_before_transcript = []
            for repeat in repeats_dict[transcript.contig]:
                rep_stop_in_interval = repeat.end < transcript.start and repeat.end > int_start
                rep_start_in_interval = repeat.start < transcript.start and repeat.start > int_start
                rep_longer_than_interval = repeat.start < int_start and repeat.end >int_stop

                if rep_stop_in_interval or rep_start_in_interval or rep_longer_than_interval:
                    repeat_before_transcript.append(repeat)
        except:
            contigs_with_no_repeats.append(transcript.contig)
            # print(f"{transcript.contig} has no repeats on it in {repeats_annot_path}")
            continue
        
        # iterate through every base in the interval
        for index, base in enumerate(range(int_start, transcript.start)):
            for repeat in repeat_before_transcript:
                if base >= repeat.start and base <= repeat.end:
                    # print(f"{repeat.start} < {base} [{index}] < {repeat.end}, {repeat}")
                    before_transcript[repeat.repeat_category][index] += 1
        
        ################
        ### fill out the second half after the coding region, same as above but post-transcript interval
        repeat_after_transcript = [repeat for repeat in repeats_dict[transcript.contig] if (repeat.start < int_stop and repeat.start>transcript.end) or (repeat.end > transcript.end and repeat.end<int_stop) or (repeat.start < transcript.end and repeat.end > int_stop)]
        for index, base in enumerate(range(transcript.end, int_stop)):
            for repeat in repeat_after_transcript:
                if base >= repeat.start and base <= repeat.end:
                    after_transcript[repeat.repeat_category][index] += 1
    
    if len(missing_in_annot_transcripts)>0:
        print(f"{len(missing_in_annot_transcripts)} transcripts (of {len(all_transcript_IDs)} total) not found in annotation and were skipped (in C. maculatus this might be due to the liftover?)")
    if len(contigs_with_no_repeats)>0:
        print(f"{len(contigs_with_no_repeats)} contigs with significant genes but no repeats on them")
    
    if not count_transcripts:
        return before_transcript, after_transcript
    elif count_transcripts:
        return all_transcripts_list


