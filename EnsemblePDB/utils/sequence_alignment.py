""" psuedo_ensemble.utils.sequence_alignment

Generalized way to easily call pairwise2 with any options desired.

Author:
    Rachael Kretsch (rkretsch@stanford.edu), 
    Siyuan Du, 
    Jacob Parres-Gold

Last edited:
    2021-08-10
"""

from Bio import pairwise2


def specific_pairwise2_align(seqA, seqB, pw_align_opts):
    '''use pairwise2 with any options
    For information on pairwise2 please see https://biopython.org/docs/1.75/api/Bio.pairwise2.html.
    Arguments:
        seqA (str): first sequence to align (ref)
        seqB (str): second sequence to align (seq)
        pw_align_opts (dict): dictionary of alignment options. Keys are:
            "alignment_type": "local" or "global". First, you can give a score for a
            a match and mismatch ("match_score" and "mismatch_score") or
            specify a dictionary that has scores for every possible pair
            ("score_match_dict"). Or you can specify a function for 
            how to score pairing on ref and seq ("score_match_func_ref" 
            and "score_match_func_seq"). Second, you can decide how
            gaps are scored. Either a set open and extend score (
            "gap_open_score" and "gap_extend_score") or you can vary
            how gaps are scored on the reference sequence and input sequence 
            ("gap_open_score_ref","gap_extend_score_ref","gap_open_score_seq",
            "gap_extend_score_seq"). Finally you can define your own
            arbitary functions ("score_gap_func_ref","score_gap_func_seq")
    Returns:
        Alignment objects???
    '''
    # Check all inputrs
    assert (pw_align_opts['alignment_type'] in ["local", "global"]
            ), "alignment_type must be global or local"
    if 'match_score' in pw_align_opts or 'mismatch_score' in pw_align_opts:
        assert (
            'mismatch_score' in pw_align_opts and 'match_score' in pw_align_opts), "Need match score and mismatch score."
        assert ('score_match_dict' not in pw_align_opts and 'score_match_func_ref' not in pw_align_opts and
                'score_match_func_seq' not in pw_align_opts), "Can only either provide a match and mismatch score, a dictionary with match scores, or functions to calculate match scores"
        match_type = 'm'
        match_input = f",{pw_align_opts['match_score']},{pw_align_opts['mismatch_score']}"
    elif 'score_match_dict' in pw_align_opts:
        assert ('score_match_func_ref' not in pw_align_opts and
                'score_match_func_seq' not in pw_align_opts), "Can only either provide a match and mismatch score, a dictionary with match scores, or functions to calculate match scores"
        match_type = 'd'
        match_input = f",{pw_align_opts['score_match_dict']}"
    elif 'score_match_func_ref' in pw_align_opts or 'score_match_func_seq' in pw_align_opts:
        assert('score_match_func_ref' in pw_align_opts and
               'score_match_func_seq' in pw_align_opts), "Need match score function for reference and sequence."
        match_type = 'c'
        match_input = f",{pw_align_opts['score_match_func_ref']},{pw_align_opts['score_match_func_seq']}"
    else:
        match_type = 'x'
        match_input = ""
    if 'gap_open_score' in pw_align_opts or 'gap_extend_score' in pw_align_opts:
        assert ('gap_open_score' in pw_align_opts and 'gap_extend_score' in pw_align_opts), "Need a score for gap opnening and extension."
        assert ('gap_open_score_ref' not in pw_align_opts and 'gap_extend_score_ref' not in pw_align_opts and
                'gap_open_score_seq' not in pw_align_opts and 'gap_extend_score_seq' not in pw_align_opts and
                'score_gap_func_ref' not in pw_align_opts and 'score_gap_func_seq' not in pw_align_opts), "Can only either provide a extend and open gap score for both strands, for each strand individual or function to calculate gap scores."
        gap_type = "s"
        gap_input = f",{pw_align_opts['gap_open_score']},{pw_align_opts['gap_extend_score']}"
    elif ('gap_open_score_ref' in pw_align_opts or 'gap_extend_score_ref' in pw_align_opts or
            'gap_open_score_seq' in pw_align_opts or 'gap_extend_score_seq' in pw_align_opts):
        assert ('gap_open_score_ref' in pw_align_opts and 'gap_extend_score_ref' in pw_align_opts and
                'gap_open_score_seq' in pw_align_opts and 'gap_extend_score_seq' in pw_align_opts), "Need a gap open and gap extend score for reference and sequence."
        assert ('score_gap_func_ref' not in pw_align_opts and 'score_gap_func_seq' not in pw_align_opts), "Can only either provide a extend and open gap score for both strands, for each strand individual or function to calculate gap scores."
        gap_type = "d"
        gap_input = f",{pw_align_opts['gap_open_score_ref']},{pw_align_opts['gap_extend_score_ref']},{pw_align_opts['gap_open_score_seq']},{pw_align_opts['gap_extend_score_seq']}"
    elif 'score_gap_func_ref' in pw_align_opts or 'score_gap_func_seq' in pw_align_opts:
        assert ('score_gap_func_ref' in pw_align_opts and 'score_gap_func_seq' in pw_align_opts), "Need a gap score function for reference and sequence."
        gap_type = "c"
        gap_input = f",{pw_align_opts['score_gap_func_ref']},{pw_align_opts['score_gap_func_seq']}"
    else:
        gap_type = "x"
        gap_input = ""

    # run alignment using correct funciton and inputs
    align = {}
    exec("a= pairwise2.align."+pw_align_opts['alignment_type']+match_type+gap_type +
         "('"+seqA+"','"+seqB+"'"+match_input+gap_input+",one_alignment_only=True)", {'pairwise2': pairwise2}, align)
    return align['a']
