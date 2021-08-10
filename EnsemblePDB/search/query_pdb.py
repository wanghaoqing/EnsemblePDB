""" pseudo_ensemble.search.query_pdb

TODO DESCRIPTION

Authors:
    Rachael Kretsch (rkretsch@stanford.edu), 
    Siyuan Du, 
    Jacob Parres-Gold

Last edited:
    2021-08-10
"""

from pathlib import Path

import pandas as pd
from Bio.pairwise2 import format_alignment
import pypdb

from EnsemblePDB.utils import file_management, table_utils, biopandas_utils, sequence_alignment


def query_by_ref_pdb(reference_pdb, reference_chains, label='MyProtein',
                     macromolecule_name=None, seq_id=0.95,
                     xray_only=True, evalue=10,
                     pairwise_align_options={
                         'alignment_type': "local", 'gap_open_score': -0.5, 'gap_extend_score': -0.1},
                     output_dir=".", min_length=20):
    '''
    Will take sequence(s) or reference pdb chain(s) and will search 
    the PDB for structures of similar sequences. 
    Args: 
        reference_pdb (str): pdb id ('XXXX') to be the reference.
        reference_chains (list of str): list of the chains of the reference
            structure to use to seqeunce search.
        label (str): Label to add to end of file names
      (optional)
        macromolecule_name (str): name to search pdb by. May add structures to list 
            not found by sequence search, or identify pdb not
            labeled with the desired name. 
            Reported in column: "From macromolecular search"
        seq_id (float): if specified uses this number as the minimum percent identity 
            of sequences to be added instead of calculating from max_mut
        xray_only (bool): only return x-ray crystallographic structures
            {default: True}
        evalue (float): 
        pairwise_align_options, {default: {'alignment_type':"local",'gap_open_score':-0.5, 'gap_extend_score':-0.1}}
            Options for paiwise sequence alignment scoring.
            Can score alignment anyway with pairwise2:
            https://biopython.org/docs/1.75/api/Bio.pairwise2.html
            Dictionary keys can be:
                alignment_type="local", match_score, mismatch_score, score_match_dict,
                score_match_func_ref, score_match_func_seq,
                gap_open_score, gap_extend_score,
                gap_open_score_ref, gap_extend_score_ref, gap_open_score_seq, gap_extend_score_seq,
                score_gap_func_ref, score_gap_func_seq,
        output_dir (string): if not working in data directory you can specify where
                        to save output {default: "."}
    Returns: 
        Pandas dataframe: summary output of the similar pdb found, 
            including the chains that aligned with to 
            reference and additionally alignment to
            reference for all chains.
        saves the returned data frame to seq_search_output*.csv
        saves the raw output from PDB sequence search to 
        sequence_search_summary*.txt and sequence_search_summary*.fasta
    '''

    ############################## Check inputs ###############################
    assert (type(reference_pdb) == str), "reference_pdb must be a string."
    assert (len(reference_pdb) ==
            4), "Check your reference_pdb, pdbs are 4 characters in length."
    assert (type(reference_chains) ==
            list), "reference_pdb_chains must be a list, eg ['A']"
    assert (all(isinstance(s, str)
                for s in reference_chains)), "Chains must be strings."
    print(f"A reference pdb was inputted. Will search PDB for structures"
          f" that align to sequence(s) of {reference_pdb} chains "
          f"{reference_chains}.")
    if macromolecule_name is not None:
        assert (type(macromolecule_name) ==
                str), "macromolecule_name must be a string."
    if seq_id is not None:
        assert (seq_id >= 0 and seq_id <=
                1), "seq_id must be a percentage (0-1)"
    assert (type(xray_only) == bool), "xray_only must be a bool."
    assert (isinstance(evalue, (int, float))
            ), "evalue should be a number."
    assert(type(pairwise_align_options) ==
           dict), "pairwise_align_options must be a dictionary."
    assert(Path(output_dir).is_dir()), "Output directory must exist"
    ############################## Check inputs ###############################

    # get list of all x-ray structures
    if not xray_only:
        print("WARNING: obtaining all structures, including structures from "
              "methods other than X-ray crystallography.")
        xray_pdbs = []
    else:
        print("Looking for X-ray crystallographic structures.")
        xray_pdbs = pypdb.Query('X-RAY DIFFRACTION',
                                query_type='ExpTypeQuery').search()
        assert (xray_pdbs is not None), "PDB search failed, try again."

    # if macrmolecule_name specified get list of pdbs with macromolecular name
    if macromolecule_name is not None:
        print(f"Looking for structures with name {macromolecule_name}")

        search_operator = pypdb.clients.search.operators.text_operators.ContainsPhraseOperator(
            value=macromolecule_name, attribute="rcsb_polymer_entity.rcsb_macromolecular_names_combined.name")
        return_type = pypdb.clients.search.search_client.ReturnType.ENTRY
        name_pdbs = pypdb.perform_search(search_operator, return_type)
        assert (name_pdbs is not None), "PDB macromolecule search failed, try again"

        # only keep structures that are x-ray structures
        if xray_only:
            pdb_list = [pdb for pdb in xray_pdbs if pdb in name_pdbs]
        else:
            pdb_list = name_pdbs

        # write dataframe
        macro_output = pd.DataFrame(pdb_list, columns=['Entry ID'])
        macro_output["From macromolecular search"] = True
    else:
        macro_output = None
        pdb_list = xray_pdbs

    # get the sequence of the reference chain(s)
    print(f"Extracting sequence(s) of {reference_pdb}.")
    structure = biopandas_utils.get_pdb_struct(reference_pdb)
    sequences = []
    for ref_chain in reference_chains:
        sequences.append(biopandas_utils.get_chain_seq(structure, ref_chain))
    if len(sequences) > 1:
        i = 0
        while i < len(sequences):
            if sequences[i] in sequences[:i]+sequences[i+1:]:
                print(
                    "WARNING: some of your sequences are redudant! You will get redudant information!")
            i += 1

    # search the pdb for the sequence and return those with >=seq_id
    seq_chain_list = None
    for i, seq in enumerate(sequences):
        search_sum_file = file_management.get_nonexistant_file(Path(output_dir, f'sequence_search_summary_{reference_chains[i]}.txt'))
        print(f"Saving sequence search raw output in {search_sum_file}.")
        temp_list = search_seq(seq, seq_id, f"{search_sum_file}", reference_chains[i], evalue=evalue)

        if seq_chain_list is None:
            seq_chain_list = temp_list
        else:
            seq_chain_list = seq_chain_list.merge(temp_list, how='outer',
                                                  on='Entry ID')

    # only take x-ray structures
    if xray_only:
        output = seq_chain_list[seq_chain_list['Entry ID'].isin(pdb_list)]
    else:
        output = seq_chain_list

    if macro_output is not None:
        output = output.merge(macro_output, how="outer", on="Entry ID")

    # get all chains and there alignments
    print("Adding all chain alignment information. For large datasets this may"
          " take some time.")
    output = get_all_chain_alignments_for_table_multi(seq_id,
                                                      output, sequences, reference_chains, pairwise_align_options, output_dir, min_length)
    output = output.applymap(table_utils.delete_all_one_value)

    # save and return output
    output_file = file_management.get_nonexistant_file(
        Path(output_dir, "seq_search_output_"+label+".csv"))

    print(f"Saved output to {output_file}. To obtain PDB summary info, run "
          f"pseudo_ensemble.search.get_data.get_pdb_data('{output_file}')")
    output.to_csv(output_file, index=False)
    return output


def query_by_sequence(sequences,
                      macromolecule_name=None,  label='MyProtein', seq_id=0.95,
                      xray_only=True, evalue=10,
                      pairwise_align_options={
                          'alignment_type': "local", 'gap_open_score': -0.5, 'gap_extend_score': -0.1},
                      output_dir=".", min_length=20):
    '''
    Will take sequence(s) or reference pdb chain(s) and will search 
    the PDB for structures of similar sequences. 
    Args: 
        sequences (list of str): list of sequence(s) to search on.
      (optional)
        macromolecule_name (str): name to search pdb by. May add structures to list 
            not found by sequence search, or identify pdb not
            labeled with the desired name. 
            Reported in column: "From macromolecular search"
        seq_id (float): if specified uses this number as the minimum percent identity 
            of sequences to be added instead of calculating from max_mut
        xray_only (bool): only return x-ray crystallographic structures
            {default: True}
        evalue (float): 
        pairwise_align_options, {default: {'alignment_type':"local",'gap_open_score':-0.5, 'gap_extend_score':-0.1}}
            Options for paiwise sequence alignment scoring.
            Can score alignment anyway with pairwise2:
            https://biopython.org/docs/1.75/api/Bio.pairwise2.html
            Dictionary keys can be:
                alignment_type="local", match_score, mismatch_score, score_match_dict,
                score_match_func_ref, score_match_func_seq,
                gap_open_score, gap_extend_score,
                gap_open_score_ref, gap_extend_score_ref, gap_open_score_seq, gap_extend_score_seq,
                score_gap_func_ref, score_gap_func_seq,
        output_dir, if not working in data directory you can specify where
                        to save output {default: "."}
    Returns: 
        Pandas dataframe: summary output of the similar pdb found, 
            including the chains that aligned with to 
            reference and additionally alignment to
            reference for all chains.
        saves the returned data frame to seq_search_output*.csv
        saves the raw output from PDB sequence search to 
        sequence_search_summary*.txt and sequence_search_summary*.fasta
    '''

    ############################## Check inputs ###############################
    assert (type(sequences) ==
            list), "Sequences must be a list even if with one element."
    assert (all(isinstance(s, str)
                for s in sequences)), "Sequences must be strings."
    print(f"Sequence(s)inputted. Will search PDB for "
          f"{len(sequences)} sequence(s).")
    if macromolecule_name is not None:
        assert (type(macromolecule_name) ==
                str), "macromolecule_name must be a string."
    if seq_id is not None:
        assert (seq_id >= 0 and seq_id <=
                1), "seq_id must be a percentage (0-1)"
    assert (type(xray_only) == bool), "xray_only must be a bool."
    assert (isinstance(evalue, (int, float))
            ), "evalue should be a number."
    assert(type(pairwise_align_options) ==
           dict), "pairwise_align_options must be a dictionary."
    assert(Path(output_dir).is_dir()), "Output directory must exist"
    ############################## Check inputs ###############################

    # get list of all x-ray structures
    if not xray_only:
        print("WARNING: obtaining all structures, including structures from "
              "methods other than X-ray crystallography.")
        xray_pdbs = []
    else:
        print("Looking for X-ray crystallographic structures.")
        xray_pdbs = pypdb.Query('X-RAY DIFFRACTION',
                                query_type='ExpTypeQuery').search()
        assert (xray_pdbs is not None), "PDB search failed, try again."

    # if macrmolecule_name specified get list of pdbs with macromolecular name
    if macromolecule_name is not None:
        print(f"Looking for structures with name {macromolecule_name}")

        search_operator = pypdb.text_operators.ContainsPhraseOperator(
            value=macromolecule_name, attribute="rcsb_polymer_entity.rcsb_macromolecular_names_combined.name")
        return_type = pypdb.ReturnType.ENTRY
        name_pdbs = pypdb.perform_search(search_operator, return_type)
        assert (name_pdbs is not None), "PDB macromolecule search failed, try again"

        # only keep structures that are x-ray structures
        if xray_only:
            pdb_list = [pdb for pdb in xray_pdbs if pdb in name_pdbs]
        else:
            pdb_list = name_pdbs

        # write dataframe
        macro_output = pd.DataFrame(pdb_list, columns=['Entry ID'])
        macro_output["From macromolecular search"] = True
    else:
        macro_output = None
        pdb_list = xray_pdbs

    # get chain names for sequence only queries
    reference_chains = []
    for i in range(len(sequences)):
        reference_chains.append("inputted_seq_" + str(i))
    if len(sequences) > 1:
        i = 0
        while i < len(sequences):
            if sequences[i] in sequences[:i]+sequences[i+1:]:
                print(
                    "WARNING: some of your sequences are redudant! You will get redudant information!")
            i += 1

    # search the pdb for the sequence and return those with >= seqid
    seq_chain_list = None
    for i, seq in enumerate(sequences):
        search_sum_file = file_management.get_nonexistant_file(Path(output_dir, f'sequence_search_summary_{reference_chains[i]}.txt'))
        print(f"Saving sequence search raw output in {search_sum_file}.")
        temp_list = search_seq(seq, seq_id, f"{search_sum_file}", reference_chains[i], evalue=evalue)

        if seq_chain_list is None:
            seq_chain_list = temp_list
        else:
            seq_chain_list = seq_chain_list.merge(temp_list, how='outer',
                                                  on='Entry ID')

    # only take x-ray structures
    if xray_only:
        output = seq_chain_list[seq_chain_list['Entry ID'].isin(pdb_list)]
    else:
        output = seq_chain_list

    if macro_output is not None:
        output = output.merge(macro_output, how="outer", on="Entry ID")

    # get all chains and there alignments
    print("Adding all chain alignment information. For large datasets this may"
          " take some time.")
    output = get_all_chain_alignments_for_table_multi(seq_id,
                                                      output, sequences, reference_chains, pairwise_align_options, output_dir, min_length)
    output = output.applymap(table_utils.delete_all_one_value)

    # save and return output
    output_file = file_management.get_nonexistant_file(
        Path(output_dir, "seq_search_output_"+label+".csv"))

    print(f"Saved output to {output_file}. To obtain PDB summary info, run "
          f"pseudo_ensemble.search_pdb.get_pdb_data('{output_file}')")
    output.to_csv(output_file, index=False)
    return output

###############################################################################
# Helper functions
###############################################################################


def search_seq(sequence, seq_id, search_output_file, chain_name, evalue=10):
    '''
    Takes in a sequence and a threshold for sequence identity and return a
    table with sequences that reach the threshold with their pdb_id, chaiin
    that matches, and alignment information.
    Arguments: 
        sequence (str): the sequence to search for
        seq_id_min (float): the threshold sequence idenity to include
        chain_name (str): name of reference chain, used only for labeling
    Returns: 
        Pandas Dataframe: with the pdb id, alignment information including
            sequence and the chain(s) that matches
    '''

    search_operator = pypdb.clients.search.operators.sequence_operators.SequenceOperator(sequence=sequence,
                                                                                         evalue_cutoff=evalue, identity_cutoff=seq_id)
    results = pypdb.search_client.perform_search_with_graph(
        query_object=search_operator, return_type=pypdb.search_client.ReturnType.POLYMER_ENTITY, return_raw_json_dict=True)

    cols = ['Entry ID', 'sequence_identity', 'evalue', 'bitscore', 'alignment_length', 'mismatches', 'gaps_opened',
            'query_beg', 'query_end', 'subject_beg', 'subject_end', 'query_length', 'subject_length', 'query_aligned_seq', 'subject_aligned_seq']
    data = pd.DataFrame(columns=cols)
    for hit in results['result_set']:
        match = hit['services'][0]["nodes"][0]["match_context"][0]
        match["Entry ID"] = hit['identifier'][:4]
        data = data.append(match, ignore_index=True)

    # Write out search results
    f = open(search_output_file, "w")
    f.write(str(results))
    f.close()

    # Make data more readable and return
    pretty_col = {"Entry ID": "Entry ID",
                  "sequence_identity": f'Search ref {chain_name}: sequence_identity',
                  "evalue": f'Search ref {chain_name}: evalue',
                  "bitscore": f'Search ref {chain_name}: bitscore',
                  "alignment_length": f'Search ref {chain_name}: alignment_length',
                  "mismatches": f'Search ref {chain_name}: mismatches',
                  "gaps_opened": f'Search ref {chain_name}: gaps_opened',
                  "query_beg": f'Search ref {chain_name}: query_beg',
                  "query_end": f'Search ref {chain_name}: query_end',
                  "subject_beg": f'Search ref {chain_name}: subject_beg',
                  "subject_end": f'Search ref {chain_name}: subject_end',
                  "query_length": f'Search ref {chain_name}: query_length',
                  "subject_length": f'Search ref {chain_name}: subject_length',
                  "query_aligned_seq": f'Search ref {chain_name}: query_aligned_seq',
                  "subject_aligned_seq": f'Search ref {chain_name}: subject_aligned_seq'}
    data = data.rename(columns=pretty_col)
    data = data.astype(str)
    data = data.groupby('Entry ID').agg(" ~ ".join).reset_index()
    return data


def get_all_chain_alignments_for_table_multi(seq_id, input_table, reference_seqs, chain_names, pairwise_align_options, output_dir, min_length):
    '''
    Gets alignments for all matching chains of all pdbs in table, given a list
    of reference_seqs and their corresponding reference chains.
    Chains that do not align to any chain listed in Non_aligned_chains
    '''
    output = input_table.copy()
    never_aligned = None
    all_chains = None

    # For each reference chain get all chain alignments
    for seq, chain in zip(reference_seqs, chain_names):
        search_sum_file = file_management.get_nonexistant_file(Path(output_dir, f'sequence_search_summary_{chain}.fasta'))
        output, non_aligned = get_all_chain_alignments_for_table(f"{search_sum_file}", seq_id, output, seq, chain, pairwise_align_options, min_length)

        # Chains not aligned to any of the references make it to never_aligned
        if never_aligned is None:
            never_aligned = pd.Series(non_aligned)
        else:
            never_aligned = never_aligned.combine(
                pd.Series(non_aligned), func=table_utils.get_overlap)
        if all_chains is None:
            all_chains = output[f"Align ref {chain}: order of chains"]
        else:
            all_chains = all_chains.combine(output[f"Align ref {chain}: order of chains"], func=table_utils.get_union)

    # Add column for chains that are not aligned to any refernce
    output["Align: Non-aligned chains"] = never_aligned
    output['Chain IDs'] = all_chains
    return output


def get_all_chain_alignments_for_table(search_output_fasta, seq_id, input_table, reference_seq, chain_name, pairwise_align_options, min_length):
    '''
    Given a table with pdbs, adds columns indicating how each chain of each
    pdb aligns to the given reference sequence.
    '''

    # COMMENT: just a reminder to check later code using this column name as key, last time I checked the names did not match, ignore if fixed
    chains = f'Align ref {chain_name}: order of chains'
    seq_align_col = f'Align ref {chain_name}: formatted alignment'
    ref_seq_aligned_col = f"Align ref {chain_name}: aligned ref sequence"
    aligned_col = f"Align ref {chain_name}: aligned sequence"
    num_matches_col = f"Align ref {chain_name}: number of AA matches to ref"
    chain_al = f"Align ref {chain_name}: chains aligned"
    ali_score = f"Align ref {chain_name}: pairwise align score"

    # For each row get all relevant chain alignments
    output = input_table.copy()
    output[chains], output[seq_align_col], output[ref_seq_aligned_col], output[aligned_col], output[ali_score], output[
        num_matches_col], output[chain_al], non_aligned = zip(*output['Entry ID'].map(lambda pdb: get_all_chain_alignments(search_output_fasta, seq_id, pdb, reference_seq, pairwise_align_options, min_length)))
    return output, non_aligned

# TODO want to extract the pdb website table into my own dictionary set to interegate it more....


def get_all_chain_alignments(search_output_fasta, seq_id, pdbid, reference_seq, pairwise_align_options, min_length):
    """
    given a pdb id and a reference sequence, returns the alignment
    of every chain of that pdb with the reference sequence
    - chain names
    - formatted alignments output
    - alignments of the reference sequence
    - alignments of the sequences
    - number of AA matches to the refrence chain
    """

    # Get sequence of the structure
    struct = biopandas_utils.get_pdb_struct(pdbid)
    if struct is None:
        return "", "", "", "", "", "", "", "", ""
    seq = struct.amino3to1()

    # set-up data-collection
    chain, output, ref_aligned, aligned, align_score, num_matches, aligned_chain, non_aligned = [
    ], [], [], [], [], [], [], []

    # For each chain if it is a aligned chain do pairwise
    # sequence alignment with reference
    for chain_id in seq['chain_id'].unique():

        # Get the chains sequence
        seqB = ''.join(
            seq.loc[seq['chain_id'] == chain_id, 'residue_name'])

        # Align according to inputted params
        alignments = sequence_alignment.specific_pairwise2_align(
            reference_seq, seqB, pairwise_align_options)

        # If alignment successful and to output
        if alignments == []:
            output.append("")
            ref_aligned.append("")
            aligned.append("")
            align_score.append("")
            num_matches.append("")
            non_aligned.append(chain_id)
            chain.append(chain_id)
            continue
        else:
            output.append(format_alignment(
                *alignments[0], full_sequences=True))
            ref_aligned.append(alignments[0].seqA)
            aligned.append(alignments[0].seqB)
            align_score.append(alignments[0].score)
            matches = format_alignment(
                *alignments[0], full_sequences=True).count("|")
            num_matches.append(matches)
            if (matches/min(len(reference_seq), max(min_length, len(seqB)))) >= seq_id:
                aligned_chain.append(chain_id)
                f = open(search_output_fasta, "a")
                f.write(f">{pdbid}_{chain_id}\n{seqB}\n")
                f.close()
            else:
                non_aligned.append(chain_id)
            chain.append(chain_id)

    # Format output and return
    chain = table_utils.get_table_entry_from_list(chain)
    output = table_utils.get_table_entry_from_list(output)
    ref_aligned = table_utils.get_table_entry_from_list(ref_aligned)
    aligned = table_utils.get_table_entry_from_list(aligned)
    align_score = table_utils.get_table_entry_from_list(align_score)
    num_matches = table_utils.get_table_entry_from_list(num_matches)
    aligned_chain = table_utils.get_table_entry_from_list(aligned_chain)
    non_aligned = table_utils.get_table_entry_from_list(non_aligned)

    return chain, output, ref_aligned, aligned, align_score, num_matches, aligned_chain, non_aligned
