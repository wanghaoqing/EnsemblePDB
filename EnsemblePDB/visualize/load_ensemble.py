''' pseudo_ensemble.visualize.load_ensemble
Insert explanation
Author:
    Rachael Kretsch (rkretsch@stanford.edu), Siyuan Du, Jacob Parres-Gold
Last edited:
    2020-01-22
Last reviewed:
TODO:
    in subset think about homodimer etc.
    save some file that has the old and new chian include split and combine changes
    # file with
    # pdb_old_chain , pdb_new_chain
    # pdb_old_chain , pdb_new_chain
    # etc.
    add input for directory to save output folders
    For sequence renumber right now its sequential allow to specify??
# TODO is this resistant to the order alignments are put in?
# temp solution naive align first
# TODO can also swap around who is reference structure!
# TODO test if ever has secondary structure visualizing issues after renumbering
# psico.editing.dssp('all') or use Stride
# not editing the sec structure (helix sheet) info as it is in the txt of
# the biopandas object, unclear how pymol maps secondary structure if the pdb
# text is followed first?
# VISUALS
# https://github.com/msabrowser/msabrowser/tree/master/examples
# https://github.com/sanderlab/alignmentviewer
# https://github.com/orangeSi/pymsaploter
# MSA
# https://github.com/benhid/Sequoya
# http://qiime.org/1.8.0/scripts/align_seqs.html
# TODO consider how deals with non-conventional AAs, gaps, and multiconformers
# non-conventional labeled as ? in Biopandas pdb and thus not renumbered
# Gaps
# in renumbering on relabel what is present in the pdb atom definitions only
# multiconformers:
# they are labeled in the alt_loc column and should be preserved here
# hetatoms
# in seq alignment just labeled as - ... in .df["HETATOM"]
# can maybe renumber based on what they are near?
# or if not match add anyway after merge hetatom and atom?
# also don't deal with anisou or ter
'''

#TODO: split chains

from pymol import cmd
from os import path
from pathlib import Path
import pandas as pd

def load_ensemble(directory, Entry_IDs = None,reference_pdb=None, split_chains=False, groups=None):
    '''
    Function that loads every pdb in a folder, reference first.
    Arguments: 
           directory (str): Path to the folder in which these
                               pdbs can be found. Note this folder should
                               not contain any pdbs you do not wish to analyse.
            Entry_IDs (list): list of entry ids to load
            reference_pdb (str): pdb id of reference to be loaded first
            groups = df with at least 2 columns, Entity ID and Subensemble
    '''
    # Checking inputs
    assert (path.isdir(directory)), directory + " is not a folder."

    # Setting up pymol
    cmd.bg_color('white')
    cmd.set('depth_cue', 0)
    cmd.set('ray_trace_fog', 0)

    # print("Using all files, *.pdb, in ", directory)

    # load reference pdb and then load and align other pdbs
    if not Entry_IDs:
        print("Using all files, *.pdb, in ", directory)
        files = list(Path(directory).glob("*.pdb"))
    else:
        print("Using all files in the given Entry ID list")
        Entry_IDs = [x.lower() for x in Entry_IDs]
        files = [x for x in list(Path(directory).glob("*.pdb")) if x[:4].lower() in Entry_IDs]
    assert (files != []), "No pdbs in " + directory


    # print(files)
    if reference_pdb is not None:
        reference_pdb = reference_pdb[:4].lower()+reference_pdb[4:]
        ref_file = Path(directory) / f"{reference_pdb}.pdb"
        assert (path.isfile(ref_file)), f"{ref_file} is not a file."
        cmd.load(ref_file)
        # print(ref_file)
        if split_chains:
            cmd.split_chains(reference_pdb)
            cmd.group(reference_pdb + "_", " ".join(
                [name for name in cmd.get_names() if name[:4] == reference_pdb]))
        files.remove(ref_file)
    for file in files:
        cmd.load(file)
        pdb = file.stem  # split("/")[-1][:-4]
        if split_chains:
            cmd.split_chains(pdb)
            cmd.group(
                pdb + "_", " ".join([name for name in cmd.get_names() if name[:4] == pdb]))

    if groups is not None:
        #JANK
        pdbs = [file.stem for file in files]
        groups = pd.read_csv(groups)
        group_pdbs = groups["Entry ID"].str.lower().tolist()
        extra_entry = groups[~groups["Entry ID"].str.lower().isin(pdbs)][
            "Entry ID"].str.lower().tolist()
        print("Some pdbs in the inputed dataframe are not in the directory", extra_entry)
        groups = groups[groups["Entry ID"].str.lower().isin(pdbs)]
        extra_pdbs = [pdb for pdb in pdbs if not pdb in group_pdbs]
        print("WARNING these pdbs were found in the directory but not in the inputted dataframe", extra_pdbs)
        for subensemble in groups["Subensemble"].unique():
            sub_group = groups[groups["Subensemble"] == subensemble][
                "Entry ID"].str.lower().tolist()
            print(subensemble, sub_group)
            if split_chains:
                length = 5
            else:
                length = 4
            cmd.group(subensemble, " ".join([name for name in cmd.get_names() if name[:4] in sub_group and len(name) == length]))

cmd.extend("load_all", load_ensemble)