""" psuedo_ensemble.utils.biopandas_utils

Basic functions in manipulating pdb in biopandas format

Author:
    Rachael Kretsch (rkretsch@stanford.edu), 
    Siyuan Du, 
    Jacob Parres-Gold

Last edited:
    2021-08-10
"""

from biopandas.pdb import PandasPdb
from pathlib import Path
import Bio
from Bio import SeqUtils
from numpy.lib.twodim_base import diagflat
import pandas as pd
import numpy as np
from tqdm import tqdm


def get_pdb_struct(reference_pdb):
    '''
    Given a 4-character pdb id fetches the Biopandas.pdb
    Arguments:
        reference_pdb (str): pdb id "XXXX"
    Returns:
        Biopandas.pdb object
    Warnings:
        Returns none if pdb id does not exist, pdb id exists but strcuture has
        not been published, or if connection is poor and it times out. If you
        know structure exists try again or at worst download pdb yourself!
    '''
    try:
        struct = PandasPdb()
        struct.fetch_pdb(reference_pdb)
    except:
        print("Could not retrieve PDB of ", reference_pdb,
              ": this may be an entry that has not published structures yet.")
        return None
    return struct


def get_PandasPDB(file):
    struct = PandasPdb()
    struct.read_pdb(str(file))
    return struct


def get_PandasPDBs(path):
    '''
    Get Pandas PDB object for a set of PDB structures (with aligned coordinates).

    Arguments: 
        path, a directory with aligned PDB files (output from PyMol)
    Returns: 
        a pandas dataframe with 'Entry ID' and 'PandasPDB'
    '''
    # print('Retrieving PDB files...')

    filenames = list(Path(path).glob('*.pdb'))
    structures = []
    IDs = []
    for file in tqdm(filenames, total=len(filenames), desc='Reading PDB files'):
        structures.append(get_PandasPDB(file))
        IDs.append(file.stem)
    data = {'Entry ID': IDs, 'PandasPDB': structures}
    PDB_df = pd.DataFrame(data=data)
    return PDB_df


def get_chain_seq(struct, chain, remove_his_tag=False):
    '''
    Gets the sequence of a chain in a structure
    Arguments: 
        struct (Biopandas.pdb)
        chain (str): identity of the chain of interest
        remove_his_tag (bool): remove the 6x his tage if present {default: False}
    Returns: 
        str: sequence of the chain with 1-letter AA code from struct
    '''
    seq = struct.amino3to1()
    sequence = ''.join(seq.loc[seq['chain_id'] == chain, 'residue_name'])

    # if requested remove hisx6 tag
    if remove_his_tag:
        his_tag = "H"*6
        if sequence[:6] == his_tag:
            sequence = sequence[6:]
        if sequence[-6:] == his_tag:
            sequence = sequence[:-6]

    return sequence


struct1 = 'Renumbered_aligned_pdbs_allbackbone_partial/1aq7.pdb'


def get_RMSD(struct1, struct2, chains, atomtypes='all'):
    '''
    Calculates the RMSD between two structures using all heavy atoms or CA only. Only considers highest occupancy conformers.
    Arguments:
        struct1, struct2 (str): Paths to .pdb files. May also be a pre-built coord dict.
        chains = list: List of renamed chains to consider for RMSD, e.g. ['A']
        atomtypes = list: Atom types to calculate RMSD on. 'all' or list, e.g. ['CA','CB']
    Returns:
        RMSD (float)
    '''
    parser = Bio.PDB.PDBParser(QUIET=True)  # Initialize parser
    if type(struct1) == str and type(struct2) == str:
        structure1 = parser.get_structure('1', struct1)[0]  # Read .pdb file
        structure2 = parser.get_structure('2', struct2)[0]  # Read .pdb file
        coords = {}
        for struct in [structure1, structure2]:
            for chain_num in range(0, len(chains)):
                for res in sorted([chain for chain in struct.get_chains()], key=lambda x: x.get_id())[chain_num].get_unpacked_list()[:
                                                                                                                                     len([x for x in struct.get_residues() if x.id[0] == ' '])]:
                    if res.get_resname() not in ['ALA', 'ARG', 'ASN', 'ASP', 'CYS',
                                                 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']:
                        continue
                    for atom in res.get_unpacked_list():
                        # Append all coordinates into coord dictionary
                        if ((atom.get_name() not in atomtypes) and (atomtypes != 'all')) or (atom.element == 'H' or atom.element == 'D'):
                            continue
                        if atom.occupancy != 1:
                            occupancies = [
                                x.occupancy for x in res.get_unpacked_list() if x.name == atom.name]
                            occupancies.sort(reverse=True)
                            if atom.occupancy != occupancies[0]:
                                continue
                        if (atom.get_parent().get_parent().get_id()+'_' +
                            str(atom.get_parent().get_id()[1])+'_' +
                            str(atom.get_parent().get_id()[2])+'_' +
                                atom.get_name()) in coords.keys():
                            coords[(atom.get_parent().get_parent().get_id()+'_' +
                                    str(atom.get_parent().get_id()[1])+'_' +
                                    str(atom.get_parent().get_id()[2])+'_' +
                                    atom.get_name())].append(atom.coord)
                        else:
                            coords[(atom.get_parent().get_parent().get_id()+'_' +
                                    str(atom.get_parent().get_id()[1])+'_' +
                                    str(atom.get_parent().get_id()[2])+'_' +
                                    atom.get_name())] = [atom.coord]
    if type(struct1) == str and type(struct2) == dict:
        structure1 = parser.get_structure('1', struct1)[0]  # Read .pdb file
        coords = struct2.copy()
        for struct in [structure1]:
            for chain_num in range(0, len(chains)):
                for res in sorted([chain for chain in struct.get_chains()], key=get_chain_name)[chain_num].get_unpacked_list()[:
                                                                                                                               len([x for x in struct.get_residues() if x.id[0] == ' '])]:
                    if res.get_resname() not in ['ALA', 'ARG', 'ASN', 'ASP', 'CYS',
                                                 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']:
                        continue
                    for atom in res.get_unpacked_list():
                        # Append all coordinates into coord dictionary
                        if ((atom.get_name() not in atomtypes) and (atomtypes != 'all')) or (atom.element == 'H' or atom.element == 'D'):
                            continue
                        if atom.occupancy != 1:
                            occupancies = [
                                x.occupancy for x in res.get_unpacked_list() if x.name == atom.name]
                            occupancies.sort(reverse=True)
                            if atom.occupancy != occupancies[0]:
                                continue
                        if (atom.get_parent().get_parent().get_id()+'_' +
                            str(atom.get_parent().get_id()[1])+'_' +
                            str(atom.get_parent().get_id()[2])+'_' +
                                atom.get_name()) in coords.keys():
                            coords[(atom.get_parent().get_parent().get_id()+'_' +
                                    str(atom.get_parent().get_id()[1])+'_' +
                                    str(atom.get_parent().get_id()[2])+'_' +
                                    atom.get_name())].append(atom.coord)
                        else:
                            coords[(atom.get_parent().get_parent().get_id()+'_' +
                                    str(atom.get_parent().get_id()[1])+'_' +
                                    str(atom.get_parent().get_id()[2])+'_' +
                                    atom.get_name())] = [atom.coord]
    if type(struct1) == dict and type(struct2) == str:
        structure2 = parser.get_structure('2', struct1)[0]  # Read .pdb file
        coords = struct1.copy()
        for struct in [structure2]:
            for chain_num in range(0, len(chains)):
                for res in sorted([chain for chain in struct.get_chains()], key=get_chain_name)[chain_num].get_unpacked_list()[:
                                                                                                                               len([x for x in struct.get_residues() if x.id[0] == ' '])]:
                    if res.get_resname() not in ['ALA', 'ARG', 'ASN', 'ASP', 'CYS',
                                                 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']:
                        continue
                    for atom in res.get_unpacked_list():
                        # Append all coordinates into coord dictionary
                        if ((atom.get_name() not in atomtypes) and (atomtypes != 'all')) or (atom.element == 'H' or atom.element == 'D'):
                            continue
                        if atom.occupancy != 1:
                            occupancies = [
                                x.occupancy for x in res.get_unpacked_list() if x.name == atom.name]
                            occupancies.sort(reverse=True)
                            if atom.occupancy != occupancies[0]:
                                continue
                        if (atom.get_parent().get_parent().get_id()+'_' +
                            str(atom.get_parent().get_id()[1])+'_' +
                            str(atom.get_parent().get_id()[2])+'_' +
                                atom.get_name()) in coords.keys():
                            coords[(atom.get_parent().get_parent().get_id()+'_' +
                                    str(atom.get_parent().get_id()[1])+'_' +
                                    str(atom.get_parent().get_id()[2])+'_' +
                                    atom.get_name())].append(atom.coord)
                        else:
                            coords[(atom.get_parent().get_parent().get_id()+'_' +
                                    str(atom.get_parent().get_id()[1])+'_' +
                                    str(atom.get_parent().get_id()[2])+'_' +
                                    atom.get_name())] = [atom.coord]
    if type(struct1) == dict and type(struct2) == dict:
        coords = struct1.copy()
        for entry in struct2.items():
            try:
                coords[entry[0]].append(entry[1])
            except:
                coords[entry[0]] = entry[1]

    # Obtain atoms shared by both structures

    coords_shared = {}
    for entry in coords.items():
        if len(entry[1]) == 2:
            coords_shared[entry[0]] = entry[1]

    deviations = dict.fromkeys(coords_shared.keys(), [])
    for atom in deviations:
        deviations[atom] = np.linalg.norm(
            coords_shared[atom][1]-coords_shared[atom][0])
    if len(deviations.keys()) == 0:
        print('No RMSD could be calculated for structures: {},{} and chains {}'.format(
            struct1, struct2, chains))
        return 0
    rmsd = np.sqrt(sum([(deviation**2) for deviation in [float(entry[1]) for entry in deviations.items()
                                                         ]])/len(deviations.keys()))

    return rmsd


def relabel_struct(struct, chain_map, residue_map, top_refs):
    '''
    Arguments:
        struct (Biopandas.pdb)
        chain_map (dict): old-to-new dictionary of chain identities
        residue_map (dict): old-to-new dictionary of residues in "chain | residue_number | isnert" format
    '''
    # TODO none-should be nans otherwise .fillna(df['col1'])

    # move top refs to ABC etc, moving the chains that were neamed that
    # print(f"before {chain_map}, to change {top_refs}")
    chain_index_ref = 65  # chr(65)==A
    for top_ref in top_refs:
        old_reassignment = chain_map[top_ref]
        for key, value in chain_map.items():
            if value == chr(chain_index_ref):
                chain_map[key] = old_reassignment
        chain_map[top_ref] = chr(chain_index_ref)
        chain_index_ref += 1

    # print(f"after {chain_map}")
    # print(struct.df['ATOM']['chain_id'].unique())
    # print(chain_map)
    struct.df['ATOM']['chain_id'] = struct.df['ATOM']['chain_id'].map(
        chain_map).fillna("Z")
    # print(struct.df['ATOM']['chain_id'].unique())
    struct.df['HETATM']['chain_id'] = struct.df['HETATM']['chain_id'].map(
        chain_map).fillna("Z")
    struct.df['ANISOU']['chain_id'] = struct.df['ANISOU']['chain_id'].map(
        chain_map).fillna("Z")

    # TODO residue_map for het and ani??
    struct.df["ATOM"]["ID"] = struct.df["ATOM"][["chain_id", "residue_number",
                                                 "insertion"]].apply(lambda x: '|'.join(x.astype(str)), axis=1)
    struct.df["ATOM"]["ID"] = struct.df["ATOM"]["ID"].map(
        residue_map).fillna(struct.df["ATOM"]["ID"])
    struct.df["ATOM"][["chain_id", "residue_number", "insertion"]] = struct.df[
        "ATOM"]["ID"].apply(lambda x: pd.Series(x.split("|")))
    struct.df["ATOM"] = struct.df["ATOM"].drop("ID", axis=1)
    return struct

# COMMENT: wrong return value??


def tranform_coords_pandaspdb(struct, transform):
    '''
    Arguments:
        struct (Biopandas.pdb)
        transform np.array matrix
    '''
    # print(struct)
    coords = biopandas_pdb_to_matrix(struct)
    new_coords = transform_coords_matrix(transform, coords)
    new_struct = replace_biopandas_coords(struct, new_coords)
    return struct


def transform_coords_matrix(T, X):
    # print(X.shape)
    trans_1 = T[0]
    rotation = T[1]
    trans_2 = T[2]
    # print(trans_1.shape)
    X = X + trans_1
    X = X @ rotation
    X = X + trans_2
    return X


def biopandas_pdb_to_matrix(struct):
    # print(struct)
    atoms = struct.df['ATOM']
    atom_coord = atoms[['x_coord', 'y_coord', 'z_coord']].to_numpy()
    return np.array(atom_coord)


def replace_biopandas_coords(struct, coords):
    struct.df['ATOM'][['x_coord', 'y_coord', 'z_coord']] = coords
    return struct

# COMBINE FUNCITON INPUT MULTIPLE CSV multiple folder multiple names for new col


def fetch_and_save_pdb(pdbid, directory):
    '''
    Download the given pdb into directory, 
    if already there do not over-write.
    Arguments:
        pdbid (str): 4 character pdbid "XXXX"
        directory (str): directory to save pdbid to
    Output:
        Saves pdb of the pdbid to the directory
    '''
    pdb_file = Path(directory, f"{pdbid.lower()}.pdb")
    if not Path(pdb_file).is_file():
        struct = get_pdb_struct(pdbid)
        if struct is not None:
            struct.to_pdb(path=pdb_file, gz=False)


def get_chain_position_list(struct, chain):
    '''
    Given a BiopandasPDB.df["ATOM"] and a chain, return a list
    of the 1 character residue names and a list of their positions
    in "chain | residue_number | insert" format. 
    Arguments:
        struct (BiopandasPDB.df["ATOM"])
        chain (str): chain to get list
    Return:
        return a list of the 1-letter residue names and a list of their positions
        in "chain | residue_number | insert" format. 
    '''
    struct['ID'] = struct[["chain_id", "residue_number", "insertion"]].apply(
        lambda x: '|'.join(x.astype(str)), axis=1)
    struct_chain = struct[struct['chain_id'] == chain]
    ref_list = struct_chain[["residue_name", "ID"]].groupby(
        ["ID"], sort=False).first().reset_index()
    ref_list["residue_name"] = ref_list["residue_name"].apply(
        str.title).map(SeqUtils.IUPACData.protein_letters_3to1)
    return ref_list["residue_name"].tolist(), ref_list["ID"].tolist()


def get_chain_position_list_from_seq(seq):
    ''' this is a function for renumber
    if desire to renumber with just a sequence, it will name it all
    chain A and just sequential numbering from there.
    Arguments:
        seq (str): str of 1-letter AA sequences
    Returns:
        the chain list in "chain | residue_number | isnert" format
        using jsut simple sequential numbering
    '''
    poss = []
    for i in range(len(seq)):
        poss.append(f"A|{i+1}|")
    return list(seq), poss


def retrieve_pdbs(PDB_list, directory):
    '''
    Given a table with pdb ids in 'Entry ID', download
    all the pdb to the directory
    '''
    tqdm.pandas(desc='Retrieve PDBs')
    PDB_list['Entry ID'].progress_apply(
        lambda pdb: fetch_and_save_pdb(pdb, directory))
