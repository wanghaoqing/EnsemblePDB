'''EnsemblePDB.analyze.rotamer

Calculates all psi, phi and chi angles of PDB files.

Authors: 
    Siyuan Du (dusiyuan@stanford.edu)

'''
import Bio.PDB
import numpy as np
import pandas as pd
import math
from pathlib import Path
from tqdm import tqdm
from glob import glob

from EnsemblePDB.utils.file_management import get_nonexistant_file

def get_all_rotamers(directory, chains,output_directory=None):
    if output_directory is None:
        output_directory = str(Path(directory).parents[0])
    all_rotamers = {}
    for chain in chains:
        rotamer = []
        all_pdbs = glob(directory+'/*.pdb')
        for pdb in tqdm(all_pdbs, total = len(all_pdbs), desc=f'Calculating rotamer angle for chain {chain}'):
            df = get_dihedrals(pdb=pdb, chain_id=chain)
            df['Entry ID'] = pdb.split('/')[-1][:-4]
            rotamer.append(df)
        all_df = pd.concat(rotamer)
        all_rotamers[chain] = all_df
        fname = get_nonexistant_file(f'{output_directory}/all_residue_rotamers_{chain}.csv')
        all_df.to_csv(fname)
        print(f'\nRotamer angle distributions saved to {fname}')
    return all_rotamers

# definition of chi angles
CHI1 = {
    'ARG': ['N','CA','CB', 'CG'],
    'ASN': ['N','CA','CB', 'CG'],
    'ASP': ['N','CA','CB', 'CG'],
    'CYS': ['N','CA','CB', 'SG'],
    'GLU': ['N','CA','CB', 'CG'],
    'GLN': ['N','CA','CB', 'CG'],
    'HIS': ['N','CA','CB', 'CG'],
    'ILE': ['N','CA','CB', 'CG1'],
    'LEU': ['N','CA','CB', 'CG'],
    'LYS': ['N','CA','CB', 'CG'],
    'MET': ['N','CA','CB', 'CG'],
    'PHE': ['N','CA','CB', 'CG'],
    'PRO': ['N','CA','CB', 'CG'],
    'SER': ['N','CA','CB', 'OG'],
    'THR': ['N','CA','CB', 'OG1'],
    'TRP':['N','CA','CB', 'CG'],
    'TYR': ['N','CA','CB', 'CG'],
    'VAL': ['N','CA','CB', 'CG1']
}
CHI2 = {
    'ARG': ['CA','CB','CG','CD'],
    'GLU': ['CA','CB','CG','CD'],
    'GLN': ['CA','CB','CG','CD'],
    'LYS': ['CA','CB','CG','CD'],
    'MET': ['CA','CB','CG','SD']
    # 'PHE': ['CA','CB','CG','CD1'], #nonrotameric
    # 'TRP': ['CA','CB','CG','CD1'], #nonrotameric
    # 'TYR': ['CA','CB','CG','CD1'] #nonrotameric
    }
CHI3 = {
    'ARG': ['CB','CG','CD','NE'],
    'LYS': ['CB','CG','CD','CE'],
    'MET': ['CB','CG','SD','CE'],
}
CHI4 = {
    'ARG': ['CG','CD','NE','CZ'],
    'LYS': ['CG','CD','CE','NZ'],   
}

def get_dihedrals(pdb, chain_id):
    '''
    For a given pdb file and a protein chain, calculate all possible psi,phi and chi angles.
    Note: for any dihedral calculation, if any of the four atom is disorderd (multiconfomer), this angle will be skipped (np.nan).
    '''
    parser = Bio.PDB.PDBParser(QUIET=True)
    # Entry_ID = pdb.split('/')[-1].split('.')[0]
    structure = parser.get_structure("_", pdb)
    chain = structure[0][chain_id]
    all_residue_angles = []
    for residue in chain.get_list():
        resid = residue.get_id()
        # skip if HETATM flag is true
        if resid[0] != " ":
            continue
        dihedral_angles = {}
        dihedral_angles['residue_number'] = resid[1]
        dihedral_angles['residue_name'] = residue.resname
        dihedral_angles['insertion'] = resid[2]
        dihedral_angles['phi'] = calcphi(chain, resid)
        dihedral_angles['psi'] = calcpsi(chain, resid)
        # calculate chi angles
        for chi in ['chi1', 'chi2','chi3','chi4']:
            dihedral_angles[chi] = np.nan
        for chi_atoms, chi_name in zip([CHI1, CHI2, CHI3, CHI4],['chi1', 'chi2','chi3','chi4']):
            if residue.resname not in chi_atoms.keys():
                continue
            atoms_to_find = chi_atoms[residue.resname]
            atoms_in_pdb = [a.name for a in residue.get_list()]
            # check if atom is in PDB & if they are disordered
            if any(list(map(lambda a: a not in atoms_in_pdb, atoms_to_find))):
                continue
            if any(list(map(lambda a: residue[a].is_disordered == 0, atoms_to_find))):
                continue
            atomobj = [residue[a] for a in atoms_to_find]
            dihedral_angles[chi_name] = calcdihedral(*atomobj)
        all_residue_angles.append(dihedral_angles)
    df = pd.DataFrame.from_records(all_residue_angles)
    return df

# def get_residue_by_seqproximity(chain, reference_resid, position):
    # '''
    # position: +1 = next, -1 = previous, etc.
    # '''
    # resid_list = []
    # for res in chain:
    #     if res.get_id()[0] != " ":
    #         continue
    #     resid_list.append(res.get_id())
    # # print(len(resid_list))
    # try:
    #     ref_index = resid_list.index(reference_resid)
    # except:
    #     print(f'{reference_resid} not found')
    #     return
    # if ref_index+position < len(resid_list) and ref_index+position >= 0:
    #     return resid_list[ref_index+position]
    # else:
    #     return
def get_residue_by_seqproximity(chain, reference_resid, position):
    '''
    position: +1 or -1 : next residue or previous residue
    '''
    atoms  = Bio.PDB.Selection.unfold_entities(chain, 'A')
    ns = Bio.PDB.NeighborSearch(atoms)
    curr_res = chain[reference_resid]
    if position == 1:
        if 'C' in [x.name for x in curr_res]:
            target_atom = curr_res['C']
        else: return
    if position == -1:
        if 'N' in [x.name for x in curr_res]:
            target_atom = curr_res['N']
        else: return
    cov_atoms = ns.search(target_atom.coord, 2)
    candidates = []
    for a in cov_atoms:
        if a.get_parent().id != reference_resid:
            candidates.append(a)
    neighbor_atom = None
    if len(candidates) == 0:
        return 
    elif len(candidates) == 1:
        neighbor_atom = candidates[0]
    else:
        distances = [np.linalg.norm(x.coord-target_atom.coord) for x in candidates]
        neighbor_atom = candidates[np.array(distances).argmin()]
    # check atom identity
    if position == 1:
        if neighbor_atom.name != 'N':
            return 
    if position == -1:
        if neighbor_atom.name != 'C':
            return
    return neighbor_atom.get_parent().id

def calcphi(chain,res_id):
    res_id__1 = get_residue_by_seqproximity(chain=chain, reference_resid= res_id, position=-1)
    if not res_id__1:
        return np.nan
    try:
        CO_atom = chain[res_id__1]['C']
        N_atom = chain[res_id]['N']
        CA_atom = chain[res_id]['CA']
        C_atom = chain[res_id]['C']
    except KeyError:
        return np.nan
    #check if disordered
    if CO_atom.is_disordered() != 0 or N_atom.is_disordered() != 0 or CA_atom.is_disordered() != 0 or C_atom.is_disordered() !=0:
        # print('DEBUG: disordered atom')
        return np.nan
    C0 = Bio.PDB.vectors.Vector(CO_atom.coord)
    N = Bio.PDB.vectors.Vector(N_atom.coord)
    CA = Bio.PDB.vectors.Vector(CA_atom.coord)
    C = Bio.PDB.vectors.Vector(C_atom.coord)
    phi = Bio.PDB.vectors.calc_dihedral(C0,N,CA,C)
    return phi * 180/math.pi

def calcpsi(chain,res_id):
    res_id_1 = get_residue_by_seqproximity(chain=chain, reference_resid=res_id, position=1)
    if not res_id_1:
        return np.nan
    try:
        N_atom = chain[res_id]['N']
        CA_atom = chain[res_id]['CA']
        C_atom = chain[res_id]['C']
        N1_atom = chain[res_id_1]['N']
    except KeyError:
        return np.nan
    if N_atom.is_disordered() != 0 or CA_atom.is_disordered() != 0 or C_atom.is_disordered() != 0 or N1_atom.is_disordered() !=0 :
        # print('DEBUG: disordered atom')
        return np.nan
    N = Bio.PDB.vectors.Vector(chain[res_id]['N'].coord)
    CA = Bio.PDB.vectors.Vector(chain[res_id]['CA'].coord)
    C = Bio.PDB.vectors.Vector(chain[res_id]['C'].coord)
    N1 = Bio.PDB.vectors.Vector(chain[res_id_1]['N'].coord)
    psi = Bio.PDB.vectors.calc_dihedral(N,CA,C,N1)
    return psi * 180/math.pi

def calcchi1(chain, res_id,chi1s = CHI1):
    atoms = CHI1[chain[res_id].resname]
    atom_coords = []
    for a in atoms:
        try:
            a_coord = Bio.PDB.vectors.Vector(chain[res_id][a].coord)
        except KeyError:
            return np.nan
        atom_coords.append(a_coord)
    chi1 = Bio.PDB.vectors.calc_dihedral(*atom_coords)
    # XG = 'init'
    # for atom in chain[res_id]:
    #     if len(atom.name) > 1 and atom.name[1] == 'G':
    #         XG = Bio.PDB.vectors.Vector(chain[res_id][atom.name].coord)
    # if XG == 'init':
    #     return np.nan
    # chi1 = Bio.PDB.vectors.calc_dihedral(N,CA,CB,XG)
    return chi1 * 180/math.pi

def calcdihedral(a1,a2,a3,a4):
    v1 = Bio.PDB.vectors.Vector(a1.coord)
    v2 = Bio.PDB.vectors.Vector(a2.coord)
    v3 = Bio.PDB.vectors.Vector(a3.coord)
    v4 = Bio.PDB.vectors.Vector(a4.coord)
    dihedral = Bio.PDB.vectors.calc_dihedral(v1,v2,v3,v4)
    return dihedral * 180/math.pi