'''EnsemblePDB.analysis.atoms

Functions to analyze deviation of atom positions within the aligned ensembles.

Authors:
    Rachael Kretsch(rkretsch@stanford.edu), 
    Siyuan Du (dusiyuan@stanford.edu),
    Jacob Parres-Gold
'''

from pathlib import Path
import pandas as pd
from tqdm import tqdm

from EnsemblePDB.utils.file_management import get_nonexistant_file
from EnsemblePDB.utils.biopandas_utils import get_PandasPDBs, get_xyz, clean_multiconformers
from EnsemblePDB.utils.analysis_utils import *


def get_MDev(directory, chains, reference_PDB, multiconformers=False,
             output_directory=None, bootstrap=True, iter=50):
    '''
    Takes .pdb files and calculates the mean deviation of atom 
    positions for all atoms. Note no alignment is perferformed.

    Arguments:
        directory (str): Folder of aligned .pdb files
        chains (list of str): List of (renamed) chains to analyze
        reference_PDB (str): Reference PDB to obtain wild-type sequence
        mutliconformers (bool): if True, count multiconformers as individual 
                        structures (with equal weight); if False, delete all 
                        multiconformer atoms {default: False}
        output_directory (str): Directory to save output MDev files to. If 
                        None, saves to the parent of given directory 
                        {default: None}
        bootstrap (bool): if bootstrap is True, perform bootstrap analysis and 
                        include standard deviation of the bootstrap 
                        distribution. Notice that this may take a long time to 
                        complete for large number of atoms {default: True}
        iter (int): number of iteration to resample data for bootstrap analysis 
                        {default: 50}

    Returns:
        DataFrames of MDev values for each atom.

    Output:
        Saves MDev files to given output directory or default.
    '''
    ppdbs = get_PandasPDBs(directory)
    atoms = get_xyz(ppdbs)
    if not multiconformers:
        atoms = clean_multiconformers(atoms)

    if output_directory is None:
        output_directory = str(Path(directory).parents[0])
    MDev = combine_save_MDev(atoms, chains, reference_PDB,
                             output_directory, bootstrap, iter)
    return MDev


def get_distance_distributions(directory, chains, multiconformers=False,
                               quantile=0.95, report_outliers=True,
                               output_directory=None):
    '''
    Get all coordinates for an atom in the ensemble and calculate the distance 
    of each from the center atom of the ensemble and from the average position.
    Also determines whether a position is an outlier of the ensemble based on 
    given quantile. 

    Arguments:
        directory (str): folder of all aligned PDB file.
        chains (list): list of selected chain to get atoms 
        multiconformers (bool): if False will remove all multiconfomers
                                {defaults: False}
        quantile (float): the quantile of data used to identify outliers 
                            {default: 0.95} 
        report_outliers (bool): if True, calculates the number of outliers in 
                            each PDB structure and saves a separate file. 
                            {default: True}
        output_directory (str): Directory to save output files to. If None, 
                        saves to the parent of given directory {default: None}

    Returns:
        dataframe with new columns of distance from center atom and average position 
        and outlier

    Output:
        Saves atomic positions and their distances from center/average position
        as a csv. Saves outlier report if user chose to get it.
    '''
    ppdbs = get_PandasPDBs(directory)
    atoms = get_xyz(ppdbs)
    if not multiconformers:
        atoms = clean_multiconformers(atoms)

    if output_directory is None:
        output_directory = str(Path(directory).parents[0])
    final = {}

    for chain in chains:
        atoms = atoms.loc[atoms['chain_id'] == chain]
        atoms = atoms.set_index(
            ['residue_number', 'insertion', 'residue_name',
             'atom_name']).sort_index()
        to_concat = []

        for i, index in tqdm(enumerate(atoms.index.unique()),
                             total=len(atoms.index.unique()),
                             desc=f'Calculating distances and' +
                             'outliers in chain {chain}'):
            atom = atoms.loc[index, :].copy()
            atom['xyz_spread'] = (calculate_msd(atom['x_coord']) +
                                  calculate_msd(atom['y_coord']) +
                                  calculate_msd(atom['z_coord']))

            # find center atom
            center = atom.iloc[atom['xyz_spread'].argmin()]
            # find average position
            average = (atom['x_coord'].mean(),
                       atom['y_coord'].mean(), atom['z_coord'].mean())

            atom['distance_from_center'] = atom.apply(lambda x: get_distance(
                x['x_coord'], x['y_coord'], x['z_coord'], center['x_coord'],
                center['y_coord'], center['z_coord']), axis=1)
            atom['outlier_from_center'] = atom['distance_from_center'].apply(
                lambda x: x > atom['distance_from_center'].quantile(quantile))

            atom['distance_from_average'] = atom.apply(lambda x: get_distance(
                x['x_coord'], x['y_coord'], x['z_coord'], average[0],
                average[1], average[2]), axis=1)
            atom['outlier_from_average'] = atom['distance_from_average'].apply(
                lambda x: x > atom['distance_from_average'].quantile(quantile))
            to_concat.append(atom)
        chain_out = pd.concat(to_concat)
        fname = get_nonexistant_file(
            f'{output_directory}/all_atom_positions_summary_chain_{chain}.csv')
        chain_out.to_csv(fname)
        print(f'\nAtomic positions and distances saved to {fname}')
        final[chain] = chain_out
    if report_outliers:
        number_of_outliers = get_outliers_report(final, output_directory)

    return final
