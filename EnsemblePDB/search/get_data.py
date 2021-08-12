""" EnsemblePDB.search.get_data

TODO 

Author:
    Rachael Kretsch (rkretsch@stanford.edu), 
    Siyuan Du, 
    Jacob Parres-Gold
"""

from pathlib import Path
from tqdm import tqdm

import pandas as pd

from EnsemblePDB.utils.file_management import get_nonexistant_file
from EnsemblePDB.utils.table_utils import get_list_from_table_entry, delete_all_one_value
from EnsemblePDB.utils.search_utils import fetch_pdb_data, combine_dict_info


###############################################################################
# Global data
###############################################################################

# data can include any of these: https://data.rcsb.org/data-attributes.html

minimum_info = {
    'citation.pdbx_database_id_doi': 'citation: doi',
    'citation.title': 'citation: title',
    'citation.year': 'citation: year',
    'diffrn.ambient_temp': 'expt: diffrn temp',
    'exptl_crystal_grow.method': 'expt: crystal grow method',
    'exptl_crystal_grow.ph': 'expt: crystal grow pH',
    'pdbx_vrpt_summary.pdbr': 'pdbx: Rfactor',
    'pdbx_vrpt_summary.pdbrfree': 'pdbx: Rfree',
    'pdbx_vrpt_summary.pdbresolution': 'pdbx: resolution',
    'pdbx_vrpt_summary.clashscore': 'pdbx: clashscore',
    'rcsb_entry_info.molecular_weight': 'rcsb-entry: molecular weight',
    'rcsb_entry_info.polymer_molecular_weight_maximum': 'rcsb-entry: polymer molecular weight max',
    'rcsb_entry_info.polymer_molecular_weight_minimum': 'rcsb-entry: polymer molecular weight min',
    'rcsb_entry_info.polymer_entity_count': 'rcsb-entry: polymer entity count',
    'rcsb_entry_info.polymer_entity_count_protein': 'rcsb-entry: polymer entity protein count',
    'rcsb_entry_info.nonpolymer_entity_count': 'rcsb-entry: nonpolymer entity count',
    'refine.biso_mean': 'refine_bfactor_mean',
    'rcsb_polymer_entity.pdbx_description': 'rscb-polymer-entity: description',
    'rcsb_polymer_entity.pdbx_ec': 'rscb-polymer-entity: EC number',
    'rcsb_polymer_entity.pdbx_mutation': 'rscb-polymer-entity: mutation',
    'rcsb_entry_container_identifiers.entity_ids': 'rscb-entry: entity ids',
    'pdbx_entity_nonpoly.comp_id': 'pdbx-entity-nonpolymer: chem component id',
    'pdbx_entity_nonpoly.entity_id': 'pdbx-entity-nonpolymer: entity id',
    'pdbx_entity_nonpoly.name': 'pdbx-entity-nonpolymer: name',
    'rcsb_nonpolymer_entity.formula_weight': 'rcsb-nonpolymer-entity: formula weight',
    'rcsb_entity_source_organism.scientific_name': 'pdbx-entity: organism scientific name',
    'entity_poly.rcsb_mutation_count': 'entity-polymer: mutation_count'}

other_data = {'cell.angle_alpha': 'unit-cell-angle_alpha (deg)',
              'cell.angle_beta': 'unit-cell-angle_beta (deg)',
              'cell.angle_gamma': 'unit-cell-angle_gamma (deg)',
              'cell.length_a': 'unit-cell-length_a (A)',
              'cell.length_b': 'unit-cell-length_b (A)',
              'cell.length_c': 'unit-cell-length_c (A)',
              'cell.zpdb': 'unit-cell_number_of_chains_z',
              'exptl_crystal_grow.pdbx_details': 'expt-crystal-grow_details',
              'pdbx_vrpt_summary.bfactor_type': 'pdbx-bfactor_type',
              'pdbx_vrpt_summary.dccr': 'pdbx-dcc_R',
              'pdbx_vrpt_summary.dccrfree': 'pdbx-dcc_Rfree',
              'pdbx_vrpt_summary.dccrefinement_program': 'pdbx-dcc_refine_program',
              'pdbx_vrpt_summary.edsr': 'pdbx-eds_Rfactor',
              'pdbx_vrpt_summary.edsresolution': 'pdbx-eds_resolution_highlimit',
              'pdbx_vrpt_summary.edsresolution_low': 'pdbx-eds_resolution_lowlimit',
              'pdbx_vrpt_summary.fo_fc_correlation': 'pdbx-fo-fc_correlation',
              'pdbx_vrpt_summary.pdbresolution_low': 'pdbx-PDB_resolution_lowlimit',
              'pdbx_vrpt_summary.angles_rmsz': 'pdbx-Molprobity-RMSZ_angles',
              'pdbx_vrpt_summary.bonds_rmsz': 'pdbx-Molprobity-RMSZ_bonds',
              'rcsb_id': 'Uniprot_id',
              'rcsb_cluster_membership.cluster_id': 'cluster_id',
              'rcsb_cluster_membership.identity': 'cluster_identity',
              'rcsb_entry_info.deposited_modeled_polymer_monomer_count': 'rcsb-entry_number_monopolymers',
              'rcsb_entry_info.deposited_nonpolymer_entity_instance_count': 'rcsb-entry_number_nonpolymer_entity',
              'rcsb_entry_info.inter_mol_covalent_bond_count': 'rcsb-entry_intermol_covalent_bond_count',
              'rcsb_entry_info.inter_mol_metalic_bond_count': 'rcsb-entry_intermol_metal_bond_count',
              'rcsb_entry_info.solvent_entity_count': 'rcsb-entry_solvent_entity_count',
              'rcsb_uniprot_feature.description': 'feature_descriptions',
              'rcsb_uniprot_feature.feature_positions.beg_seq_id': 'feature_beg_seq_ids',
              'rcsb_uniprot_feature.feature_positions.end_seq_id': 'feature_end_seq_ids',
              'rcsb_uniprot_feature.type': 'feature_types',
              'refine-details': 'refine_details'}

other_data.update(minimum_info)

cols_to_group = [['rcsb_uniprot_feature.feature_positions.beg_seq_id', 'rcsb_uniprot_feature.type',
                  'rcsb_uniprot_feature.feature_positions.end_seq_id', 'rcsb_uniprot_feature.description']]


def get_pdb_data(ensemble_csv, desired_cols=minimum_info, label='MyProtein', info_list=["entry", "entity", "chem"], cols_to_group=[], output_dir=None):
    '''
    Given a csv with PDB IDs in the 'Entry ID' column, downloads
    and merges summary reports for those PDB ids from the PDB.

    Arguments: 
        ensemble_csv (str): location of csv file with the PDB IDs
        desired_cols (dict): keys are the name of the PDB information, values are
            what you would like this column to be names
        label (str): Label to add to output file name
        info_list (list of str): list of what type of attributes you want
            note "entry" must be one of them
        cols_to_group (list of lists): Used if there are columns expected to have more than one entry, such as sequence annotations, where each
            site is expected to have a set of data: "type", "description," etc. These data will need to be paired per site, potentially 
            inserting "None" (str) if no data is available for a site to keep everything aligned. Each list contains strings of the column 
            names that will be included in that group. 
            Example: cols_to_group = [['rcsb_uniprot_feature.type','rcsb_uniprot_feature.description'], ...]

    Returns: 
        dataframe of the summary_report
        the saved dataframe (summary_report*.csv)
    '''
    # get pdb ids from csvs
    data = pd.read_csv(ensemble_csv)
    assert ('Entry ID' in data.columns), 'the inputted csv must contain a Entry ID column with pdb ids'

    # Make sure we add the columns they need to get data outside of entry data.
    if "entry" not in info_list:
        print("ERROR: entry needs to be included to get rest of data")
        return
    if "entity" in info_list or "entity: uniprot" in info_list:
        if 'rcsb_entry_container_identifiers.entity_ids' not in desired_cols:
            print("'rcsb_entry_container_identifiers.entity_ids' needed to obtain entity information adding to desired columns")
            desired_cols['rcsb_entry_container_identifiers.entity_ids'] = 'rcsb_entry_container_identifiers.entity_ids'
    if "instance" in info_list:
        if 'rcsb_entry_info.deposited_nonpolymer_entity_instance_count' not in desired_cols:
            print("'rcsb_entry_info.deposited_nonpolymer_entity_instance_count' needed to obtain entity information adding to desired columns")
            desired_cols['rcsb_entry_info.deposited_nonpolymer_entity_instance_count'] = 'rcsb_entry_info.deposited_nonpolymer_entity_instance_count'
        if 'rcsb_entry_info.deposited_polymer_entity_instance_count' not in desired_cols:
            print("'rcsb_entry_info.deposited_polymer_entity_instance_count' needed to obtain entity information adding to desired columns")
            desired_cols['rcsb_entry_info.deposited_polymer_entity_instance_count'] = 'rcsb_entry_info.deposited_polymer_entity_instance_count'
    if "assembly" in info_list:
        if 'rcsb_entry_container_identifiers.assembly_ids' not in desired_cols:
            print("'rcsb_entry_container_identifiers.assembly_ids' needed to obtain entity information adding to desired columns")
            desired_cols['rcsb_entry_container_identifiers.assembly_ids'] = 'rcsb_entry_container_identifiers.assembly_ids'
    if "chem" in info_list or "chem: drugbank" in info_list:
        assert("entity" in info_list), "Need to collect chem names from entity. Please include entity in info_list."
        if 'pdbx_entity_nonpoly.comp_id' not in desired_cols:
            print(
                "'pdbx_entity_nonpoly.comp_id' needed to obtain entity information adding to desired columns")
        if 'rcsb_entry_info.nonpolymer_entity_count' not in desired_cols:
            print("'rcsb_entry_info.nonpolymer_entity_count' needed to determine whether chem data exists, add to desired columns")
            desired_cols['pdbx_entity_nonpoly.comp_id'] = 'pdbx_entity_nonpoly.comp_id'

    info_not_needed = ["entry", "entry: pubmed", "entity",
                       "entity: uniprot", "instance", "assembly", "chem", "chem: drugbank"]
    for info_desired in info_list:
        if info_desired in info_not_needed:
            info_not_needed.remove(info_desired)
        else:
            print("ERROR: your info_list contains an invalid option.")
            return

    # get pdb_id list
    pdb_ids = data['Entry ID'].values.tolist()
    all_info = {}
    for pdb in tqdm(pdb_ids, total=len(pdb_ids), desc='Adding data'):

        if "entry" not in info_not_needed:
            info = fetch_pdb_data(
                'pdb_entry', pdb, desired_cols, cols_to_group)
        if "entry: pubmed" not in info_not_needed:
            pubmed_info = fetch_pdb_data('pubmed_entry',
                                         pdb, desired_cols, cols_to_group)
            info = combine_dict_info(info, pubmed_info, [])

        if "entity" not in info_not_needed or "entity: uniprot" not in info_not_needed:
            old_info = list(info.keys())
            entities = get_list_from_table_entry(
                info['rcsb_entry_container_identifiers.entity_ids'])
            for entity in entities:
                if "entity" not in info_not_needed:
                    entity_info = fetch_pdb_data('pdb_entity',
                                                 f'{pdb}/{entity}', desired_cols, cols_to_group)
                    info = combine_dict_info(info, entity_info, old_info)
                if "entity: uniprot" not in info_not_needed:
                    uniprot_info = fetch_pdb_data('uniprot_entity',
                                                  f'{pdb}/{entity}', desired_cols, cols_to_group)
                    info = combine_dict_info(info, uniprot_info, old_info)

        if "instance" not in info_not_needed:
            old_info = list(info.keys())
            num_chains = int(info['rcsb_entry_info.deposited_nonpolymer_entity_instance_count'])+int(
                info['rcsb_entry_info.deposited_polymer_entity_instance_count'])
            chains = [chr(i) for i in range(65, 65+int(num_chains))]
            # these chains are pdb, author may have other, but pdb just list alphabetically.
            for chain in chains:
                chain_info = fetch_pdb_data('pdb_instance',
                                            f'{pdb}/{chain}', desired_cols, cols_to_group)
                info = combine_dict_info(info, chain_info, old_info)

        if "assembly" not in info_not_needed:
            old_info = list(info.keys())
            assemblies = get_list_from_table_entry(
                info['rcsb_entry_container_identifiers.assembly_ids'])
            for assembly in assemblies:
                assembly_info = fetch_pdb_data('assembly',
                                               f'{pdb}/{assembly}', desired_cols, cols_to_group)
                info = combine_dict_info(info, assembly_info, old_info)

        if ("chem" not in info_not_needed or "chem: drugbank" not in info_not_needed) and\
                ("pdbx_entity_nonpoly.comp_id" in info.keys()):
            old_info = list(info.keys())
            if int(info['rcsb_entry_info.nonpolymer_entity_count']) > 0:
                chems = get_list_from_table_entry(
                    info['pdbx_entity_nonpoly.comp_id'])
                for chem in chems:
                    if "chem" not in info_not_needed:
                        chem_info = fetch_pdb_data('chemcomp',
                                                   chem, desired_cols, cols_to_group)
                        info = combine_dict_info(info, chem_info, old_info)
                    if "chem: drugbank" not in info_not_needed:
                        drugbank_info = fetch_pdb_data('drugbank_chemcomp',
                                                       chem, desired_cols, cols_to_group)
                        info = combine_dict_info(
                            info, drugbank_info, old_info)

        all_info[pdb] = info
    df = pd.DataFrame.from_dict(all_info, orient='index')
    df['Entry ID'] = df.index
    df = data.merge(df, how="outer", on="Entry ID")
    if not output_dir:
        output_file = get_nonexistant_file(
            str(Path(ensemble_csv).parents[0])+'/summary_report_'+label+'.csv')
    else:
        output_file = output_dir + '/summary_report_'+label+'.csv'
    print(f"\nSaving summary table {output_file}. Continue to filter.")
    missing_cols = [col for col in desired_cols if col not in df.columns]
    if len(missing_cols) > 0:
        print("WARNING: No data found for", missing_cols)
    df = df.applymap(delete_all_one_value)
    df = df.rename(columns=desired_cols)
    df.to_csv(output_file, index=False)
    return df
