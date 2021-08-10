""" pseudo_ensemble.refine.filter_ensemble

Functions to filter ensemble summary reports by organisms, keywords, mutations, resolution.

Authors:
    Rachael Kretsch (rkretsch@stanford.edu)
    Siyuan Du (dusiyuan@stanford.edu)
    Jacob Parres-Gold

Last edited:
    2021-08-10
"""

import numpy as np
import pandas as pd
from pathlib import Path

from EnsemblePDB.utils import file_management


def filter_pdb_data(summary_report_csv, protein_name, exclude_terms=None, max_res=None, organism=None, max_muts=None, output_directory=None):
    '''
    Takes a csv from a pdb search and filters according to user-specified
    metrics. Saves a csv with columns for the filter results and a csv
    with only the structures that passed all filters.
    Arguments: 
        summary_report_csv (str): The path to csv file of summary report
        protein_names (str): Name of the protein
        (optional)
        exclude_terms (list of strings): terms to exclude from entity descriptions {default: None}
        max_res (float): Maximum high-limit resolution to allow
        organism (list of strings): names of source organism {default: None}
        max_muts (int): maximum number of mutations allowed for the protein {default: None}
        refs_for_gaps = reference for sequence input should be input_# {default: None}
        output_directory (str): path to save filtered csvs to. Default saves to the same directory that summary_report_csv is in {default: None}.
    Returns:
        Dataframe, filtered summary report
    '''
    # Check user inputs
    assert(Path(summary_report_csv).is_file()), f"Summary report, "\
        f"{summary_report_csv}, is not a file."
    assert (isinstance(protein_name, str)), "protein_name must be a string"
    if exclude_terms is not None:
        assert (type(exclude_terms) ==
                list), "exclude_terms must be a list of string"
        assert (all(isinstance(s, str)
                    for s in exclude_terms)), "exclude terms must be strings"
    if max_res is not None:
        assert (isinstance(max_res, (int, float))), "max_res must be a number"
    if max_muts is not None:
        assert (isinstance(max_muts, (int, float))), "max_res must be a number"
    if organism is not None:
        assert (isinstance(organism, list)), "organism must be a list"

    # Import CSV
    summary_report = pd.read_csv(summary_report_csv)

    # Filters
    summary_report["Filtering out"] = False
    summary_report['Correct protein'] = summary_report.apply(
        lambda x: contains_keyword(x, protein_name), axis=1)
    summary_report['Filtering out'] = summary_report['Filtering out'] | (
        ~summary_report['Correct protein'])
    if exclude_terms is not None:
        for term in exclude_terms:
            summary_report['Entity description contains excluded words'] = summary_report.apply(
                lambda x: contains_keyword(x, term), axis=1)
            summary_report["Filtering out"] = summary_report['Filtering out'] | summary_report["Entity description contains excluded words"]
    if max_res is not None:
        summary_report["Too low resolution"] = summary_report.apply(
            lambda x: True if x['pdbx: resolution'] > max_res else False, axis=1)
        summary_report["Filtering out"] = summary_report['Filtering out'] | summary_report["Too low resolution"]
    if organism is not None:
        summary_report["Wrong organism"] = summary_report.apply(
            lambda x: wrong_organism(x, protein_name, organism), axis=1)
        summary_report["Filtering out"] = summary_report['Filtering out'] | summary_report["Wrong organism"]
    if max_muts is not None:
        summary_report["Too many mutations"] = summary_report.apply(
            lambda x: check_mutations(x, protein_name, max_muts), axis=1)
        summary_report["Filtering out"] = summary_report['Filtering out'] | summary_report["Too many mutations"]

    filtered = summary_report.loc[~summary_report['Filtering out']]\
        .drop(['Filtering out', 'Correct protein', 'Entity description contains excluded words',    'Too low resolution', 'Wrong organism', 'Too many mutations'], axis=1)
    if not output_directory:
        parent_dir = str(Path(summary_report_csv).parents[0])
        output_file_full = file_management.get_nonexistant_file(
            parent_dir+'/summary_report_'+protein_name+'_with_filter.csv')
        output_file = file_management.get_nonexistant_file(
            parent_dir+'/summary_report_'+protein_name+'_filtered.csv')
    else:
        output_file_full = file_management.get_nonexistant_file(
            output_directory+'/summary_report_'+protein_name+'_with_filter.csv')
        output_file = file_management.get_nonexistant_file(
            output_directory+'/summary_report_'+protein_name+'_filtered.csv')
    filtered.to_csv(output_file)
    summary_report.to_csv(output_file_full)
    print(f'Saved filtered summary report to {output_file}')
    print(f'Saved full summary report with filters to {output_file_full}')
    return filtered


def contains_keyword(x, name, col='rscb-polymer-entity: description', index=0):
    '''
    Apply on rows of dataframes, filter any entry (row) where description of chain does not contain the input content.
    Assumes the first polymer entity is always the major protein (not ligand).
    Case insentitive. 
    Arguments:
        x (row of dataframe)
        name (str): the keyword to filter with.
        col (str): column to filter the keyword on    
    Returns:
        bool
    '''
    words = x[col].split('~')[0].lower()
    if name.lower() in words:
        return True
    else:
        return False


def wrong_organism(x, name, orgs):
    '''
    Apply on rows of dataframes to check whether the protein is from a given organism
    Arguments:
        x (row of dataframe)
        name (str): protein name.
        org (list): organism scientific names
    Returns:
        bool
    '''
    organisms = x['pdbx-entity: organism scientific name'].lower().split('~')
    wrong_org = True
    for org in orgs:
        if org.lower() in organisms[0].lower():
            wrong_org = False
    return wrong_org


def check_mutations(x, name, max_muts):
    '''
    Apply on rows of dataframes to check whether the protein has mutation exceeding (>=) max mutation limit
    Arguments:
        x (row of dataframe)
        name (str): protein name.
        max_muts: maximum mutation number
    Returns:
        bool
    '''
    mutations = str(x['entity-polymer: mutation_count']).split('~')
    if int(mutations[0]) >= max_muts:
        return True
    else:
        return False
