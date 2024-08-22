#!/usr/bin/env python
# coding: utf-8
"""
Module for useful python functions/classes used in the pipeline
"""
import pickle
import scanpy as sc
import anndata
import pandas as pd
import convert_species

class KeyDict(dict):
    """
    A dictionary that returns its key if the key is not present
    """
    def __missing__(self, key):
        return key

def str2bool(v):
    """Parse boolean string arguments for use by argparse command line inputs."""
    return str(v).lower() in ("yes", "true", "t", "y", "1")

def apply_model(adata: anndata.AnnData,
                model_config: str,
                homolog_table_path: str):
    """
    Apply a machine learning model to the given AnnData object.

    Parameters:
        adata (AnnData): The input AnnData object.
        model_config (dict): Configuration for the machine learning model.
        homolog_table_path (str): Path to the homolog table file.

    Returns:
        pd.Series: The predicted values for each cell.
    """
    # Load the machine learning model from the configuration
    with open(model_config['path'], 'rb') as model_file:
        clf = pickle.load(model_file)

    # Get the taxon ID from the model configuration
    model_taxon = model_config['taxon']

    # Check if the model taxon matches the sample taxon
    if model_taxon != adata.uns['sample_taxon']:
        print("Warning: model taxon is not the same as the sample taxon, converting sample taxon to model taxon for prediction")

        # Load the homolog table
        homolog_table = pd.read_csv(homolog_table_path, sep="\t")

        # Get the ID-to-symbols dictionary and compatible genes
        id2symbols, _, _ = convert_species.get_id2symbols_dict(homolog_table)
        compatible_genes, reverse_translate_dict = convert_species.get_compatible_genes(id2symbols)

        # Get the one-to-one dictionary for gene conversion
        one2one_dict = convert_species.get_one_to_one_only(reverse_translate_dict, compatible_genes)

        # Convert the AnnData object to the model taxon
        adata = convert_species.convert_species_one2one_adata(adata, one2one_dict, target_genes=clf.feature_names_in_)
    else:
        # Get the genes that are not in the AnnData object
        genes_not_in_adata = set(clf.feature_names_in_).difference(set(adata.var_names))

        # Create a new AnnData object with the missing genes
        adata_df = adata.to_df()[list(set(adata.var_names).intersection(set(clf.feature_names_in_)))]
        adata_df[list(genes_not_in_adata)] = 0
        adata_df = adata_df[clf.feature_names_in_]
        new_adata = anndata.AnnData(adata_df)
        adata = new_adata

    # Normalize and log-transform the AnnData object
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata, base=2)

    # Make predictions using the machine learning model
    return pd.Series(clf.predict(adata.to_df()), index=adata.obs_names.copy())

def filter_threshold_combinations(threshold_combinations, adata, var_colname=None):
    """
    Filter threshold combinations based on gene presence in adata.

    Args:
        threshold_combinations (dict): Dictionary of threshold combinations.
        adata (AnnData): AnnData object.
        var_colname (str, optional): Column name to use for var names. Defaults to None.

    Returns:
        list: Filtered threshold combinations.
    """
    filtered_combinations = []  # Initialize an empty list to store filtered combinations
    for criteria in threshold_combinations:  # Iterate over each combination group
        for criteria_name, criteria_config in criteria.items():  # Iterate over each criteria config
            if criteria_config["type"] in {"gene", "bins"}:  # Check if the criteria type is "gene" or "bins"
                if criteria_config["type"] == "gene":
                    gene = criteria_name  # Set the gene to the criteria name if type is "gene"
                elif criteria_config["type"] == "bins":
                    gene = criteria_config["gene"]  # Set the gene to the criteria config gene if type is "bins"
                if not var_colname and gene not in adata.var_names:  # Check if gene is not in adata var names if var_colname is not provided
                    print(f"Warning: gene {gene} not found in the input data, skipping")  # Print a warning if gene is not found
                    continue  # Skip to the next iteration without including the combination group
                if var_colname and gene not in adata.var[var_colname].astype(str):  # Check if gene is not in adata var column if var_colname is provided
                    print(f"Warning: gene {gene} not found in the input data, skipping")  # Print a warning if gene is not found
                    continue  # Skip to the next iteration without including the combination group
                filtered_combinations.append(criteria)  # Append the combination group to the filtered combinations if no genes were skipped
            # May want to add checks for other types here
            else:
                filtered_combinations.append(criteria)
                        

    return filtered_combinations  # Return the filtered threshold combinations
