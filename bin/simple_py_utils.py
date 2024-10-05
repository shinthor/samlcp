#!/usr/bin/env python
# coding: utf-8
"""
Utility functions with no uncommon library dependencies
"""

class KeyDict(dict):
    """
    A dictionary that returns its key if the key is not present
    """
    def __missing__(self, key):
        return key

def str2bool(v):
    """Parse boolean string arguments for use by argparse command line inputs."""
    return str(v).lower() in ("yes", "true", "t", "y", "1")

def filter_threshold_combinations(threshold_combinations, var_names, var_colname=None):
    """
    Filter threshold combinations based on gene presence in var_names.
    Args:
    threshold_combinations (dict): Dictionary of threshold combinations.
    var_names (list): List of variable names.
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
                if not var_colname and gene not in var_names:  # Check if gene is not in var_names if var_colname is not provided
                    print(f"Warning: gene {gene} not found in the input data, skipping")  # Print a warning if gene is not found
                    continue  # Skip to the next iteration without including the combination group
                if var_colname and gene not in [str(x) for x in var_names]:  # Check if gene is not in var_names if var_colname is provided
                    print(f"Warning: gene {gene} not found in the input data, skipping")  # Print a warning if gene is not found
                    continue  # Skip to the next iteration without including the combination group
                filtered_combinations.append(criteria)  # Append the combination group to the filtered combinations if no genes were skipped
            else:
                filtered_combinations.append(criteria)
    return filtered_combinations  # Return the filtered threshold combinations
