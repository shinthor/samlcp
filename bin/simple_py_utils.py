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

def filter_threshold_combinations(threshold_combinations, adata, var_colname=None):
    """
    Filter threshold combinations based on gene presence in var_names.
    Args:
    threshold_combinations (dict): Dictionary of threshold combinations.
    adata (list): adata for which we want to verify contains the necessary variable names.
    var_colname (str, optional): Column name to use for var names. Defaults to None.
    Returns:
    list: Filtered threshold combinations.
    """
    if var_colname:
        var_names = set(adata.var[var_colname])  # Get variable names from var column
    else:
        var_names = set(adata.var_names)
    filtered_combinations = []  # Initialize an empty list to store filtered combinations
    for criteria in threshold_combinations:  # Iterate over each combination group
        for criteria_name, criteria_config in criteria.items():  # Iterate over each criteria config
            if criteria_config["type"] in {"gene", "bins"}:  # Check if the criteria type is "gene" or "bins"
                gene = None
                if criteria_config["type"] == "gene":
                    gene = criteria_name  # Set the gene to the criteria name if type is "gene"
                elif criteria_config["type"] == "bins":
                    gene = criteria_config["gene"]  # Set the gene to the criteria config gene if type is "bins"
                if gene is not None and gene not in var_names:  # Check if gene is not in var_names 
                    print(f"Warning: gene {gene} not found in the input data, skipping")  # Print a warning if gene is not found
                    continue  # Skip to the next iteration without including the combination group
                filtered_combinations.append(criteria)  # Append the combination group to the filtered combinations if no genes were skipped
            else:
                filtered_combinations.append(criteria)
    return filtered_combinations  # Return the filtered threshold combinations

def convert_new_to_old_threshold_combinations(new_form_threshold_combinations):
    """
    Convert new format threshold combinations to the old format.

    The new format consists of a dictionary with two keys: 'definitions' and 'combinations'.
    'definitions' contains the configuration details for each criteria, and 'combinations'
    contains a list of lists, where each inner list represents a combination of criteria
    names.

    The old format is a list of dictionaries, where each dictionary represents a combination
    of criteria with their respective configurations.

    Args:
    new_form_threshold_combinations (dict): A dictionary in the new format with keys 'definitions' and 'combinations'.

    Returns:
    list: A list of dictionaries in the old format, where each dictionary represents a combination of criteria.
    """
    # Extract definitions and combinations from the new threshold combinations format
    definitions = new_form_threshold_combinations['definitions']
    combinations = new_form_threshold_combinations['combinations']
    
    # Initialize the old format list
    old_format = []
    
    # Iterate through each combination
    for combination in combinations:
        # Create a new dictionary for this combination
        combination_dict = {}
        
        # Add each definition to the combination dictionary
        for config_name in combination:
            combination_dict[config_name] = definitions[config_name]
        
        # Add the combination dictionary to the old format list
        old_format.append(combination_dict)
    
    return old_format