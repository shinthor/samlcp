#!/usr/bin/env python
# coding: utf-8
"""
This script preprocesses the data for graphing and saves the results to a copy of the original
scanpy file with metadata columns added to adata.obs.
"""
import argparse
import pickle
import anndata
from scipy.sparse import issparse
from collections import namedtuple
import pyjson5 as json
import numpy as np
import pandas as pd
import scanpy as sc
import sensig_score
import convert_species


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
    for _, combination_group in enumerate(threshold_combinations):  # Iterate over each combination group
        for criteria in combination_group:  # Iterate over each criteria in the combination group
            for criteria_name, criteria_config in criteria.items():  # Iterate over each criteria config
                if criteria_config["type"] in {"gene", "bins"}:  # Check if the criteria type is "gene" or "bins"
                    if criteria_config["type"] == "gene":
                        gene = criteria_name  # Set the gene to the criteria name if type is "gene"
                    elif criteria_config["type"] == "bins":
                        gene = criteria_config["gene"]  # Set the gene to the criteria config gene if type is "bins"
                    if not var_colname and gene not in adata.var_names:  # Check if gene is not in adata var names if var_colname is not provided
                        print(f"Warning: gene {gene} not found in the input data, skipping")  # Print a warning if gene is not found
                        continue  # Skip to the next iteration without including the combination group
                    if var_colname and gene not in adata.var[var_colname]:  # Check if gene is not in adata var column if var_colname is provided
                        print(f"Warning: gene {gene} not found in the input data, skipping")  # Print a warning if gene is not found
                        continue  # Skip to the next iteration without including the combination group
                    else:
                        filtered_combinations.append(combination_group)  # Append the combination group to the filtered combinations if no genes were skipped

    return filtered_combinations  # Return the filtered threshold combinations

def preprocess_data(input_file_path,
                    output_file_path,
                    threshold_combinations,
                    name_to_add,
                    sample_taxon,
                    resources,
                    column2="cell_ontology_class",
                    colnamestart="combination",
                    use_raw=False,
                    var_colname=None,
                    layer_name=None
):
    # Load the input .h5ad file
    adata = anndata.read_h5ad(input_file_path)

    homolog_table_path = resources.homolog_table_path
    cell_cycle_genes_path = resources.cell_cycle_genes_path
    
    # Add the name of the dataset to the .uns attribute
    adata.uns['name'] = name_to_add
    adata.uns['sample_taxon'] = sample_taxon
    
    # If use_raw is True, use the raw data from the input file
    if use_raw:
        if adata.raw is not None:
            adata.X = adata.raw.X
            adata.var_names = adata.raw.var_names
        else:
            print("No raw data found in the input file but use_raw set to true. Using the original data.")
    
    # If a layer name is provided, use it to select the data layer
    if layer_name:
        adata.X = adata.layers[layer_name]
    
    # If a variable column name is provided, use it to select the variable names
    if var_colname:
        adata.var_names = adata.var[var_colname].astype(str)
    
    # Filter the threshold combinations based on the input data
    filtered_combinations = filter_threshold_combinations(threshold_combinations, adata, var_colname)
    # Ensure the data is dense and has unique gene names
    if issparse(adata.X):
        adata.X = np.asarray(adata.X.todense())
    adata.var_names = adata.var_names.astype(str)
    gene_counts = adata.var_names.tolist()
    gene_counts_series = pd.Series(gene_counts)
    duplicated_minus_first = gene_counts_series.duplicated(keep="first")
    duplicated_genes = gene_counts_series[duplicated_minus_first].unique()
    for gene in duplicated_genes:
        avg_expr = adata[:, adata.var_names == gene].X.mean(axis=1)
        gene_idx = np.where(adata.var_names == gene)[0]
        print("Gene idx:", gene_idx)
        adata.X[:, gene_idx] = avg_expr[:, np.newaxis]
    adata = adata[:, ~duplicated_minus_first]
    
    adata.layers["counts_pipe"] = adata.X
    sc.pp.normalize_total(adata, target_sum=1e4)
    adata.layers["norm_pipe"] = adata.X
    sc.pp.log1p(adata, base=2)
    adata.layers["log_pipe"] = adata.X
    sc.pp.scale(adata)
    adata.layers["scale_pipe"] = adata.X


    # Iterate over the threshold combinations
    adata_df = adata.to_df()
    for _, combination_group in enumerate(filtered_combinations):
        # Generate a column name based on the combination and the criteria name
        column_name = f"{colnamestart}"
        curr_series = pd.Series("", index=adata_df.index, dtype=str)
        for criteria in combination_group:
            for criteria_name, criteria_config in criteria.items():
                if criteria_config["type"] == "gene":
                    # Apply gene-based criteria
                    curr_series = curr_series.str.cat((adata_df[criteria_name] > criteria_config["threshold"]).map({True: f"{criteria_name}_pos", False: f"{criteria_name}_neg"}), sep="_")
                elif criteria_config["type"] == "model":
                    # Apply model-based criteria
                    curr_series = curr_series.str.cat(apply_model(adata, criteria_config, homolog_table_path=homolog_table_path), sep="_")
                elif criteria_config["type"] == "bins":
                    # Apply bin-based criteria
                    curr_adata = adata
                    curr_adata_index = curr_adata.obs.index
                    sc.pp.normalize_total(curr_adata, target_sum=1e4)
                    unique_col2_vals = curr_adata.obs[column2].unique()
                    gene_subset = curr_adata[:, criteria_config["gene"]]
                    idx_list = []
                    bins_list = []
                    for curr_col2_val in unique_col2_vals:
                        curr_col2_subset = gene_subset[gene_subset.obs[column2] == curr_col2_val]
                        curr_col2_subset_idx = curr_col2_subset.obs.index
                        curr_exp = curr_col2_subset.X.flatten()
                        upper_val = np.percentile(curr_exp, criteria_config["upper_limit_percentile"])
                        if upper_val == 0:
                            upper_val = np.max(curr_exp)
                        lower_val = np.percentile(curr_exp, criteria_config["lower_limit_percentile"])
                        if lower_val == 0 and upper_val > 0:
                            lower_val = np.min(curr_exp[np.nonzero(curr_exp)])
                        bin_edges = np.linspace(lower_val, upper_val, criteria_config["num_bins"])
                        curr_bins = np.digitize(curr_exp, bin_edges)
                        bins_list += list(curr_bins)
                        idx_list += list(curr_col2_subset_idx)
                    curr_series = curr_series.str.cat(criteria_name + "_lvl_" + pd.Series(bins_list, index=idx_list).astype(str).loc[curr_adata_index], sep="_")
                elif criteria_config["type"] == "cell_cycle":
                    # Apply cell cycle-based criteria
                    with open(cell_cycle_genes_path, encoding="utf-8") as ccgf:
                        cell_cycle_genes = [x.strip() for x in ccgf]
                    s_genes = cell_cycle_genes[:43]
                    g2m_genes = cell_cycle_genes[43:]
                    if sample_taxon == 10090:
                        # Convert human cell cycle genes to mouse cell cycle genes
                        homolog_table = pd.read_csv(homolog_table_path, sep="\t")
                        id2symbols, _, _ = convert_species.get_id2symbols_dict(homolog_table)
                        compatible_genes, reverse_translate_dict = convert_species.get_compatible_genes(id2symbols)
                        s_genes = [compatible_genes[x][0] for x in s_genes if x in compatible_genes]
                        g2m_genes = [compatible_genes[x][0] for x in g2m_genes if x in compatible_genes]
                    s_genes = [x for x in s_genes if x in adata.var_names]
                    g2m_genes = [x for x in g2m_genes if x in adata.var_names]
                    adata_copy = adata.copy()
                    sc.pp.normalize_total(adata_copy, target_sum=1e4)
                    sc.pp.log1p(adata_copy, base=2)
                    sc.pp.scale(adata_copy)
                    sc.tl.score_genes_cell_cycle(adata_copy, s_genes=s_genes, g2m_genes=g2m_genes)
                    curr_series = curr_series.str.cat(adata_copy.obs["phase"].astype(str), sep="_")
                elif criteria_config["type"] == "sensig_score":
                    curr_adata = adata.copy()
                    curr_adata_index = curr_adata.obs.index
                    scorer = sensig_score.gen_sensig_scorer(taxon_id=sample_taxon)
                    sensig_scores = scorer.sensig_score(curr_adata)
                    unique_col2_vals = curr_adata.obs[column2].unique()
                    idx_list = []
                    bins_list = []
                    for curr_col2_val in unique_col2_vals:
                        curr_col2_subset = curr_adata[curr_adata.obs[column2] == curr_col2_val]
                        curr_col2_subset_idx = curr_col2_subset.obs.index
                        lower_val = np.percentile(curr_col2_subset.obs["sensig_score"], criteria_config["lower_limit_percentile"])
                        upper_val = np.percentile(curr_col2_subset.obs["sensig_score"], criteria_config["upper_limit_percentile"])
                        bin_edges = np.linspace(lower_val, upper_val, criteria_config["num_bins"])
                        curr_bins = np.digitize(curr_col2_subset.obs["sensig_score"], bin_edges)
                        bins_list += list(curr_bins)
                        idx_list += list(curr_col2_subset_idx)
                    curr_series = curr_series.str.cat("sensig_score" + pd.Series(bins_list, index=idx_list).astype(str).loc[curr_adata_index], sep="_")
                elif criteria_config["type"] == "comparative_score":
                    curr_adata = adata.copy()
                    curr_adata_index = curr_adata.obs.index
                    if "comparison_name" in criteria_config:
                        scorer = sensig_score.gen_comparative_scorer(taxon_id=sample_taxon, comparison_name=criteria_config["comparison_name"])
                        name_base = "{}_score".format(criteria_config["comparison_name"])
                    else:
                        scorer = sensig_score.gen_comparative_scorer(taxon_id=sample_taxon)
                        name_base = "paracrine_score"
                    sensig_scores = scorer.sensig_score(curr_adata)
                    unique_col2_vals = curr_adata.obs[column2].unique()
                    idx_list = []
                    bins_list = []
                    for curr_col2_val in unique_col2_vals:
                        curr_col2_subset = curr_adata[curr_adata.obs[column2] == curr_col2_val]
                        curr_col2_subset_idx = curr_col2_subset.obs.index
                        lower_val = np.percentile(curr_col2_subset.obs[name_base], criteria_config["lower_limit_percentile"])
                        upper_val = np.percentile(curr_col2_subset.obs[name_base], criteria_config["upper_limit_percentile"])
                        bin_edges = np.linspace(lower_val, upper_val, criteria_config["num_bins"])
                        curr_bins = np.digitize(curr_col2_subset.obs[name_base], bin_edges)
                        bins_list += list(curr_bins)
                        idx_list += list(curr_col2_subset_idx)
                    curr_series = curr_series.str.cat(name_base + pd.Series(bins_list, index=idx_list).astype(str).loc[curr_adata_index], sep="_")
                column_name += f"_{criteria_name}"
        adata.obs[column_name] = curr_series.str.lstrip("_")

    # Save the modified data to the output file path
    if output_file_path is not None:
        adata.write_h5ad(output_file_path)
    return adata


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process some data.')
    parser.add_argument('--input-file', type=str, help='path to input file')
    parser.add_argument('--threshold-combinations', type=str,
                        help='path to threshold combinations JSON file')
    parser.add_argument('--output-file', type=str, default='processed_data.h5ad',
                        help='path to output file')
    parser.add_argument('--name-to-add', type=str,
                        help='name to add to the processed data')
    parser.add_argument('--sample-taxon', type=int,
                        help='sample taxon ID')
    parser.add_argument('--homolog-table-path', type=str,
                        help='path to homolog table file HOM_MouseHumanSequence.rpt from Jackson labs')
    parser.add_argument('--cell-cycle-genes-path', type=str,
                        help='path to cell cycle genes file from Aviv Regev lab')
    parser.add_argument('--column2', type=str, default='cell_ontology_class', help='column to use for second level of grouping')
    parser.add_argument('--use-raw', default=False, help='use scanpy raw data', type=bool)
    parser.add_argument('--var-colname', type=str, default=None, help='column name to use for var names')
    parser.add_argument('--layer-name', type=str, default=None, help='layer name to use for data')
    args = parser.parse_args()
    
    with open(args.threshold_combination, encoding="utf-8") as f:
        threshold_combinations = json.load(f)
    
    Resources = namedtuple('Resources', ['homolog_table_path', 'cell_cycle_genes_path'])
    resources = Resources(args.homolog_table_path, args.cell_cycle_genes_path)
    preprocess_data(args.input_file,
                    args.output_file,
                    threshold_combinations,
                    args.name_to_add,
                    args.sample_taxon,
                    resources,
                    args.column2,
                    use_raw=args.use_raw,
                    var_colname=args.var_colname,
                    layer_name=args.layer_name
                    )
