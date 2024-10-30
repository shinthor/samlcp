#!/usr/bin/env python
# coding: utf-8
"""
This script preprocesses the data for graphing and saves the results to a tsv file and optionally
a copy of the original scanpy file with metadata columns added to adata.obs.
"""
import argparse
import os
import gc
from collections import namedtuple
import anndata
from scipy.sparse import issparse
import pyjson5 as json
import numpy as np
import pandas as pd
import scanpy as sc
from cachetools import LRUCache, cached, Cache
from cachetools.keys import hashkey
import sensig_score
import convert_species
import simple_py_utils
import py_pipe_utils as myutils

@cached( # Decorator to cache the outputs of this function in memory
        LRUCache(maxsize=128, getsizeof=Cache.getsizeof), # default getsizeof returns 1 for any object, so max size means number of items
        key=lambda adata, criteria_name, criteria_config, column2, homolog_table_path, aid_adata, aid_adata_normalized : hashkey(criteria_name, json.dumps(criteria_config))) # just cache based on the criteria and ignore the unhashable args
def preprocess_category(adata, criteria_name, criteria_config, column2, homolog_table_path, aid_adata, aid_adata_normalized):
    if criteria_config["type"] == "gene":
        curr_gene_idx = adata.var_names.get_loc(criteria_name)
        # Apply gene-based criteria
        return pd.Series(adata.X[:, curr_gene_idx] > criteria_config["threshold"],
                                                    index=adata.obs_names, copy=True).map(
                                                        {True: f"{criteria_name}_pos", False: f"{criteria_name}_neg"}), None
    elif criteria_config["type"] == "model":
        # Apply model-based criteria
        model_predictions = myutils.apply_model(adata, criteria_config, homolog_table_path=homolog_table_path)
        if "rename_classes" in criteria_config:
            # Make a dict from the rename_classes dict that renames classes included and keeps the rest the same
            map_dict = simple_py_utils.KeyDict(**criteria_config["rename_classes"])
            model_predictions = model_predictions.map(map_dict)
        if "overwrite_classes" in criteria_config:
            # Map the predictions to the desired classes. Unincluded classes will be set to NaNs so they will not override the original classes.
            return model_predictions.map(criteria_config["overwrite_classes"]), "overwrite_class"
        return model_predictions, None
    elif criteria_config["type"] == "model_select":
        # Apply model-based criteria
        model_predictions = myutils.apply_model(adata, criteria_config, homolog_table_path=homolog_table_path)
        # Use the predictions to keep only the specified classes and set the rest to NaN.
        model_predictions = model_predictions.map({x: "" for x in criteria_config["keep_classes"]})
        return model_predictions, None
    elif criteria_config["type"] == "bins":
        # Apply bin-based criteria
        curr_adata = adata
        curr_adata_index = curr_adata.obs.index
        adata.X = aid_adata_normalized.X[()]
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
        adata.X = aid_adata.X[()]
        return (criteria_name + "_lvl_" + pd.Series(bins_list, index=idx_list).astype(str).loc[curr_adata_index]), None
    elif criteria_config["type"] == "cell_cycle":
        return adata.obs["phase"].astype(str), None
    elif criteria_config["type"] == "sensig_score":
        curr_adata = adata
        curr_adata_index = curr_adata.obs.index
        scorer = sensig_score.gen_sensig_scorer(sensig_params=criteria_config["params"])
        scorer.sensig_score(curr_adata)
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
        return ("sensig_score" + pd.Series(bins_list, index=idx_list).astype(str).loc[curr_adata_index]), None
    elif criteria_config["type"] == "comparative_score":
        curr_adata = adata
        curr_adata_index = curr_adata.obs.index
        scorer = sensig_score.gen_comparative_scorer(scorer_main_params=criteria_config["main_params"], scorer_competitor_params=criteria_config["competitor_params"])
        name_base = criteria_config["main_params"]["new_score_column"]
        scorer.sensig_score(curr_adata)
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
        return (name_base + pd.Series(bins_list, index=idx_list).astype(str).loc[curr_adata_index]), None
def preprocess_data(input_file_path: str,
                    output_file_path: str,
                    threshold_combinations,
                    name_to_add,
                    sample_taxon,
                    resources,
                    column2="cell_ontology_class",
                    colnamestart="combination",
                    use_raw=False,
                    var_colname=None,
                    layer_name=None,
                    save_h5ad=False,
                    debug=False
) -> anndata.AnnData:
    """
    Preprocess the data for graphing and save the results to a tsv file and optionally
    a copy of the original scanpy file with metadata columns added to adata.obs.
    Parameters:
        input_file_path: path to the input .h5ad file
        output_file_path: path to the output .tsv file
        threshold_combinations: list of threshold combinations to use
        name_to_add: name to add to the .uns attribute
        sample_taxon: taxon of the sample
        resources: resources object
        column2: name of the column to use for the second column in the threshold combinations
        colnamestart: name to use for the first column in the threshold combinations
        use_raw: whether to use the raw data
        var_colname: name of the column to use for the variable names
        layer_name: name of the layer to use
        save_h5ad: whether to save a copy of the original scanpy file with metadata columns added to adata.obs
    Returns:
        adata: the original scanpy file with metadata columns added to adata.obs
    """
    # Load the input .h5ad file
    adata = anndata.read_h5ad(input_file_path)

    homolog_table_path = resources.homolog_table_path
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
        # Set old index to a column in the variable table
        adata.var['old_index'] = adata.var.index
        # Set var_colname as the new index
        adata.var_names = adata.var[var_colname].astype(str)
        # Get rid of the new index column name to avoid name conflicts
        adata.var.index.name = None
    
    # Filter the threshold combinations based on the input data
    filtered_combinations = simple_py_utils.filter_threshold_combinations(threshold_combinations, adata, var_colname)

    # Ensure the data is dense and has unique gene names
    if issparse(adata.X):
        adata.X = np.asarray(adata.X.todense())
    adata.var_names = adata.var_names.astype(str)
    gene_counts = adata.var_names.tolist()
    gene_counts_series = pd.Series(adata.var_names.copy(), index=adata.var_names.copy())
    duplicated_minus_first = gene_counts_series.duplicated(keep="first")
    duplicated_genes = gene_counts_series[duplicated_minus_first].unique()
    for gene in duplicated_genes:
        avg_expr = adata[:, adata.var_names == gene].X.mean(axis=1)
        gene_idx = np.where(adata.var_names == gene)[0]
        print("Gene idx:", gene_idx)
        adata.X[:, gene_idx] = avg_expr[:, np.newaxis]
    adata = adata[:, ~duplicated_minus_first]
    aid_adata_path = f'{name_to_add}_aid_adata.h5ad'
    aid_adata_normalized_path = f'{name_to_add}_aid_adata_normalized.h5ad'

    adata.write_h5ad(aid_adata_path)
    aid_adata = sc.read_h5ad(aid_adata_path, backed="r")
    sc.pp.normalize_total(adata, target_sum=1e4)
    adata.write_h5ad(aid_adata_normalized_path)
    aid_adata_normalized = sc.read_h5ad(aid_adata_normalized_path, backed="r")

    # Generate the cell cycle scores ahead of time
    sc.pp.log1p(adata, base=10)
    sc.pp.scale(adata)
    with open(resources.cell_cycle_genes_path, encoding="utf-8") as ccgf:
        cell_cycle_genes = [x.strip() for x in ccgf]
    s_genes = cell_cycle_genes[:43]
    g2m_genes = cell_cycle_genes[43:]
    if sample_taxon == 10090:
        # Convert human cell cycle genes to mouse cell cycle genes
        homolog_table = pd.read_csv(homolog_table_path, sep="\t")
        id2symbols, _, _ = convert_species.get_id2symbols_dict(homolog_table)
        compatible_genes, _ = convert_species.get_compatible_genes(id2symbols)
        s_genes = [compatible_genes[x][0] for x in s_genes if x in compatible_genes]
        g2m_genes = [compatible_genes[x][0] for x in g2m_genes if x in compatible_genes]
        if debug:
            print(f"s_genes converted: {s_genes}")
            print(f"g2m_genes converted: {g2m_genes}")
        
    s_genes = [x for x in s_genes if x in adata.var_names]
    g2m_genes = [x for x in g2m_genes if x in adata.var_names]
    if debug:
        print(f"s_genes overlap: {s_genes}")
        print(f"g2m_genes overlap: {g2m_genes}")
    sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes, use_raw=False)

    # Restore counts to adata.X
    adata.X = aid_adata.X[()]

    # Iterate over the threshold combinations
    original_adata = adata
    for criteria in filtered_combinations:
        try:
            # Generate a column name based on the combination and the criteria name
            column_name = f"{colnamestart}"
            curr_series = pd.Series("", index=adata.obs.index, dtype=str, copy=True)
            adata = original_adata
            for criteria_name, criteria_config in criteria.items():
                series_add, series_directive = preprocess_category(adata, criteria_name, criteria_config, column2, homolog_table_path, aid_adata, aid_adata_normalized)
                if series_directive == "overwrite_class":
                    curr_series[series_add.notnull()] = series_add
                else:
                    curr_series = curr_series.str.cat(series_add, sep="_")
                column_name += f"_{criteria_name}"
                gc.collect()
            adata.obs[column_name] = curr_series.str.strip("_")
        except Exception as e:
            print(f"Error in {criteria}: {e}")

    # Save the modified data to the output file path
    if output_file_path is not None:
        adata.var.index.name = None
        # Get full output_file_path without the extension
        output_file_path_no_ext = os.path.splitext(output_file_path)[0]
        # Save the adata.obs to tsv
        adata.obs.to_csv(output_file_path_no_ext + ".tsv.gz", sep="\t", compression="gzip")
        # Save the data to the output file path as full h5ad if the user wants to save it
        if save_h5ad:
            adata.write_h5ad(output_file_path_no_ext + ".h5ad", compression="gzip")
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
    parser.add_argument('--save-h5ad', type=simple_py_utils.str2bool, default=False,
                        nargs='?', const=True,
                        help='whether to save h5ad file')
    parser.add_argument('--debug', type=simple_py_utils.str2bool, default=False,
                        nargs='?', const=True,
                        help='Debug mode')
    args = parser.parse_args()

    with open(args.threshold_combinations, encoding="utf-8") as f:
        threshold_combinations_dict = json.load(f)
    if "definitions" in threshold_combinations_dict:
        threshold_combinations_dict = simple_py_utils.convert_new_to_old_threshold_combinations(threshold_combinations_dict)
    Resources = namedtuple('Resources', ['homolog_table_path', 'cell_cycle_genes_path'])
    gene_resources = Resources(args.homolog_table_path, args.cell_cycle_genes_path)
    preprocess_data(args.input_file,
                    args.output_file,
                    threshold_combinations_dict,
                    args.name_to_add,
                    args.sample_taxon,
                    gene_resources,
                    args.column2,
                    use_raw=args.use_raw,
                    var_colname=args.var_colname,
                    layer_name=args.layer_name,
                    save_h5ad=args.save_h5ad,
                    debug=args.debug
                    )
