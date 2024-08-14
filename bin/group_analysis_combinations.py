#!/usr/bin/env python
# coding: utf-8
"""
Add gene combinations from an AnnData object and prepare data for graphing.

This script takes an AnnData object as input, takes gene combinations based on column names,
and prepares data for graphing by calculating proportions of each level3 category for each
combination of level1 category and level2 category.

The resulting dataframes are saved as TSV files and can be used for graphing in a separate script.
"""
import argparse
import os
import pandas as pd
import scanpy as sc


def process_data_for_graphing(graph_df: pd.DataFrame,
                              level1_category: str = "Age",
                              level2_category: str = "Cell Type",
                               level3_category: str = "Predicted Senescence Class") -> pd.DataFrame:
    """
    Groups data by level1_category, level2_category, level3_category and returns a dataframe with the counts of each
    level2_category and level3_category to be used for plotting a 3 category grid of pie charts.
    """
    # Calculate unique counts for each combination of level1_category, level2_category, and level3_category
    graph_df_unique_counts = graph_df.groupby(
        [level1_category, level2_category, level3_category]).size().reset_index().rename(columns={0: "count"})
    print("graph_df_unique_counts")
    print(graph_df_unique_counts.head())
    # Calculate the proportion of each level2_category for each level1_category
    count_age_celltype_df = graph_df.groupby([level1_category, level2_category]).agg(
        {level2_category: "count"}).rename(columns={level2_category: "Count"}).reset_index()
    count_age_celltype_df[f"Proportion_{level2_category}_Per_{level1_category}"] = count_age_celltype_df["Count"] / \
                                                                                    count_age_celltype_df.groupby(
                                                                                        level1_category,
                                                                                        axis="index")[
                                                                                        "Count"].transform("sum")
    print("count_age_celltype_df")
    print(count_age_celltype_df.head())
    # Set multiindex for count_age_celltype_df
    multiindex = pd.MultiIndex.from_frame(count_age_celltype_df[[level1_category, level2_category]])
    count_age_celltype_df = count_age_celltype_df.drop(columns=[level1_category, level2_category]).set_index(multiindex)

    # Set multiindex for graph_df_unique_counts
    multiindex_main_grouper = pd.MultiIndex.from_frame(
        graph_df_unique_counts[[level1_category, level2_category, level3_category]])
    graph_df_unique_counts = graph_df_unique_counts.drop(
        columns=[level1_category, level2_category, level3_category]).set_index(multiindex_main_grouper)
    graph_df_unique_counts["Proportion_{1}_Per_{0}".format(level1_category, level2_category)]  = count_age_celltype_df["Proportion_{1}_Per_{0}".format(level1_category, level2_category)] 

    # Calculate proportion of each level3_category for each combination of level1_category and level2_category
    graph_df_unique_counts[f"Proportion_{level3_category}_Per_{level1_category}_AND_{level2_category}_Combination"] = \
        graph_df_unique_counts["count"] / graph_df_unique_counts.groupby(
            [level1_category, level2_category], axis="index")["count"].transform("sum")

    graph_df_unique_counts = graph_df_unique_counts.fillna(0)
    return graph_df_unique_counts

def age_column_transform(age: str) -> int:
    """Find number in string and return it as an int"""
    return int("".join(filter(str.isdigit, age)))
def main(adata: sc.AnnData,
         output_path_base: str,
         colnamestart: str="combination_",
         level1_category: str = "age",
         level2_category: str = "cell_ontology_class",
         name_to_add: str = "Tabula"):
    """
    Create dataframes with the proportion of each level3_category for each combination of level1_category and level2_category,
    useful for graphing in a separate script.
    """
    # Define gene combinations based on the column names added during preprocessing
    gene_combinations = [col for col in adata.obs.columns if col.startswith(colnamestart)]

    # Use process_data_for_graphing with the new columns added
    for column_name in gene_combinations:
        curr_obs = adata.obs
        # Add mask to drop rows where curr_obs[column_name] is null
        obs_mask = adata.obs[column_name].notnull()
        df_dict = {
            level1_category: curr_obs[level1_category][obs_mask],
            level2_category: curr_obs[level2_category][obs_mask],
            column_name[len(colnamestart):]: curr_obs[column_name][obs_mask]
        }
        if level1_category == "age" or level1_category == "development_stage":
            df_dict[level1_category] = df_dict[level1_category].apply(age_column_transform)
        out_df = process_data_for_graphing(pd.DataFrame(df_dict), level1_category=level1_category, level2_category=level2_category, level3_category=column_name[len(colnamestart):])
        out_df.to_csv(
            os.path.join(output_path_base,
                         name_to_add + "_" + column_name[len(colnamestart):]+"_gene_combinations.tsv"), sep="\t"
        )

if __name__ == "__main__":
    # Args from nextflow are $processed_data "${params.output_dir}/${input_file_name}/gene_combinations/${input_file_name}_"
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_path", help="Path to the input file")
    parser.add_argument("--output_path_base", help="Base path for output files")
    parser.add_argument("--uns_name", help="Name to add to output files")
    parser.add_argument("--level1_category", help="Level 1 category for graphing")
    parser.add_argument("--level2_category", help="Level 2 category for graphing")
    args = parser.parse_args()

    adata = sc.read_h5ad(args.input_path)
    output_path_base = args.output_path_base
    uns_name = args.uns_name
    level1_category = args.level1_category
    level2_category = args.level2_category
    # Load the AnnData object from the output of the preprocessing step
    main(adata, output_path_base, level1_category = level1_category, level2_category=level2_category, name_to_add=uns_name)