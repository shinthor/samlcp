#!/usr/bin/env python
# coding: utf-8
"""
This module defines classes and methods to generate "Senescence Signature scores"
(SenSig scores) for a given single-cell RNA sequencing dataset.
This is an implementation of the logic described in the paper:
"Transfer learning in a biomaterial fibrosis model identifies in vivo senescence heterogeneity and contributions to vascularization and matrix production across species and diverse pathologies"
(Cherry et al, 2023, GeroScience)
Higher scores mean the cell is more senescence-like.
"""
import numpy as np
import pandas as pd
import anndata as ad
from typing import Dict

class SensigScorer():
    """
    A scorer class to generate Senescence Signature scores.
    
    Parameters:
    sensig_df (pd.DataFrame): Differential expression analysis between senescent and non-senescent cells.
    threshold (float, optional): Threshold for filtering genes. Defaults to None.
    filter_column (str, optional): Column to filter genes. Defaults to "FDR".
    gene_symbol_column (str, optional): Column for gene symbols. Defaults to "HGNC.symbol".
    log_fc_column (str, optional): Column for log fold change. Defaults to "logFC".
    adata_new_column_default (str, optional): Default column name for AnnData object. Defaults to "sensig_score".
    """
    def __init__(self, sensig_df, threshold=None, filter_column="FDR", gene_symbol_column="HGNC.symbol", log_fc_column="logFC", adata_new_column_default="sensig_score"):
        self.log_fc_column = log_fc_column
        self.adata_new_column_default = adata_new_column_default
        sensig_df = self.filter_threshold(sensig_df, threshold=threshold, filter_column=filter_column)
        sensig_df = self.resolve_duplicates(sensig_df, gene_symbol_column=gene_symbol_column)
        self.sensig_sign = np.sign(sensig_df[self.log_fc_column])
        self.sensig_sign.index=sensig_df[gene_symbol_column]
        self.sensig_sign.dropna(inplace=True)
    def filter_threshold(self, sensig_df, threshold=0.05, filter_column="FDR"):
        if threshold == None:
            return sensig_df
        else:
            return sensig_df[sensig_df[filter_column]<threshold]
    def resolve_duplicates(self, sensig_df, gene_symbol_column="HGNC.symbol"):
        duplicate_df = sensig_df[sensig_df[gene_symbol_column].duplicated(keep=False)]
        no_duplicate_df = sensig_df.drop_duplicates(subset=gene_symbol_column, keep=False, inplace=False)
        unique_duplicated_genes = duplicate_df[gene_symbol_column].unique()
        rows_to_add = []  # Accumulator list for rows to be added
        for gene in unique_duplicated_genes:
            gene_df = duplicate_df[duplicate_df[gene_symbol_column] == gene]
            # Remove this check if you want to include genes with only one row
            if len(gene_df) > 0:  # Ensure there are rows for this gene
                gene_df = gene_df.sort_values(by=self.log_fc_column, ascending=False)
                row_to_add = gene_df.iloc[0]
                row_to_add[self.log_fc_column] = gene_df[self.log_fc_column].mean()
                rows_to_add.append(row_to_add) 
            # Concatenate the list of rows to the DataFrame
        if rows_to_add:
            no_duplicate_df = pd.concat([no_duplicate_df] + rows_to_add, axis=0, ignore_index=True)

        return no_duplicate_df
    def get_sensig_score_parts(self, query_df):
        sensig_genes_mask = np.isin(query_df.columns, self.sensig_sign.index)
        stdev = query_df.std(ddof=1)
        # degrees of freedom, 1, this is the (n-1) in the denominator for std which is appropriate for a "sample std" rather than population std.
        means = query_df.mean(axis=0)
        all_z_scores = ((query_df - means)/stdev)
        sensig_z_scores = all_z_scores.loc[:, all_z_scores.columns[sensig_genes_mask]]
        sensig_sign = self.sensig_sign.loc[all_z_scores.columns[sensig_genes_mask]]
        return sensig_z_scores.multiply(sensig_sign)
    def sensig_score(self, query_df, adata_new_column=None):
        """
        Calculate Senescence Signature scores for a given query dataset.
        
        Parameters:
        query_df (pd.DataFrame): Query dataset to calculate Senescence Signature scores.
        adata_new_column (str, optional): Column name for AnnData object. Defaults to None.
        
        Returns:
        pd.Series: Senescence Signature scores.
        """
        is_adata = False
        query_adata = None
        if isinstance(query_df, ad.AnnData):
            query_adata = query_df
            query_df = query_df.to_df()
            is_adata = True
        if adata_new_column is None:
            adata_new_column = self.adata_new_column_default
        sensig_genes_mask = np.isin(query_df.columns, self.sensig_sign.index)
        stdev = query_df.std(ddof=1)
        # degrees of freedom, 1, this is the (n-1) in the denominator for std which is appropriate for a "sample std" rather than population std.
        means = query_df.mean(axis=0)
        all_z_scores = ((query_df - means)/stdev)
        sensig_z_scores = all_z_scores.loc[:, all_z_scores.columns[sensig_genes_mask]]
        sensig_sign = self.sensig_sign.loc[all_z_scores.columns[sensig_genes_mask]]
        sensig_score = sensig_z_scores.multiply(sensig_sign).sum(axis=1)
        if is_adata:
            query_adata.obs[adata_new_column] = sensig_score
        return sensig_score

class ComparativeScorer(SensigScorer):
    """
    A comparative scorer class to generate Senescence Signature scores relative to another signature.
    
    Parameters:
    sensig_df (pd.DataFrame): Differential expression analysis between senescent and non-senescent cells.
    competitor_scorer (SensigScorer): Another Senescence Signature scorer to compare with.
    threshold (float, optional): Threshold for filtering genes. Defaults to None.
    filter_column (str, optional): Column to filter genes. Defaults to "FDR".
    gene_symbol_column (str, optional): Column for gene symbols. Defaults to "HGNC.symbol".
    log_fc_column (str, optional): Column for log fold change. Defaults to "logFC".
    adata_new_column_default (str, optional): Default column name for AnnData object. Defaults to "notch_score".
    """
    def __init__(self, sensig_df, competitor_scorer: SensigScorer, threshold=None, filter_column="FDR", gene_symbol_column="HGNC.symbol", log_fc_column="logFC", adata_new_column_default="notch_score"):
        super().__init__(sensig_df, threshold=threshold, filter_column=filter_column, gene_symbol_column=gene_symbol_column, log_fc_column=log_fc_column, adata_new_column_default=adata_new_column_default)
        # Remove the score for genes that are in the competitor_scorer.sensig_sign with the same sign
        # The remaining genes should be genes that are not in common between the scorers, or are opposite signs
        mask_non_overlap = ~self.sensig_sign.index.isin(competitor_scorer.sensig_sign.index)
        overlap_index = self.sensig_sign.index[self.sensig_sign.index.isin(competitor_scorer.sensig_sign.index)]
        mask_opposite_sign = (self.sensig_sign.loc[overlap_index] != competitor_scorer.sensig_sign.loc[overlap_index]).values
        opposite_sign_index = overlap_index[mask_opposite_sign]
        non_overlap_index = self.sensig_sign.index[mask_non_overlap]
        # Use union of opposite_sign_index and non_overlap_index to set new self.sensig_sign
        new_index = opposite_sign_index.union(non_overlap_index)
        self.sensig_sign = self.sensig_sign.loc[new_index]

def gen_sensig_scorer(sensig_params: Dict[str, str]) -> SensigScorer:
    sensig_df = pd.read_csv(sensig_params["filepath"], sep="\t")
    return SensigScorer(sensig_df,
                        threshold=sensig_params["threshold"],
                        filter_column=sensig_params["filter_column"],
                        gene_symbol_column=sensig_params["gene_symbol_column"],
                        log_fc_column=sensig_params["log_fc_column"],
                        adata_new_column_default=sensig_params["new_score_column"]
                        )

def gen_comparative_scorer(
                           scorer_main_params: Dict[str, str],
                           scorer_competitor_params: Dict[str, str]
                           ) -> ComparativeScorer:
    main_df = pd.read_csv(scorer_main_params["filepath"], sep="\t", index_col=0)
    competitor_df = pd.read_csv(scorer_competitor_params["filepath"], sep="\t", index_col=0)
    competitor_scorer = SensigScorer(
        competitor_df,
        threshold=scorer_competitor_params["threshold"],
        filter_column=scorer_competitor_params["filter_column"],
        gene_symbol_column=scorer_competitor_params["gene_symbol_column"],
        log_fc_column=scorer_competitor_params["log_fc_column"])
    return ComparativeScorer(
        main_df,
        competitor_scorer,
        threshold=scorer_main_params["threshold"],
        filter_column=scorer_main_params["filter_column"],
        gene_symbol_column=scorer_main_params["gene_symbol_column"],
        log_fc_column=scorer_main_params["log_fc_column"],
        adata_new_column_default=scorer_main_params["new_score_column"])