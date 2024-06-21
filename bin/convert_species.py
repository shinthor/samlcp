import anndata
import pandas as pd
from collections import defaultdict 

def get_id2symbols_dict(homolog_table,
                         entry_id_column="DB Class Key",
                         symbol_column="Symbol",
                         taxon_id_column="NCBI Taxon ID"):
    """
    Create dictionaries mapping class keys to symbols and sizes.
    Args:
        homolog_table (pandas DataFrame): Table containing homologous genes.
        entry_id_column (str, optional): Column name for class keys. Defaults to "DB Class Key".
        symbol_column (str, optional): Column name for gene symbols. Defaults to "Symbol".
        taxon_id_column (str, optional): Column name for taxon IDs. Defaults to "NCBI Taxon ID".
    Returns:
        id2symbols (dict): Dictionary mapping class keys to dictionaries of taxon IDs to lists of symbols.
        id2sizes (dict): Dictionary mapping class keys to dictionaries of taxon IDs to symbol counts.
        max_size (int): Maximum symbol count across all taxons and class keys.
    """
    id2symbols = {}
    id2sizes = {}
    max_size = 0

    # Extract columns from homolog table
    class_keys = homolog_table[entry_id_column]
    taxon_ids = homolog_table[taxon_id_column]
    symbols = homolog_table[symbol_column]

    # Iterate over each row in the homolog table
    for i in range(homolog_table.shape[0]):
        curr_key = class_keys.iloc[i]
        curr_taxon = taxon_ids.iloc[i]
        curr_symbol = symbols.iloc[i]

        # Update id2symbols and id2sizes dictionaries
        if curr_key in id2symbols:
            if curr_taxon in id2symbols[curr_key]:
                id2symbols[curr_key][curr_taxon].append(curr_symbol)
                id2sizes[curr_key][curr_taxon] += 1
                id2sizes[curr_key]["total"] += 1
                if id2sizes[curr_key][curr_taxon] > max_size:
                    max_size = id2sizes[curr_key][curr_taxon]
            else:
                id2symbols[curr_key][curr_taxon] = [curr_symbol]
                id2sizes[curr_key][curr_taxon] = 1
                id2sizes[curr_key]["total"] += 1
        else:
            id2symbols[curr_key] = {curr_taxon: [curr_symbol]}
            id2sizes[curr_key] = {"total": 1, curr_taxon: 1}

    return id2symbols, id2sizes, max_size

def get_compatible_genes(id2symbols, taxon_we_want=9606, compatibility_taxon=10090):
    """
    Get compatible genes between two species.
    Args:
        id2symbols (dict): Dictionary mapping class keys to dictionaries of taxon IDs to lists of symbols.
        taxon_we_want (int, optional): Taxon ID of the species we want to convert to. Defaults to 9606.
        compatibility_taxon (int, optional): Taxon ID of the species we are converting from. Defaults to 10090.
    Returns:
        compatible_genes (dict): Dictionary mapping symbols to lists of compatible genes.
        reverse_translate_dict (dict): Dictionary mapping compatible genes to lists of symbols.
    """
    compatible_genes = dict()
    reverse_translate_dict = dict()
    for curr_id in id2symbols:
        if taxon_we_want in id2symbols[curr_id] and compatibility_taxon in id2symbols[curr_id]:
            for curr_symbol in id2symbols[curr_id][taxon_we_want]:
                if curr_symbol in compatible_genes:
                    compatible_genes[curr_symbol] += id2symbols[curr_id][compatibility_taxon]
                else:
                    compatible_genes[curr_symbol] = id2symbols[curr_id][compatibility_taxon]
            for curr_symbol in id2symbols[curr_id][compatibility_taxon]:
                if curr_symbol in reverse_translate_dict:
                    reverse_translate_dict[curr_symbol] += id2symbols[curr_id][taxon_we_want]
                else:
                    reverse_translate_dict[curr_symbol] = id2symbols[curr_id][taxon_we_want]
    return compatible_genes, reverse_translate_dict

def get_one_to_one_only(compatible_genes, reverse_translate_dict):
    """
    Get a one-to-one dictionary of compatible genes.
    Args:
        compatible_genes (dict): Dictionary mapping symbols to lists of compatible genes.
        reverse_translate_dict (dict): Dictionary mapping compatible genes to lists of symbols.

    Returns:
        one2one_dict (dict): Dictionary mapping one-to-one compatible genes.
    """
    one2one_dict = dict()
    for curr_gene in compatible_genes:
        if len(compatible_genes[curr_gene]) <= 1 and len(reverse_translate_dict[compatible_genes[curr_gene][0]]) <= 1:
            one2one_dict[curr_gene] = compatible_genes[curr_gene][0]
    return one2one_dict

def convert_species_one2one_adata(adata, one2one_dict, 
                                 old_var_name="MouseGeneSymbol", 
                                 new_var_name="HumanGeneSymbol", 
                                 target_genes=None,
                                 var_names_column=None) -> anndata.AnnData:
    """
    Convert the species of an AnnData object using a one-to-one dictionary.

    Args:
        adata (anndata.AnnData): AnnData object to convert.
        one2one_dict (dict): One-to-one dictionary of gene conversions.
        old_var_name (str, optional): Name of the old species gene symbol column. Defaults to "MouseGeneSymbol".
        new_var_name (str, optional): Name of the new species gene symbol column. Defaults to "HumanGeneSymbol".
        target_genes (list, optional): List of target genes to subset to. Defaults to None.
        var_names_column (str, optional): Name of the column in adata.var that contains the gene symbols. Defaults to None.

    Returns:
        anndata.AnnData: AnnData object with species converted.
    """
    # Create a copy of the AnnData object
    adata = adata.copy()
    # If var_names_column is provided, set adata.var_names accordingly
    if var_names_column:
        adata.var_names = adata.var[var_names_column]
    # Set the old species gene symbol column
    adata.var[old_var_name] = adata.var_names
    # Convert the gene symbols to the new species using the one-to-one dictionary
    adata.var[new_var_name] = adata.var_names.map(one2one_dict)
    # Remove genes that were not converted
    adata = adata[:, adata.var[new_var_name].notna()]
    # Set the new gene symbols as the var_names
    adata.var_names = adata.var[new_var_name]
    # Create a new AnnData object
    new_adata = adata
    # If target genes are provided, subset the data to those genes
    if target_genes is not None:
        print("Warning: not all target genes were found in the dataset")
        genes_not_in_adata = set(target_genes).difference(set(adata.var_names))
        adata_df = adata.to_df()[list(set(adata.var_names).intersection(set(target_genes)))]
        for gene in genes_not_in_adata:
            adata_df[gene] = 0
        adata_df = adata_df[target_genes]
        new_adata = anndata.AnnData(adata_df)
        new_adata.obs = adata.obs
        new_adata.var_names = adata_df.columns
        new_adata.uns = adata.uns
    return new_adata

def convert_species_one2one_df(df, one2one_dict, 
                              old_var_name="MouseGeneSymbol", 
                              new_var_name="HumanGeneSymbol", 
                              target_genes=None) -> pd.DataFrame:
    """
    Convert the species of a DataFrame using a one-to-one dictionary.

    Args:
        df (pd.DataFrame): DataFrame to convert.
        one2one_dict (dict): One-to-one dictionary of gene conversions.
        old_var_name (str, optional): Name of the old species gene symbol column. Defaults to "MouseGeneSymbol".
        new_var_name (str, optional): Name of the new species gene symbol column. Defaults to "HumanGeneSymbol".
        target_genes (list, optional): List of target genes to subset to. Defaults to None.

    Returns:
        pd.DataFrame: DataFrame with species converted.
    """
    # Create a copy of the DataFrame
    df = df.copy()
    # Convert the gene symbols to the new species using the one-to-one dictionary
    translated_series = df.columns.map(one2one_dict)
    df.columns = translated_series
    # Remove genes that were not converted
    df = df.loc[:, df.columns.notna()]
    # If target genes are provided, subset the data to those genes
    if target_genes is not None:
        genes_not_in_df = set(target_genes).difference(set(df.columns))
        for gene in genes_not_in_df:
            df[gene] = 0
        df = df[target_genes]
    
    return df

def full_convert_species_adata(adata,
                               homolog_table_path,
                               target_geneset=None,
                               old_var_name="MouseGeneSymbol",
                               new_var_name="HumanGeneSymbol",
                               var_names_column=None,
                               taxon_we_want=9606,
                               compatibility_taxon=10090) -> anndata.AnnData:
    """
    Convert species of adata to a new species using a homolog table
    Args:
        adata: anndata object
        homolog_table_path: path to homolog table
        target_geneset: geneset to convert/subset to
        old_var_name: name of the old species gene symbol column
        new_var_name: name of the new species gene symbol column
        var_names_column: name of the column in adata.var that contains the gene symbols
        taxon_we_want: taxon id of the species we want to convert to
        compatibility_taxon: taxon id of the species we are converting from
    Returns:
        anndata object with species converted
    """
    homolog_table = pd.read_csv(homolog_table_path, sep="\t")
    id2symbols, _, _ = get_id2symbols_dict(homolog_table)
    compatible_genes, reverse_translate_dict = get_compatible_genes(id2symbols)
    if taxon_we_want == 9606 and compatibility_taxon == 10090:
        one2one_dict = get_one_to_one_only(reverse_translate_dict, compatible_genes)
    else:
        one2one_dict = get_one_to_one_only(compatible_genes, reverse_translate_dict)
    adata = convert_species_one2one_adata(adata, one2one_dict, target_genes = target_geneset)
    return adata

def full_convert_species_series(input_gene_series,
                                homolog_table_path,
                                old_var_name="MouseGeneSymbol",
                                new_var_name="HumanGeneSymbol",
                                taxon_we_want=9606,
                                compatibility_taxon=10090
                                ) -> list:
    """
    Convert species of a pandas series of genes to a new species using a homolog table
    Args:
        input_gene_series: pandas series with list of genes to convert
        homolog_table_path: path to homolog table
        old_var_name: name of the old species gene symbol column
        new_var_name: name of the new species gene symbol column
        taxon_we_want: taxon id of the species we want to convert to
        compatibility_taxon: taxon id of the species we are converting from
    Returns:
        list of genes converted
    """
    homolog_table = pd.read_csv(homolog_table_path, sep="\t")
    id2symbols, _, _ = get_id2symbols_dict(homolog_table)
    compatible_genes, reverse_translate_dict = get_compatible_genes(id2symbols)
    if taxon_we_want == 9606 and compatibility_taxon == 10090:
        one2one_dict = get_one_to_one_only(reverse_translate_dict, compatible_genes)
    else:
        one2one_dict = get_one_to_one_only(compatible_genes, reverse_translate_dict)
    one2one_default_dict = defaultdict(lambda: pd.NA, one2one_dict)
    converted_genes = input_gene_series.map(one2one_default_dict)
    # remove genes that were not converted
    return converted_genes[converted_genes.notna()]