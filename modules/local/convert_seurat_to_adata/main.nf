process CONVERT_SEURAT_TO_ADATA{
    label "process_low"
    label "process_high_memory"
    conda "${moduleDir}/environment.yml"

    input:
    val(original_input_file)

    output:
    path("*.h5ad"), emit: converted_h5ad_file

    script:
    """
    Rscript \
    '${moduleDir}/convert_seurat_adata_sceasy.R' \
    '${original_input_file}'
    """
}