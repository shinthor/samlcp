process CONVERT_ANNDATA_TO_SEURAT{
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
    '${moduleDir}/convert_anndata_seurat_sceasy.R' \
    '${original_input_file}'
    """
}