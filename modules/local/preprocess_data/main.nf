process PREPROCESS_DATA {
    label "process_high"
    label "process_high_memory"
    conda "${moduleDir}/environment.yml"

    input:
    path(original_input_file)
    val(threshold_combinations)
    val(uns_name)
    val(sample_taxon)
    val(column2)
    val(homolog_table_path)
    val(cell_cycle_genes_path)
    val(use_raw)
    val(var_colname)
    val(layer_name)

    output:
    path("*_processed_data.h5ad"), emit: processed_data_file

    script:
    """
    python \
    '${workflow.projectDir}/bin/preprocess_data.py' \
    --input-file='${original_input_file}' \
    --threshold-combinations='${threshold_combinations}' \
    --output-file='${uns_name}_processed_data.h5ad' \
    --name-to-add='${uns_name}' \
    --sample-taxon='${sample_taxon}' \
    --column2='${column2}' \
    --homolog-table-path='${homolog_table_path}' \
    --cell-cycle-genes-path='${cell_cycle_genes_path}' \
    --use-raw='${use_raw}' \
    --var-colname='${var_colname}' \
    --layer-name='${layer_name}'
    """
}