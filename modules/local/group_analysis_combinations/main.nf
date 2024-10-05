process GROUP_ANALYSIS_COMBINATIONS {
    label "process_medium"
    maxRetries 4
    conda "${moduleDir}/environment.yml"

    input:
    path(processed_data_file)
    val(uns_name)
    val(column1)
    val(column2)
    val(sample_column)

    output:
    path("*_gene_combinations.tsv"), emit: gene_combinations
    path("*_gene_combinations_with_sample.tsv"), emit: gene_combinations_with_sample

    script:
    """
    python \
    ${workflow.projectDir}/bin/group_analysis_combinations.py \
    --input_path=$processed_data_file \
    --output_path_base="." \
    --uns_name="${uns_name}" \
    --level1_category="${column1}" \
    --level2_category="${column2}" \
    --sample_category_column_name="${sample_column}"
    """
}