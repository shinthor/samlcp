process GROUP_ANALYSIS_COMBINATIONS {
    label "process_medium"
    maxRetries 4
    conda "${moduleDir}/environment.yml"

    input:
    path(processed_data_file)
    val(uns_name)
    val(column1)
    val(column2)

    output:
    path("*_gene_combinations.tsv")


    script:
    """
    python \
    ${workflow.projectDir}/bin/group_analysis_combinations.py \
    --input_path=$processed_data_file \
    --output_path_base="." \
    --uns_name="${uns_name}" \
    --level1_category="${column1}" \
    --level2_category="${column2}"
    """
}