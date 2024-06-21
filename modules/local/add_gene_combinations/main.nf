process ADD_GENE_COMBINATIONS {
    label "process_medium"
    label "process_high_memory"
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
    ${workflow.projectDir}/bin/add_gene_combinations.py \
    --input_path=$processed_data_file \
    --output_path_base="." \
    --uns_name="${params.uns_name}" \
    --level1_category="${params.column1}" \
    --level2_category="${params.column2}"
    """
}