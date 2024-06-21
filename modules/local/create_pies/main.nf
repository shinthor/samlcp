process CREATE_PIE_CHARTS {
    label "process_low"
    conda "${moduleDir}/environment.yml"

    input:
    path(processed_data_file)
    val(uns_name)
    val(column1)
    val(column2)

    output:
    path("*.png")

    script:
    """
    Rscript \
    ${workflow.projectDir}/bin/create_pie_charts.R \
    $gene_combination \
    "${uns_name}" \
    "${column1}" \
    "${column2}" \
    "."
    """
}