process CREATE_BARS{
    label "process_low"
    conda "${moduleDir}/environment.yml"

    input:
    path(processed_data_file)
    val(uns_name)
    val(column1)
    val(column2)

    output:
    path("*.{png,pdf}")

    script:
    """
    Rscript \
    ${workflow.projectDir}/bin/create_bar_charts.R \
    $processed_data_file \
    "${uns_name}" \
    "${column1}" \
    "${column2}" \
    "."
    """
}