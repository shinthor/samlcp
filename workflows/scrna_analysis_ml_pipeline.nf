/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-validation'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_scrna_analysis_ml_pipeline_pipeline'
include { PREPROCESS_ADATA        } from '../modules/local/preprocess_adata'
// include { ANALYZE_SINGLE_CATEGORY } from '../modules/local/analyze_single_category'
// include { COMBINE_CATEGORIES      } from '../modules/local/combine_categories'

include { GROUP_ANALYSIS_COMBINATIONS } from  '../modules/local/group_analysis_combinations'
include { CREATE_PIES             } from '../modules/local/create_pies'
include { CREATE_BARS             } from '../modules/local/create_bars'
include { CONVERT_SEURAT_TO_ANNDATA   } from '../modules/local/convert_seurat_to_anndata'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SCRNA_ANALYSIS_ML_PIPELINE {

    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    ch_folder_paths = ch_samplesheet.map { it[1] }
    print("Folder paths before: ${ch_folder_paths.collect()}")
    // Convert Seurat to AnnData for RDS files
    if (params.seurat_input) {
        CONVERT_SEURAT_TO_ANNDATA(ch_folder_paths)
        ch_folder_paths = CONVERT_SEURAT_TO_ANNDATA.out
    }
    print("Folder paths: ${ch_folder_paths.collect()}")

    homolog_table_path = Channel.of("${workflow.projectDir}/resources/HOM_MouseHumanSequence.rpt")
    cell_cycle_genes_path = Channel.of("${workflow.projectDir}/resources/regev_lab_cell_cycle_genes.txt")

    // RUN PREPROCESS_DATA
    PREPROCESS_ADATA(
        ch_folder_paths,
        params.threshold_combinations,
        params.uns_name,
        params.sample_taxon,
        params.column2,
        homolog_table_path,
        cell_cycle_genes_path,
        params.use_raw,
        params.var_colname,
        params.layer_name,
        params.save_h5ad
        )

    // RUN GROUP_ANALYSIS_COMBINATIONS
    GROUP_ANALYSIS_COMBINATIONS(
        PREPROCESS_ADATA.out.processed_tsv.flatten(), 
        params.uns_name,
        params.column1,
        params.column2,
        params.sample_column
        )

    // RUN CREATE_PIES
    CREATE_PIES(
        GROUP_ANALYSIS_COMBINATIONS.out.gene_combinations.flatten(),
        params.uns_name,
        params.column1,
        params.column2
        )
    
    CREATE_BARS(
        GROUP_ANALYSIS_COMBINATIONS.out.gene_combinations_with_sample.flatten(),
        params.uns_name,
        params.column1,
        params.column2
        )


    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_pipeline_software_mqc_versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))

    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
