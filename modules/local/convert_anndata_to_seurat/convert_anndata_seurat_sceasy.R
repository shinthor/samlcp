library(sceasy)
library(reticulate)
library(Seurat)
args <- commandArgs(trailingOnly=TRUE)
input_file_path <- args[1]
# Get the input file name minus the extension
input_file_name <- tools::file_path_sans_ext(basename(input_file_path))

sceasy::convertFormat(
        input_file_path,
        assay="RNA",
        from="anndata",
        to="seurat",
        main_layer="counts",
        outFile=paste0(input_file_name, ".rds"))