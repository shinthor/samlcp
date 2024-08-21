args <- commandArgs(trailingOnly = TRUE)
input_file_path <- args[1]
# Get the input file name minus the extension
input_file_name <- tools::file_path_sans_ext(basename(input_file_path))

# Set conda env path to the one provided if provided and not NULL/empty string
if (length(args) > 1 && args[2] != "" && !is.null(args[2])) {
  reticulate::use_condaenv(args[2])
}

library(sceasy)
library(reticulate)
library(Seurat)

sceasy::convertFormat(
        input_file_path,
        assay="RNA",
        from="anndata",
        to="seurat",
        main_layer="counts",
        outFile=paste0(input_file_name, ".rds"))