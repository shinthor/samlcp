library(sceasy)
library(reticulate)
library(Seurat)
args <- commandArgs(trailingOnly=TRUE)
input_file_path <- args[1]
# Get the input file name minus the extension
input_file_name <- tools::file_path_sans_ext(basename(input_file_path))

seurat_object <- readRDS(input_file_path)
object_ver <- Version(seurat_object)

# Fix v5 object to v4 compatibility
# sceasy already checks for version lower than 3.0 and runs UpdateSeuratObject
if (object_ver$major > 4) {
    seurat_object[["RNA3"]] <- as(object = seurat_object[["RNA"]], Class = "Assay")
    DefaultAssay(seurat_object) <- "RNA3"
    seurat_object[["RNA"]] <- NULL
    seurat_object <- RenameAssays(object = seurat_object, RNA3 = 'RNA')
    sceasy::convertFormat(
        seurat_object,
        assay="RNA",
        from="seurat",
        to="anndata",
        main_layer="counts",
        drop_single_values=FALSE,
        outFile=paste0(input_file_name, ".h5ad"))
} else {
    sceasy::convertFormat(
        seurat_object,
        assay="RNA",
        from="seurat",
        to="anndata",
        main_layer="counts",
        drop_single_values=FALSE,
        outFile=paste0(input_file_name, ".h5ad"))
}