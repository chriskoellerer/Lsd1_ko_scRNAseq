########## miscellaneous QoL helpers

# get DR & 2D embeddings in Seurat format
get_embeddings_for_seurat <- function(file_path) {
  read.csv(file_path,
           row.names = 1) %>%
    setNames(seq_len(ncol(.))) %>%
    as.matrix()
}

# run quickndirty seurat flow
run_standard_seurat <- function(seurat_obj, npcs = 30) {
  seurat_obj %>%
    NormalizeData() %>%
    FindVariableFeatures() %>%
    ScaleData() %>%
    RunPCA(npcs = npcs) %>%
    FindNeighbors(dims = 1:npcs) %>%
    FindClusters() %>%
    RunUMAP(dims = 1:npcs)
}

# tag outlying data points by man absolute deviation (MAD)
mad_tagging <- function(object, category, nmads){
  data <- object@meta.data %>% pull(category)
  median_val <- median(data, na.rm = TRUE)
  mad_val <- mad(data, na.rm = TRUE)
  is_outlier <- (data < (median_val - nmads * mad_val)) | (data > (median_val + nmads * mad_val))
}


# save all the building blocks of a seurat object (ie in case of switching over to scanpy)
save_seurat_comps <- function(object, save_path){
  Matrix::writeMM(
    object[["RNA"]]@layers$counts,
    file.path(save_path, "matrix.mtx")
  )
  write.table(
    as.data.frame(rownames(object)),
    file.path(save_path, "features.tsv"),
    col.names = FALSE,
    row.names = FALSE,
    quote = FALSE)
  write.table(
    as.data.frame(colnames(object)),
    file.path(save_path, "barcodes.tsv"),
    col.names = FALSE,
    row.names = FALSE,
    quote = FALSE)
  write.csv(
    object@meta.data,
    file.path(save_path, "metadata.csv")
  )
}
