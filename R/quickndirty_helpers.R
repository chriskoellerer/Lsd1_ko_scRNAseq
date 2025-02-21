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

# for a given pair of idents, run FindMarkers()
find_markers_for_pair <- function(
    pair,
    seurat_obj,
    gene_of_interest = NULL,
    min_pct = 0.25,
    min_diff_pct = 0.1,
    out_dir = "marker_outputs"
) {
  # make sure dir exists
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  # extract idents
  ident1 <- pair[1]
  ident2 <- pair[2]
  
  # Seurat::FindMarkers()
  markers <- FindMarkers(
    object = seurat_obj,
    ident.1 = ident1,
    ident.2 = ident2,
    min.pct = min_pct,
    min.diff.pct = min_diff_pct
  )
  
  # export top100 markers
  comp_label <- paste0(gsub(" ", "_", ident1), "_vs_", gsub(" ", "_", ident2))
  csv_file <- file.path(out_dir, paste0("Top100Markers_", comp_label, ".csv"))
  
  top100 <- head(markers[order(markers$p_val_adj, decreasing = FALSE), ], 100)
  write.csv(top100, file = csv_file)
  
  # if gene of interest present, export as well
  if (!is.null(gene_of_interest) && gene_of_interest %in% rownames(markers)) {
    fav_file <- file.path(out_dir, paste0("FavoriteGene_", gene_of_interest, "_", comp_label, ".txt"))
    write.table(
      markers[gene_of_interest, , drop = FALSE],
      file = fav_file,
      sep = "\t",
      quote = FALSE,
      col.names = NA
    )
  }
  
  return(markers)
}