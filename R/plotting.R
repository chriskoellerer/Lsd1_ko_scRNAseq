
# general figure params 
# ~4:3
std_width = 10.54
std_height = 8.10

# call & save scCustomize::Dimplot_scCustom()
save_dimplot <- function(object, 
                         reduction = "umap", 
                         category = c("seurat_clusters"), 
                         num_columns = 2, 
                         title = NULL, 
                         no_legend = FALSE,
                         folder = NULL,
                         output_formats = c("pdf", "png")) {
  
  # generate plot
  plot <- scCustomize::DimPlot_scCustom(
    seurat_object = object,
    reduction = reduction,
    group.by = category,
    num_columns = num_columns
  )
  
  # optional: remove legend for publication
  if (no_legend) {
    plot <- plot +
      NoLegend()
  }
  
  # save formats
  for (format in output_formats) {
    
    if(length(category) > 1){
      filename <- paste0(reduction, "_", title, ".", format)
    } else {
      filename <- paste0(reduction, "_", category, ".", format) 
    }
    
    # wether to use results_folder
    if(is.null(folder)) {
      ggsave(filename, plot = plot, bg = "white", width = num_columns*std_width)
    } else {
      ggsave(file.path(folder, filename), plot = plot, bg = "white", width = num_columns*std_width)
    }
  }
  
  # return da plot
  return(plot)
  
}


# call & save scCustomize::FeaturePlot_scCustom()
save_featureplot_big <- function(object, 
                                 features, 
                                 reduction = "umap", 
                                 num_columns = 2, 
                                 title, 
                                 pt_size = NULL,
                                 order = FALSE,
                                 folder = NULL,
                                 no_legend = FALSE,
                                 output_formats = c("pdf", "png")) {
  
  # generate plot
  plot <- scCustomize::FeaturePlot_scCustom(
    seurat_object = object,
    features = features,
    reduction = reduction,
    pt.size = pt_size,
    order = order,
    num_columns = num_columns
  )
  
  # optional: remove legend for publication
  if (no_legend) {
    plot <- plot +
      NoLegend()
  }
  
  # save formats
  for (format in output_formats) {
    # filename
    filename <- paste0(reduction, "_", title, ".", format)
    # folder
    if(is.null(folder)) {
      ggsave(filename, plot = plot, bg = "white")
    } else {
      ggsave(file.path(folder, filename), plot = plot, bg = "white")
    }
  }
  # get da plot
  return(plot)
}
