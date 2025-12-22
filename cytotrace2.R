library(CytoTRACE2) 
library(argparse)
library(Seurat)
library(ggplot2)
library(patchwork)
library(readr)
library(tibble)

parser <- ArgumentParser(description = 'Integrate multiple scRNA-seq datasets into one Seurat object')
parser$add_argument("--rds", type="character", help="File path to joined_integrated_seurat_object.rds")
args <- parser$parse_args()

data <- readRDS(args$rds)

if (length(Layers(data[["RNA"]], search = "counts")) > 1) {
    cat("Joining counts layers...\n")
    data <- JoinLayers(data, assay = "RNA", layers = "counts")
}

batches <- unique(data$orig.ident)
for (batch in batches) {
    cat("Processing batch:", batch, "\n")
    
    # Subset to this batch
    batch_cells <- colnames(data)[data$orig.ident == batch]
    batch_data <- data[, batch_cells]

    # Extract expression data
    expression_data <- GetAssayData(batch_data, assay = "RNA", layer = "counts")
    
    # Run CytoTRACE2
    cytotrace2_result <- cytotrace2(expression_data, species = "human")
    
    
    # Prepare annotation
    annotation <- data.frame(Phenotype = batch_data$singleR.labels)
    
    # Write stats
    cytotrace2_result_with_rownames <- merge(cytotrace2_result, annotation, by = 0, all = TRUE)
    names(cytotrace2_result_with_rownames)[1] <- "barcode"
    write_csv(cytotrace2_result_with_rownames, paste0("cytotrace_stats_", batch, ".csv"))

    # Generate plots
    plots <- plotData(
        cytotrace2_result = cytotrace2_result,
        annotation = annotation,
        expression_data = expression_data
    )
    
    # Calculate dimensions
    n_plots <- length(plots)
    ncol <- 2
    nrow <- ceiling(n_plots / ncol)
    plot_width <- 8
    plot_height <- 6
    total_width <- plot_width * ncol
    total_height <- plot_height * nrow
    
    # Combine and save
    combined_plot <- wrap_plots(plots, ncol = ncol)
    ggsave(paste0("cytotrace_plots_", batch, ".pdf"), 
           combined_plot, 
           width = total_width, 
           height = total_height, 
           dpi = 300)
    
    cat("Saved: cytotrace_plots_", batch, ".pdf\n")
}

