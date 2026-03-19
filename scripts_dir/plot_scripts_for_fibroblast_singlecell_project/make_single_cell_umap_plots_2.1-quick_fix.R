library(argparse)
library(Seurat)
library(ggplot2)
library(tidyverse)
library(pheatmap)
library(purrr)
library(patchwork)

sanitize_label_for_filename <- function(label) {
    gsub("[/\\:*?\"<>|]", "_", label)
}

make_plots <- function(integrated_data, cell_metadata, outdir, plot_order = NULL){
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
    write.csv(table(integrated_data[[cell_metadata]]), file.path(outdir, "cellcounts.csv"))

    umap_plot <- DimPlot(
        integrated_data,
        reduction = 'umap.cca',
        group.by = cell_metadata,
    )
    umap_plot_path <- file.path(outdir, paste0("all_", cell_metadata, "_umap.pdf"))
    ggsave(umap_plot_path, umap_plot, width = 10, height = 5, dpi = 300)


    umap_coords <- as.data.frame(Embeddings(integrated_data, "umap.cca"))
    colnames(umap_coords) <- c("UMAP1", "UMAP2")
    umap_coords$cell_metadata <- integrated_data@meta.data[[cell_metadata]]

    if (is.null(plot_order)){
        cell_metadata_uniq <- unique(umap_coords$cell_metadata)
    }else{
        cell_metadata_uniq <- plot_order
    }

    lapply(cell_metadata_uniq, function(cm) {
        print(paste(cell_metadata, cm))
        umap_coords$highlight <- ifelse(umap_coords$cell_metadata == cm, "Selected", "Other")
        umap_coords <- umap_coords[order(umap_coords$highlight, decreasing = FALSE), ]
        plot <- ggplot(umap_coords, aes(UMAP1, UMAP2, color = highlight)) +
            geom_point(size = 0.3) +
            scale_color_manual(values = c("Other" = "lightgrey", "Selected" = "red")) +
            ggtitle(cm) +
            theme_minimal() +
            theme(legend.position = "none")

        if (cm == "n/a"){
            ggsave(file.path(outdir, "na_umap.pdf"),
                plot,
                width = 5,
                height = 5,
                limitsize = FALSE
            )
        }else{
            ggsave(file.path(outdir, paste0(cm, "_umap.pdf")),
                plot,
                width = 5,
                height = 5,
                limitsize = FALSE
            )
        }
    })
}

parser <- ArgumentParser(description = 'Integrate multiple scRNA-seq datasets into one Seurat object')
parser$add_argument("rds", help="File path to rds of choice")
parser$add_argument("--outdir", help="directory to write to", required = TRUE)
args <- parser$parse_args()

integrated_data <- readRDS(args$rds)
integrated_data$singleR.labels_main <- sanitize_label_for_filename(integrated_data$singleR.labels_main)
integrated_data$singleR.labels_fine <- sanitize_label_for_filename(integrated_data$singleR.labels_fine)

#integrated_data$singleR.labels_fine_plus_diagnosis <- paste0(integrated_data$singleR.labels_fine, "_", integrated_data$diagnosis)

#make_plots(integrated_data, "sample_id", file.path(args$outdir, "sample_id"))
#make_plots(integrated_data, "diagnosis", file.path(args$outdir, "diagnosis"))
#make_plots(integrated_data, "cca_clusters", file.path(args$outdir, "cca_clusters"))
#make_plots(integrated_data, "dataset", file.path(args$outdir, "dataset"))
#make_plots(integrated_data, "singleR.labels_fine", file.path(args$outdir, "singleR.labels_fine"))
#make_plots(integrated_data, "cell_group", file.path(args$outdir, "cell_group"))

for (diag in unique(integrated_data$diagnosis)) {
  col_name <- paste0("singleR_", diag)
  just_sub <- subset(integrated_data, diagnosis == diag)
  make_plots(just_sub, "singleR.labels_fine", file.path(args$outdir, col_name))
}
