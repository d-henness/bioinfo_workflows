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
parser$add_argument("--diagnosis", help="directory to write to", required = TRUE)
args <- parser$parse_args()

joined_integrated_data <- readRDS(args$rds)
joined_integrated_data$singleR.labels_main <- sanitize_label_for_filename(joined_integrated_data$singleR.labels_main)
joined_integrated_data$singleR.labels_fine <- sanitize_label_for_filename(joined_integrated_data$singleR.labels_fine)


wanted_groups = c("fibroblasts", "pericyte")
just_subset <- subset(joined_integrated_data, cell_group %in% wanted_groups)
just_subset <- subset(just_subset, diagnosis == args$diagnosis)

rm(joined_integrated_data)
gc()

cell_counts <- table(just_subset$sample_id)
keep_samples <- names(cell_counts[cell_counts >= 100])
rejected_samples <- names(cell_counts[cell_counts < 100])
print("rejected samples:")
print(rejected_samples)
just_subset <- subset(just_subset, sample_id %in% keep_samples)

just_subset <- JoinLayers(just_subset, assay = "RNA")
just_subset[["RNA"]] <- split(just_subset[["RNA"]], f = just_subset$sample_id)

just_subset <- NormalizeData(just_subset)
just_subset <- FindVariableFeatures(just_subset)
just_subset <- ScaleData(just_subset)
just_subset <- RunPCA(just_subset)

just_subset <- FindNeighbors(just_subset, dims = 1:30, reduction = "pca")

integrated_data <- IntegrateLayers(
  object = just_subset, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.cca",
  verbose = FALSE
)
# get rid of big objects
rm(just_subset)
gc()

integrated_data <- FindNeighbors(integrated_data, reduction = "integrated.cca", dims = 1:30)
integrated_data <- FindClusters(integrated_data, cluster.name = "cca_clusters")

integrated_data <- RunUMAP(integrated_data, reduction = "integrated.cca", dims = 1:30, reduction.name = "umap.cca")
integrated_data <- JoinLayers(integrated_data, assay = "RNA")

make_plots(integrated_data, "singleR.labels_fine", file.path(args$outdir, args$diagnosis, "singleR.labels_fine"))
