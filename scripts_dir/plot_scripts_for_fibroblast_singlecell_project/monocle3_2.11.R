library(monocle3)
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(argparse)
library(patchwork)
library(magrittr)
library(lubridate)

mapping <- read.csv("cell_type_mapping.csv")

recluster <- function(object){
    cell_counts <- table(object$sample_id)
    keep_samples <- names(cell_counts[cell_counts >= 100])
    rejected_samples <- names(cell_counts[cell_counts < 100])
    print("rejected samples:")
    print(rejected_samples)
    object <- subset(object, sample_id %in% keep_samples)

    object[["RNA"]] <- split(object[["RNA"]], f = object$sample_id)
    object <- NormalizeData(object)
    object <- FindVariableFeatures(object)
    object <- ScaleData(object)
    object <- RunPCA(object)

    object <- FindNeighbors(object, dims = 1:30, reduction = "pca")
    object <- FindClusters(object, resolution = 0.5, cluster.name = "unintegrated_clusters")
    object <- RunUMAP(object, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")

    integrated_data <- IntegrateLayers(
      object = object, method = CCAIntegration,
      orig.reduction = "pca", new.reduction = "integrated.cca",
      verbose = FALSE
    )

    integrated_data <- FindNeighbors(integrated_data, reduction = "integrated.cca", dims = 1:30)
    integrated_data <- FindClusters(integrated_data, cluster.name = "cca_clusters")

    integrated_data <- RunUMAP(integrated_data, reduction = "integrated.cca", dims = 1:30, reduction.name = "umap.cca")
    integrated_data <- JoinLayers(integrated_data, assay = "RNA")
    return(integrated_data)
}

make_multiplots <- function(cds, cell_metadata, outdir, filenm, level_order = NULL){
    umap_coords <- as.data.frame(reducedDims(cds)[["UMAP"]])
    colnames(umap_coords) <- c("UMAP1", "UMAP2")
    umap_coords$cell_metadata <- colData(cds)[[cell_metadata]]
    umap_coords$cluster <- cds@clusters$UMAP$clusters

    if (is.null(level_order)){
        cell_metadata_uniq <- unique(umap_coords$cell_metadata)
    }else{
        cell_metadata_uniq <- level_order
    }
    highlight_plots <- lapply(cell_metadata_uniq, function(cm) {
      umap_coords$highlight <- ifelse(umap_coords$cell_metadata == cm, "Selected", "Other")
      # Plot "Other" first so highlighted cells are drawn on top
      umap_coords <- umap_coords[order(umap_coords$highlight, decreasing = FALSE), ]
      ggplot(umap_coords, aes(UMAP1, UMAP2, color = highlight)) +
        geom_point(size = 0.3) +
        scale_color_manual(values = c("Other" = "lightgrey", "Selected" = "red")) +
        ggtitle(cm) +
        theme_minimal() +
        theme(legend.position = "none")
    })

    highlight_plot <- wrap_plots(highlight_plots, ncol = min(c(length(cell_metadata_uniq), 4)))

    ggsave(file.path(outdir, filenm), 
        highlight_plot,
        width = 5 * min(c(length(cell_metadata_uniq), 4)),
        height = 5 * ceiling(length(cell_metadata_uniq) / 4)
    )
}

run_monocle3 <- function(seurat_object, outdir){
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

    cds <- as.cell_data_set(seurat_object)
    reducedDims(cds)[["UMAP"]] <- Embeddings(seurat_object, "umap.cca")
    cds <- cluster_cells(cds)
    plot1 <- plot_cells(cds, color_cells_by = "m_labels", show_trajectory_graph = FALSE, label_cell_groups = FALSE)
    plot2 <- plot_cells(cds, color_cells_by = "new_cat_labels", show_trajectory_graph = FALSE, label_cell_groups = FALSE)
    plot3 <- plot_cells(cds, color_cells_by = "CytoTRACE2_Potency", show_trajectory_graph = FALSE, label_cell_groups = FALSE)
    plot <- wrap_plots(plot1, plot2, plot3, ncol = 3)

    filenm <- file.path(outdir, "labels_and_potency.pdf")
    ggsave(filenm, plot, width = 24, height = 8)

    make_multiplots(cds, "m_labels", outdir, "celltype_highlights.pdf")
    make_multiplots(cds, "new_cat_labels", outdir, "new_cat_labels_highlights.pdf")
    make_multiplots(cds, "CytoTRACE2_Potency", outdir, "CytoTRACE2_Potency_highlights.pdf", c("Multipotent", "Oligopotent", "Unipotent", "Differentiated"))

    cds <- learn_graph(cds)
    plot1 <- plot_cells(
        cds,
        color_cells_by = "m_labels",
        label_groups_by_cluster = FALSE,
        label_leaves = FALSE,
        label_branch_points = FALSE,
        label_cell_groups = FALSE
    )

    plot2 <- plot_cells(
        cds,
        color_cells_by = "new_cat_labels",
        label_groups_by_cluster = FALSE,
        label_leaves = FALSE,
        label_branch_points = FALSE,
        label_cell_groups = FALSE
    )

    plot3 <- plot_cells(
        cds,
        color_cells_by = "CytoTRACE2_Potency",
        label_groups_by_cluster = FALSE,
        label_leaves = FALSE,
        label_branch_points = FALSE,
        label_cell_groups = FALSE
    )

    plot <- wrap_plots(plot1, plot2, plot3, ncol = 3)
    filenm <- file.path(outdir, "principal_graphs.pdf")
    ggsave(filenm, plot, width = 15, height = 5)

    max_ct2 <- which.max(unlist(FetchData(seurat_object, "CytoTRACE2_Score")))
    max_ct2 <- colnames(seurat_object)[max_ct2]
    cds <- order_cells(cds, root_cells = max_ct2)
    plot1 <- plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE, label_branch_points = FALSE)
    plot2 <- plot_cells(cds, color_cells_by = "CytoTRACE2_Score", label_cell_groups = FALSE, label_leaves = FALSE, label_branch_points = FALSE)

    plot <- wrap_plots(plot1, plot2)
    filenm <- file.path(outdir, "principal_graph_cytotrace_score.pdf")
    ggsave(filenm, plot, width = 10, height = 5)
}

date_and_time <- format(now(), "%Y%m%d_%H%M%S")
parser <- ArgumentParser(description = 'Integrate multiple scRNA-seq datasets into one Seurat object')
parser$add_argument("--outdir", help="directory to write to", required = TRUE)
args <- parser$parse_args()
print(args)


just_subset <- readRDS("rds_files/cytotrace_seurat.rds")
wanted_groups = c("fibroblasts", "pericyte")
just_subset <- subset(just_subset, cell_group %in% wanted_groups)
just_subset <- subset(just_subset, !is.na(CytoTRACE2_Potency))

just_subset$m_labels <- mapping$m_labels[match(just_subset$singleR.labels_fine, mapping$singleR.labels_fine)]
just_subset$new_cat_labels <- mapping$new_cat_labels[match(just_subset$singleR.labels_fine, mapping$singleR.labels_fine)]

for (diag in unique(just_subset$diagnosis)){
    just_diag <- subset(just_subset, diagnosis == diag)
    just_diag <- recluster(just_diag)
    outdir <- file.path(args$outdir, diag)
    run_monocle3(just_diag, outdir)
}
