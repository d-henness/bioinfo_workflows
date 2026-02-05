library(monocle3)
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(argparse)
library(patchwork)
library(magrittr)
library(lubridate)

fibroblast_celltypes <- c(
    "Chondrocytes:MSC-derived",
    "Fibroblasts:breast",
    "iPS_cells:adipose_stem_cells",
    "iPS_cells:CRL2097_foreskin",
    "iPS_cells:fibroblasts",
    "iPS_cells:PDB_fibroblasts",
    "Osteoblasts",
    "Smooth_muscle_cells:bronchial",
    "Smooth_muscle_cells:bronchial:vit_D",
    "Tissue_stem_cells:BM_MSC:TGFb3",
    "Tissue_stem_cells:BM_MSC"
)

make_multiplots <- function(cds, cell_metadata, outdir, filenm, order = NULL){
    umap_coords <- as.data.frame(reducedDims(cds)[["UMAP"]])
    colnames(umap_coords) <- c("UMAP1", "UMAP2")
    umap_coords$cell_metadata <- colData(cds)[[cell_metadata]]
    umap_coords$cluster <- cds@clusters$UMAP$clusters

    if (is.null(order)){
        cell_metadata_uniq <- unique(umap_coords$cell_metadata)
    }else{
        cell_metadata_uniq <- order
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
    plot1 <- plot_cells(cds, color_cells_by = "singleR.labels_fine", show_trajectory_graph = FALSE, label_cell_groups = FALSE)
    plot2 <- plot_cells(cds, color_cells_by = "CytoTRACE2_Potency", show_trajectory_graph = FALSE, label_cell_groups = FALSE)
    plot <- wrap_plots(plot1, plot2)

    filenm <- file.path(outdir, "labels_and_potency.pdf")
    ggsave(filenm, plot, width = 16, height = 8)

    make_multiplots(cds, "singleR.labels_fine", outdir, "celltype_highlights.pdf")
    make_multiplots(cds, "CytoTRACE2_Potency", outdir, "CytoTRACE2_Potency_highlights.pdf", c("Multipotent", "Oligopotent", "Unipotent", "Differentiated"))
    make_multiplots(cds, "diagnosis", outdir, "diagnosis_highlights.pdf")

    cds <- learn_graph(cds)
    plot1 <- plot_cells(
        cds,
        color_cells_by = "singleR.labels_fine",
        label_groups_by_cluster = FALSE,
        label_leaves = FALSE,
        label_branch_points = FALSE,
        label_cell_groups = FALSE
    )

    plot2 <- plot_cells(
        cds,
        color_cells_by = "CytoTRACE2_Potency",
        label_groups_by_cluster = FALSE,
        label_leaves = FALSE,
        label_branch_points = FALSE,
        label_cell_groups = FALSE
    )

    plot3 <- plot_cells(
        cds,
        color_cells_by = "diagnosis",
        label_groups_by_cluster = FALSE,
        label_leaves = FALSE,
        label_branch_points = FALSE,
        label_cell_groups = FALSE
    )
    plot <- wrap_plots(plot1, plot2, plot3)
    filenm <- file.path(outdir, "principal_graph.pdf")
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
parser$add_argument("joined_integrated_seurat_object", help="File paths to joined_integrated_seurat_object.rds")
args <- parser$parse_args()
print(args)


joined_integrated_data <- readRDS(args$joined_integrated_seurat_object)
joined_integrated_data$diagnosis <- trimws(joined_integrated_data$diagnosis)

just_fibs <- subset(joined_integrated_data, singleR.labels_fine %in% fibroblast_celltypes)
just_control <- subset(just_fibs, diagnosis == 'control')
just_SSC <- subset(just_fibs, diagnosis == 'SSC')

run_monocle3(just_fibs, paste0("monocle3_subset_fibs_using_origional_clustering_", date_and_time))
run_monocle3(just_SSC, paste0("monocle3_subset_fibs_subset_SSC_using_origional_clustering_", date_and_time))
run_monocle3(just_control, paste0("monocle3_subset_fibs_subset_control_using_origional_clustering_", date_and_time))

just_fibs <- FindNeighbors(just_fibs, reduction = "integrated.cca", dims = 1:30)
just_fibs <- FindClusters(just_fibs, cluster.name = "cca_clusters")
just_fibs <- RunUMAP(just_fibs, reduction = "integrated.cca", dims = 1:30, reduction.name = "umap.cca")
just_control <- subset(just_fibs, diagnosis == 'control')
just_SSC <- subset(just_fibs, diagnosis == 'SSC')

run_monocle3(just_fibs, paste0("monocle3_subset_fibs_and_reclustered_", date_and_time))
run_monocle3(just_SSC, paste0("monocle3_subset_fibs_and_reclustered_subset_SSC_", date_and_time))
run_monocle3(just_control, paste0("monocle3_subset_fibs_and_reclustered_subset_control", date_and_time))

just_control <- FindNeighbors(just_control, reduction = "integrated.cca", dims = 1:30)
just_control <- FindClusters(just_control, cluster.name = "cca_clusters")
just_control <- RunUMAP(just_control, reduction = "integrated.cca", dims = 1:30, reduction.name = "umap.cca")

just_SSC <- FindNeighbors(just_SSC, reduction = "integrated.cca", dims = 1:30)
just_SSC <- FindClusters(just_SSC, cluster.name = "cca_clusters")
just_SSC <- RunUMAP(just_SSC, reduction = "integrated.cca", dims = 1:30, reduction.name = "umap.cca")

run_monocle3(just_SSC, paste0("monocle3_subset_fibs_subset_SSC_and_reclustered_", date_and_time))
run_monocle3(just_control, paste0("monocle3_subset_fibs_subset_control_and_reclustered_", date_and_time))
