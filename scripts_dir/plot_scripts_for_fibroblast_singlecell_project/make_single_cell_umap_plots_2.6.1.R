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

cell_type_mapping <- c(
    "Fibroblasts:breast" = "Superficial APCDD1",
    "iPS_cells:fibroblasts" = "Superficial APCDD1",
    "iPS_cells:PDB_fibroblasts" = "Superficial APCDD1",
    "iPS_cells:skin_fibroblast" = "Superficial APCDD1",
    "iPS_cells:adipose_stem_cells" = "Superficial APCDD1",
    "Chondrocytes:MSC-derived" = "Reticular_PI16",
    "Smooth_muscle_cells:bronchial" = "Reticular_PI16",
    "Smooth_muscle_cells:bronchial:vit_D" = "Reticular_PI16",
    "Tissue_stem_cells:iliac_MSC" = "Adipogenic_FABP4",
    "Fibroblasts:foreskin" = "Perivascular",
    "Smooth_muscle_cells:vascular:IL-17" = "Perivascular",
    "Osteoblasts" = "Hair_follicle_like",
    "Osteoblasts:BMP2" = "Hair_follicle_like",
    "iPS_cells:CRL2097_foreskin" = "Schwann-like",
    "MSC" = "Mitotic",
    "Tissue_stem_cells:adipose-derived_MSC_AM3" = "Mitotic",
    "Tissue_stem_cells:BM_MSC:BMP2" = "Stressed/Contractile",
    "Tissue_stem_cells:BM_MSC" = "Stressed/Contractile",
    "Tissue_stem_cells:BM_MSC:TGFb3" = "Stressed/Contractile"
)

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
            ggsave(file.path(outdir, paste0(sanitize_label_for_filename(cm), "_umap.pdf")),
                plot,
                width = 5,
                height = 5,
                limitsize = FALSE
            )
        }
    })
}
#joined_integrated_data <- readRDS("rds_files/joined_integrated_seurat_object_with_group_data_just_fibroblasts_pericytes_reclustered.rds")

parser <- ArgumentParser(description = 'Integrate multiple scRNA-seq datasets into one Seurat cbject')
parser$add_argument("--outdir", help="directory to write to", required = TRUE)
parser$add_argument("--diagnosis", help="directory to write to", required = TRUE)
args <- parser$parse_args()
print(args)

joined_integrated_data <- readRDS("rds_files/joined_integrated_seurat_object_with_group_data.rds")


wanted_groups = c("fibroblasts", "pericyte")
just_subset <- subset(joined_integrated_data, cell_group %in% wanted_groups)
just_subset <- subset(just_subset, diagnosis == args$diagnosis)
just_subset$new_cat_labels <- unname(cell_type_mapping[just_subset$singleR.labels_fine])

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

joined_integrated_data$singleR.labels_fine <- sanitize_label_for_filename(joined_integrated_data$singleR.labels_fine)
make_plots(integrated_data, "new_cat_labels", file.path(args$outdir, args$diagnosis, "new_cat_labels"))
