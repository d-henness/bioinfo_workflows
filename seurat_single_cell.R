library(argparse)
library(Seurat)
library(ggplot2)
library(SingleR)
library(celldex)
library(tidyverse)
library(pheatmap)

fibroblast_celltypes <- c(
    "Chondrocytes:MSC-derived",
    "Fibroblasts:breast",
    "iPS_cells:adipose_stem_cells",
    "iPS_cells:fibroblasts",
    "iPS_cells:PDB_fibroblasts",
    "Neurons:Schwann_cell",
    "Osteoblasts",
    "Smooth_muscle_cells:bronchial",
    "Smooth_muscle_cells:bronchial:vit_D",
    "iPS_cells:CRL2097_foreskin",
    "Tissue_stem_cells:BM_MSC",
    "Tissue_stem_cells:BM_MSC:TGFb3"
)

sanitize_label_for_filename <- function(label) {
    gsub("[/\\:*?\"<>|]", "_", label)
}

# Create an ArgumentParser object
parser <- ArgumentParser(description = 'Integrate multiple scRNA-seq datasets into one Seurat object')
parser$add_argument("--datasets", nargs='+', type="character", help="File paths to datasets")
parser$add_argument("--min_cells", type="integer", help="Min cells per celltype for diff express", default = 50)

# Parse the command-line arguments
args <- parser$parse_args()

saved_objects_filename <- "scRNA_seurat_objects.rds"
if (file.exists(saved_objects_filename)) {
    seurat_objects <- readRDS(saved_objects_filename)
}else{
    # Read and process each dataset
    seurat_objects <- lapply(args$datasets, function(dataset_path) {
        # Load dataset
        print(dataset_path)
        data <- Read10X(data.dir = dataset_path)

        # Create Seurat object and filter
        print(basename(dataset_path))
        seurat_object <- CreateSeuratObject(counts = data, project = basename(dataset_path))
        seurat_object$mitoPercent <- PercentageFeatureSet(seurat_object, pattern = '^MT-')
        seurat_object.filtered <- subset(seurat_object, subset = nCount_RNA > 800 &
            nFeature_RNA > 500 &
            mitoPercent < 10)

        # Pre-processing: Normalization, variable features, scaling
        seurat_object.filtered <- NormalizeData(seurat_object.filtered) %>% FindVariableFeatures() %>% ScaleData()

        return(seurat_object.filtered)
    })
  saveRDS(seurat_objects, saved_objects_filename)
}

# Integration steps
# Find integration anchors
saved_integration_anchors_filename <- "scrna_integration_anchors.rds"
if (file.exists(saved_integration_anchors_filename)) {
    integration_anchors <- readRDS(saved_integration_anchors_filename)
}else{
    integration_anchors <- FindIntegrationAnchors(object.list = seurat_objects, dims = 1:30)
    saveRDS(integration_anchors, saved_integration_anchors_filename)
}
print("here")
# Integrate data

saved_integrated_data <- "integrated_seurat_object.rds"
if (file.exists(saved_integrated_data)){
    integrated_data <- readRDS(saved_integrated_data)
}else{
    integrated_data <- IntegrateData(anchorset = integration_anchors, dims = 1:30)

    # Post-integration steps (scaling, PCA, UMAP, clustering)
    integrated_data <- ScaleData(integrated_data) %>%
        RunPCA() %>%
        FindNeighbors(reduction = "pca", dims = 1:30) %>%
        FindClusters(resolution = 0.5) %>%
        RunUMAP(reduction = "pca", dims = 1:30) %>%
        RunTSNE(reduction = "pca", dims = 1:30)
    saveRDS(integrated_data, file = saved_integrated_data)
}

print(length(unique(integrated_data$seurat_clusters)))
print(unique(integrated_data$seurat_clusters))

# Visualization
plot <- DimPlot(integrated_data, reduction = "umap", label = TRUE, group.by = c("orig.ident", "seurat_clusters"))
ggsave("integrated_data.pdf", plot, width = 16, height = 8, dpi = 300)

ref <- celldex::HumanPrimaryCellAtlasData()

print(integrated_data[["RNA"]])

joined_integrated_data <- JoinLayers(integrated_data, assay = "RNA", layers = "data")
print(head(joined_integrated_data))

print(DefaultAssay(joined_integrated_data))
DefaultAssay(joined_integrated_data) <- "RNA"
print(DefaultAssay(joined_integrated_data))

counts <- GetAssayData(joined_integrated_data, assay = "RNA", layer = "data")

#pred <- SingleR(
#    test = counts,
#    ref = ref,
#    labels = ref$label.main
#)
#
#
#joined_integrated_data$singleR.labels_main <- pred$labels[match(rownames(joined_integrated_data@meta.data), rownames(pred))]
#print(paste0("Number of main labels:  ", length(unique(joined_integrated_data$singleR.labels_main))))
#write.csv(table(joined_integrated_data$singleR.labels_main), "cell_type_counts_main.csv")
#
#plot <- DimPlot(joined_integrated_data, reduction = 'umap', group.by = c('seurat_clusters', 'singleR.labels_main'))
#ggsave("joined_integrated_data_tissue_labels_main.pdf", plot, width = 16, height = 8, dpi = 300)
#
#plot <- plotScoreHeatmap(pred)
#ggsave("score_heatmap_main.pdf", plot, width = 16, height = 8, dpi = 300)
#
#plot <- plotDeltaDistribution(pred)
#ggsave("delta_distribution_main.pdf", plot, width = 16, height = 8, dpi = 300)

pred <- SingleR(
    test = counts,
    ref = ref,
    labels = ref$label.fine,
)

joined_integrated_data$singleR.labels_fine <- pred$labels[match(rownames(joined_integrated_data@meta.data), rownames(pred))]
print(paste0("Number of fine labels:  ", length(unique(joined_integrated_data$singleR.labels_fine))))
write.csv(table(joined_integrated_data$singleR.labels_fine), "cell_type_counts_fine.csv")

plot <- DimPlot(joined_integrated_data, reduction = 'umap', group.by = c('seurat_clusters', 'singleR.labels_fine'))
ggsave("joined_integrated_data_tissue_labels_fine.pdf", plot, width = 16, height = 8, dpi = 300)

# Make plots for each cell type

cell_types <- unique(joined_integrated_data$singleR.labels_fine)
umap_plot_dir <- "celltype_umap_plots"
tsne_plot_dir <- "celltype_tsne_plots"
dir.create(umap_plot_dir, showWarnings = FALSE)
dir.create(tsne_plot_dir, showWarnings = FALSE)

for (cell_type in cell_types){
    # Create a new column that highlights only this cell type
    joined_integrated_data$highlight <- ifelse(
        joined_integrated_data$singleR.labels_fine == cell_type,
        cell_type,
        "Other"
    )
    
    # Reorder so "Other" cells are plotted first (behind highlighted cells)
    joined_integrated_data$highlight <- factor(
        joined_integrated_data$highlight,
        levels = c("Other", cell_type)
    )
    
    # Create the plot
    umap_plot <- DimPlot(
        joined_integrated_data,
        reduction = 'umap',
        group.by = 'highlight',
        cols = c("Other" = "lightgrey", setNames("#E41A1C", cell_type)),
        order = cell_type  # Plot highlighted cells on top
    )

    tsne_plot <- DimPlot(
        joined_integrated_data,
        reduction = 'tsne',
        group.by = 'highlight',
        cols = c("Other" = "lightgrey", setNames("#E41A1C", cell_type)),
        order = cell_type  # Plot highlighted cells on top
    )

    cell_type_label <- sanitize_label_for_filename(cell_type)
    umap_plot_path <- file.path(umap_plot_dir, paste0(cell_type_label, "_vs_all_umap.pdf"))
    tsne_plot_path <- file.path(tsne_plot_dir, paste0(cell_type_label, "_vs_all_tsne.pdf"))
    ggsave(umap_plot_path, umap_plot, width = 16, height = 8, dpi = 300)
    ggsave(tsne_plot_path, tsne_plot, width = 16, height = 8, dpi = 300)
}

umap_plot_dir <- "umap_fibroblast_only_plots"
tsne_plot_dir <- "tsne_fibroblast_only_plots"
dir.create(umap_plot_dir, showWarnings = FALSE)
dir.create(tsne_plot_dir, showWarnings = FALSE)

# Subset to only fibroblast cells
fibroblast_cells <- subset(
    joined_integrated_data,
    singleR.labels_fine %in% fibroblast_celltypes
)

umap_plot <- DimPlot(
    fibroblast_cells,
    reduction = 'umap',
    group.by = 'singleR.labels_fine'
)

tsne_plot <- DimPlot(
    fibroblast_cells,
    reduction = 'tsne',
    group.by = 'singleR.labels_fine'
)

umap_plot_path <- file.path(umap_plot_dir, "all_fibroblast_celltypes_umap.pdf")
tsne_plot_path <- file.path(tsne_plot_dir, "all_fibroblast_celltypes_tsne.pdf")
ggsave(umap_plot_path, umap_plot, width = 16, height = 8, dpi = 300)
ggsave(tsne_plot_path, tsne_plot, width = 16, height = 8, dpi = 300)


for (cell_type in fibroblast_celltypes){
    fibroblast_cells$highlight <- ifelse(
        fibroblast_cells$singleR.labels_fine == cell_type,
        cell_type,
        "Other"
    )

    fibroblast_cells$highlight <- factor(
        fibroblast_cells$highlight,
        levels = c("Other", cell_type)
    )

    umap_plot <- DimPlot(
        fibroblast_cells,
        reduction = 'umap',
        group.by = 'highlight',
        cols = c("Other" = "lightgrey", setNames("#E41A1C", cell_type)),
        order = cell_type  # Plot highlighted cells on top
    )

    tsne_plot <- DimPlot(
        fibroblast_cells,
        reduction = 'tsne',
        group.by = 'highlight',
        cols = c("Other" = "lightgrey", setNames("#E41A1C", cell_type)),
        order = cell_type  # Plot highlighted cells on top
    )

    cell_type_label <- sanitize_label_for_filename(cell_type)
    umap_plot_path <- file.path(umap_plot_dir, paste0(cell_type_label, "_vs_all_umap.pdf"))
    tsne_plot_path <- file.path(tsne_plot_dir, paste0(cell_type_label, "_vs_all_tsne.pdf"))
    ggsave(umap_plot_path, umap_plot, width = 16, height = 8, dpi = 300)
    ggsave(tsne_plot_path, tsne_plot, width = 16, height = 8, dpi = 300)
}

print("umap plots done")


plot <- plotScoreHeatmap(pred)
ggsave("score_heatmap_fine.pdf", plot, width = 16, height = 8, dpi = 300)

plot <- plotDeltaDistribution(pred)
ggsave("delta_distribution_fine.pdf", plot, width = 16, height = 8, dpi = 300)

saveRDS(joined_integrated_data, file = "joined_integrated_seurat_object.rds")


#############################################################################################
## comment this out if not expecting fibroblasts
#
#fibroblast_barcodes <- rownames(joined_integrated_data@meta.data)[
#    joined_integrated_data$singleR.labels == "Fibroblasts"
#]
#
#fibroblast_subset <- subset(joined_integrated_data, cells = fibroblast_barcodes)
#
#print(paste("Number of fibroblast cells:", length(fibroblast_barcodes)))
#
#counts_fibroblasts <- GetAssayData(fibroblast_subset, assay = "RNA", layer = "data")
#
#pred_fibroblasts_fine <- SingleR(
#    test = counts_fibroblasts,
#    ref = ref,
#    labels = ref$label.fine  # More fine-grained labels
#)
#
## Add the fine-grained labels to the subset
#fibroblast_subset$singleR.labels.fine <- pred_fibroblasts_fine$labels[
#    match(rownames(fibroblast_subset@meta.data), rownames(pred_fibroblasts_fine))
#]
#
## Visualize
#plot <- DimPlot(fibroblast_subset, reduction = 'umap', 
#                group.by = c('seurat_clusters', 'singleR.labels.fine'))
#ggsave("fibroblast_subset_fine_labels.pdf", plot, width = 16, height = 8, dpi = 300)
#
## Score heatmap for fibroblast subtypes
#plot <- plotScoreHeatmap(pred_fibroblasts_fine)
#ggsave("fibroblast_score_heatmap_fine.pdf", plot, width = 16, height = 8, dpi = 300)
#
#plot <- plotDeltaDistribution(pred_fibroblasts_fine)
#ggsave("fibroblast_delta_distribution.pdf", plot, width = 16, height = 8, dpi = 300)
#############################################################################################
#


## Diff express between cell types
#Idents(joined_integrated_data) <- "singleR.labels_main"
#print(length(unique(joined_integrated_data$singleR.labels_main)))
#print(unique(joined_integrated_data$singleR.labels_main))
#for (label in unique(joined_integrated_data$singleR.labels_main)){
#    n_cells <- sum(joined_integrated_data$singleR.labels_main == label)
#    if (n_cells < args$min_cells){
#        print(paste0(label, " had less than ", args$min_cells," cells"))
#    }else{
#        print(label)
#        markers <- FindMarkers(joined_integrated_data, ident.1 = label, ident.2 = NULL)
#        write.csv(markers, paste0(sanitize_label_for_filename(label), "_vs_all_diffexpress_main.csv"))
#    }
#}
#
#markers <- FindMarkers(joined_integrated_data, ident.1 = "Fibroblasts", ident.2 = "Chondrocytes")
#write.csv(markers, "Fibroblasts_vs_Chondrocyte _diffexpress.csv")

outdir <- "diffexpress_one_vs_all"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

Idents(joined_integrated_data) <- "singleR.labels_fine"
print(length(unique(joined_integrated_data$singleR.labels_fine)))
print(unique(joined_integrated_data$singleR.labels_fine))
for (label in unique(joined_integrated_data$singleR.labels_fine)){
    n_cells <- sum(joined_integrated_data$singleR.labels_fine == label)
    if (n_cells < args$min_cells){
        print(paste0(label, " had less than ", args$min_cells," cells"))
    }else{
        print(label)
        markers <- FindMarkers(joined_integrated_data, ident.1 = label, ident.2 = NULL)
        file_path <- file.path(outdir, paste0(sanitize_label_for_filename(label), "_vs_all_diffexpress.csv"))
        write.csv(markers, file_path)
    }
}

outdir <- "diffexpress_fibroblasts_only_one_vs_all"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

Idents(fibroblast_cells) <- "singleR.labels_fine"
print(length(unique(fibroblast_cells$singleR.labels_fine)))
print(unique(fibroblast_cells$singleR.labels_fine))
for (label in unique(fibroblast_cells$singleR.labels_fine)){
    n_cells <- sum(fibroblast_cells$singleR.labels_fine == label)
    if (n_cells < args$min_cells){
        print(paste0(label, " had less than ", args$min_cells," cells"))
    }else{
        print(label)
        markers <- FindMarkers(fibroblast_cells, ident.1 = label, ident.2 = NULL)
        file_path <- file.path(outdir, paste0(sanitize_label_for_filename(label), "_vs_all_diffexpress_fibroblasts_only.csv"))
        write.csv(markers, file_path)
    }
}
