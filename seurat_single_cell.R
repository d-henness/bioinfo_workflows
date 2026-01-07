library(argparse)
library(Seurat)
library(ggplot2)
library(SingleR)
library(celldex)
library(tidyverse)
library(pheatmap)

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
        RunUMAP(reduction = "pca", dims = 1:30)
    saveRDS(integrated_data, file = saved_integrated_data)
}

print(length(unique(integrated_data$seurat_clusters)))
print(unique(integrated_data$seurat_clusters))

# Visualization
plot <- DimPlot(integrated_data, reduction = "umap", label = TRUE, group.by = c("orig.ident", "seurat_clusters"))
ggsave("integrated_data.pdf", plot, width = 16, height = 8, dpi = 300)

ref <- celldex::HumanPrimaryCellAtlasData()

joined_integrated_data <- JoinLayers(integrated_data, assay = "RNA", layers = "data")
print(head(joined_integrated_data))

counts <- GetAssayData(joined_integrated_data, assay = "RNA", layer = "data")

pred <- SingleR(
    test = counts,
    ref = ref,
    labels = ref$label.main
)


joined_integrated_data$singleR.labels_main <- pred$labels[match(rownames(joined_integrated_data@meta.data), rownames(pred))]
print(paste0("Number of main labels:  ", length(unique(joined_integrated_data$singleR.labels_main))))
write.csv(table(joined_integrated_data$singleR.labels_main), "cell_type_counts_main.csv")

plot <- DimPlot(joined_integrated_data, reduction = 'umap', group.by = c('seurat_clusters', 'singleR.labels_main'))
ggsave("joined_integrated_data_tissue_labels_main.pdf", plot, width = 16, height = 8, dpi = 300)

plot <- plotScoreHeatmap(pred)
ggsave("score_heatmap_main.pdf", plot, width = 16, height = 8, dpi = 300)

plot <- plotDeltaDistribution(pred)
ggsave("delta_distribution_main.pdf", plot, width = 16, height = 8, dpi = 300)

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


# Diff express between cell types
Idents(joined_integrated_data) <- "singleR.labels_main"
print(length(unique(joined_integrated_data$singleR.labels_main)))
print(unique(joined_integrated_data$singleR.labels_main))
for (label in unique(joined_integrated_data$singleR.labels_main)){
    n_cells <- sum(joined_integrated_data$singleR.labels_main == label)
    if (n_cells < args$min_cells){
        print(paste0(label, " had less than ", args$min_cells," cells"))
    }else{
        print(label)
        markers <- FindMarkers(joined_integrated_data, ident.1 = label, ident.2 = NULL)
        write.csv(markers, paste0(sanitize_label_for_filename(label), "_vs_all_diffexpress_main.csv"))
    }
}

markers <- FindMarkers(joined_integrated_data, ident.1 = "Fibroblasts", ident.2 = "Chondrocytes")
write.csv(markers, "Fibroblasts_vs_Chondrocyte _diffexpress.csv")

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
        write.csv(markers, paste0(sanitize_label_for_filename(label), "_vs_all_diffexpress.csv"))
    }
}

markers <- FindMarkers(joined_integrated_data, ident.1 = "Fibroblasts", ident.2 = "Chondrocytes")
write.csv(markers, "Fibroblasts_vs_Chondrocyte _diffexpress_fine.csv")

