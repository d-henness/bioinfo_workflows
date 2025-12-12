library(argparse)
library(Seurat)
library(ggplot2)
library(SingleR)
library(celldex)
library(tidyverse)
library(pheatmap)

# Create an ArgumentParser object
parser <- ArgumentParser(description = 'Integrate multiple scRNA-seq datasets into one Seurat object')

# Add arguments to specify the locations of datasets
parser$add_argument("--datasets", nargs='+', type="character", help="File paths to datasets")

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

# Visualization
plot <- DimPlot(integrated_data, reduction = "umap", label = TRUE, group.by = c("orig.ident", "seurat_clusters"))
ggsave("integrated_data.pdf", plot, width = 16, height = 8, dpi = 300)



ref <- celldex::HumanPrimaryCellAtlasData()
print(head(as.data.frame(colData(ref))))

joined_integrated_data <- JoinLayers(integrated_data, assay = "RNA", layers = "data")
print(head(joined_integrated_data))

counts <- GetAssayData(joined_integrated_data, assay = "RNA", layer = "data")

pred <- SingleR(
    test = counts,
    ref = ref,
    labels = ref$label.main
)

joined_integrated_data$singleR.labels <- pred$labels[match(rownames(joined_integrated_data@meta.data), rownames(pred))]
plot <- DimPlot(joined_integrated_data, reduction = 'umap', group.by = c('seurat_clusters', 'singleR.labels'))
ggsave("joined_integrated_data_tissue_labels.pdf", plot, width = 16, height = 8, dpi = 300)

plot <- plotScoreHeatmap(pred)
ggsave("score_heatmap.pdf", plot, width = 16, height = 8, dpi = 300)

plot <- plotDeltaDistribution(pred)
ggsave("delta_distribution.pdf", plot, width = 16, height = 8, dpi = 300)

saveRDS(joined_integrated_data, file = "joined_integrated_seurat_object.rds")
