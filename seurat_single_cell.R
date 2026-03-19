library(argparse)
library(Seurat)
library(ggplot2)
library(SingleR)
library(celldex)
library(tidyverse)
library(pheatmap)
library(purrr)
library(EnhancedVolcano)
library(patchwork)


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


sanitize_label_for_filename <- function(label) {
    gsub("[/\\:*?\"<>|]", "_", label)
}

sanitize_label_for_AggregateExpression <- function(label) {
    gsub("[/\\:*?\"<>|_]", "-", label)
}

#make_multiplots <- function(integrated_data, cell_metadata, outdir, filenm, order = NULL){
#    umap_coords <- Embeddings(integrated_data, "umap.cca")
#    colnames(umap_coords) <- c("UMAP1", "UMAP2")
#    umap_coords$cell_metadata <- colData(integrated_data)[[cell_metadata]]
#    umap_coords$cluster <- integrated_data@clusters$UMAP$clusters
#
#    if (is.null(order)){
#        cell_metadata_uniq <- unique(umap_coords$cell_metadata)
#    }else{
#        cell_metadata_uniq <- order
#    }
#    highlight_plots <- lapply(cell_metadata_uniq, function(cm) {
#      umap_coords$highlight <- ifelse(umap_coords$cell_metadata == cm, "Selected", "Other")
#      # Plot "Other" first so highlighted cells are drawn on top
#      umap_coords <- umap_coords[order(umap_coords$highlight, decreasing = FALSE), ]
#      ggplot(umap_coords, aes(UMAP1, UMAP2, color = highlight)) +
#        geom_point(size = 0.3) +
#        scale_color_manual(values = c("Other" = "lightgrey", "Selected" = "red")) +
#        ggtitle(cm) +
#        theme_minimal() +
#        theme(legend.position = "none")
#    })
#
#    highlight_plot <- wrap_plots(highlight_plots, ncol = min(c(length(cell_metadata_uniq), 4)))
#
#    ggsave(file.path(outdir, filenm), 
#        highlight_plot,
#        width = 5 * min(c(length(cell_metadata_uniq), 4)),
#        height = 5 * ceiling(length(cell_metadata_uniq) / 4)
#    )
#}


make_multiplots <- function(integrated_data, cell_metadata, outdir, filenm, plot_order = NULL){
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

    umap_coords <- as.data.frame(Embeddings(integrated_data, "umap.cca"))
    colnames(umap_coords) <- c("UMAP1", "UMAP2")
    umap_coords$cell_metadata <- integrated_data@meta.data[[cell_metadata]]

    if (is.null(plot_order)){
        cell_metadata_uniq <- unique(umap_coords$cell_metadata)
    }else{
        cell_metadata_uniq <- plot_order
    }

    highlight_plots <- lapply(cell_metadata_uniq, function(cm) {
        umap_coords$highlight <- ifelse(umap_coords$cell_metadata == cm, "Selected", "Other")
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
        height = 5 * ceiling(length(cell_metadata_uniq) / 4),
        limitsize = FALSE
    )
}

# Create an ArgumentParser object
parser <- ArgumentParser(description = 'Integrate multiple scRNA-seq datasets into one Seurat object')
parser$add_argument("--mapping_file", help="File paths to mapping file with sample_id,path,diagnosis as header")
parser$add_argument("--min_cells", type="integer", help="Min cells per celltype for diff express", default = 8)

# Parse the command-line arguments
args <- parser$parse_args()

saved_objects_filename <- "scRNA_seurat_objects.rds"
if (file.exists(saved_objects_filename)) {
    seurat_objects <- readRDS(saved_objects_filename)
}else{
    # Read and process each dataset
    mapping_df <- read.csv(args$mapping_file)

    seurat_objects <- pmap(mapping_df, function(sample_id, path, diagnosis, dataset) {
        print(sample_id)
        print(path)
        print(diagnosis)
        print(dataset)

        data <- Read10X(data.dir = path)

        seurat_object <- CreateSeuratObject(counts = data)
        seurat_object <- AddMetaData(seurat_object, metadata = sample_id, col.name = "sample_id")
        seurat_object <- AddMetaData(seurat_object, metadata = diagnosis, col.name = "diagnosis")
        seurat_object <- AddMetaData(seurat_object, metadata = dataset, col.name = "dataset")
        print(head(seurat_object@meta.data))

        print(paste0("Number of barcodes before filtering ", length(seurat_object$orig.ident)))
        seurat_object$mitoPercent <- PercentageFeatureSet(seurat_object, pattern = '^MT-')
        seurat_object_filtered <- subset(seurat_object, subset = nCount_RNA > 800 &
            nFeature_RNA > 500 &
            mitoPercent < 10)
        print(paste0("Number of barcodes after filtering ", length(seurat_object_filtered$orig.ident)))

        return(seurat_object_filtered)
    })
  saveRDS(seurat_objects, saved_objects_filename)
}


# Integration steps
# Find integration anchors
saved_integrated_data <- "integrated_seurat_object.rds"
if (file.exists(saved_integrated_data)){
    # get rid of big objects
    rm(seurat_objects)
    gc()
    integrated_data <- readRDS(saved_integrated_data)
}else{
    object <- merge(seurat_objects[[1]], seurat_objects[-1])

    # get rid of big objects
    rm(seurat_objects)
    gc()

    print(object)
    object <- JoinLayers(object, assay = "RNA")
    object[["RNA"]] <- split(object[["RNA"]], f = object$sample_id)
    print(object)

    object <- NormalizeData(object)
    object <- FindVariableFeatures(object)
    object <- ScaleData(object)
    object <- RunPCA(object)

    object <- FindNeighbors(object, dims = 1:30, reduction = "pca")
    object <- FindClusters(object, resolution = 0.5, cluster.name = "unintegrated_clusters")

    object <- RunUMAP(object, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
    # visualize by batch and cell type annotation
    # cell type annotations were previously added by Azimuth
    plot <- DimPlot(object, reduction = "umap.unintegrated", group.by = c("sample_id", "diagnosis"))
    ggsave("unintegrated_data.pdf", plot, width = 16, height = 8, dpi = 300)

    integrated_data <- IntegrateLayers(
      object = object, method = CCAIntegration,
      orig.reduction = "pca", new.reduction = "integrated.cca",
      verbose = FALSE
    )


    # get rid of big objects
    rm(object)
    gc()

    integrated_data <- FindNeighbors(integrated_data, reduction = "integrated.cca", dims = 1:30)
    integrated_data <- FindClusters(integrated_data, cluster.name = "cca_clusters")

    integrated_data <- RunUMAP(integrated_data, reduction = "integrated.cca", dims = 1:30, reduction.name = "umap.cca")
    plot <- DimPlot(
      integrated_data,
      reduction = "umap.cca",
      group.by = c("sample_id", "diagnosis", "cca_clusters"),
      combine = FALSE, label.size = 2
    )

    ggsave("umap_cca.pdf", plot, width = 16, height = 8, dpi = 300)

    integrated_data <- RunTSNE(integrated_data, reduction = "integrated.cca", dims = 1:30, reduction.name = "tsne.cca")
    plot <- DimPlot(
      integrated_data,
      reduction = "tsne.cca",
      group.by = c("sample_id", "diagnosis", "cca_clusters"),
      combine = FALSE, label.size = 2
    )
    saveRDS(integrated_data, file = saved_integrated_data)

    outdir = "umaps"
    ggsave("tsne_cca.pdf", plot, width = 16, height = 8, dpi = 300)
    make_multiplots(integrated_data, "sample_id", outdir, "sample_id_ump.pdf")
    make_multiplots(integrated_data, "diagnosis", outdir, "diagnosis_umap.pdf")
    make_multiplots(integrated_data, "cca_clusters", outdir, "seurat_cca_clusters_umap.pdf")
    make_multiplots(integrated_data, "dataset", outdir, "dataset_umap.pdf")

}
print("here")

ref <- readRDS("hpca_ref.rds")

print(integrated_data[["RNA"]])

singleR_rds <- "singleR.rds"
if (file.exists(singleR_rds)){
    joined_integrated_data <- readRDS(singleR_rds)
}else{
    joined_integrated_data <- JoinLayers(integrated_data, assay = "RNA")
    print(head(joined_integrated_data))

    print(DefaultAssay(joined_integrated_data))
    DefaultAssay(joined_integrated_data) <- "RNA"
    print(DefaultAssay(joined_integrated_data))
    print(Layers(joined_integrated_data))

    counts <- GetAssayData(joined_integrated_data, assay = "RNA", layer = "data")

    pred <- SingleR(
    test = counts,
    ref = ref,
    labels = ref$label.fine,
    )

    joined_integrated_data$singleR.labels_fine <- pred$labels[match(rownames(joined_integrated_data@meta.data), rownames(pred))]

    pred <- SingleR(
    test = counts,
    ref = ref,
    labels = ref$label.main,
    )

    joined_integrated_data$singleR.labels_main <- pred$labels[match(rownames(joined_integrated_data@meta.data), rownames(pred))]

    print(paste0("Number of fine labels:  ", length(unique(joined_integrated_data$singleR.labels_fine))))
    write.csv(table(joined_integrated_data$singleR.labels_fine), "cell_type_counts_fine.csv")
    saveRDS(joined_integrated_data, file = singleR_rds)

    plot <- plotScoreHeatmap(pred)
    ggsave("score_heatmap_fine.pdf", plot, width = 16, height = 8, dpi = 300)

    plot <- plotDeltaDistribution(pred)
    ggsave("delta_distribution_fine.pdf", plot, width = 16, height = 8, dpi = 300)
}


# Make plots for each cell type

cell_types <- unique(joined_integrated_data$singleR.labels_fine)
umap_plot_dir <- "celltype_umap_plots"
tsne_plot_dir <- "celltype_tsne_plots"
dir.create(umap_plot_dir, showWarnings = FALSE)
dir.create(tsne_plot_dir, showWarnings = FALSE)

for (cell_type in cell_types){
    print(paste0("ploting ", cell_type))
    # Create a new column that highlights only this cell type
    joined_integrated_data$highlight <- ifelse(
        joined_integrated_data$singleR.labels_fine == cell_type,
        paste0(cell_type, "_", joined_integrated_data$diagnosis),
        "Other"
    )

    highlight_levels <- unique(joined_integrated_data$highlight)
    highlight_levels <- highlight_levels[highlight_levels != "Other"]
    
    # Reorder so "Other" cells are plotted first (behind highlighted cells)
    joined_integrated_data$highlight <- factor(
        joined_integrated_data$highlight,
        levels = c("Other", highlight_levels)
    )

    highlight_colors <- c("Other" = "lightgrey")
    # Add colors for each cell_type_diagnosis combination
    palette <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00")  # extend if needed
    for (i in seq_along(highlight_levels)) {
        highlight_colors[highlight_levels[i]] <- palette[i]
    }
    
    # Create the plot
    umap_plot <- DimPlot(
        joined_integrated_data,
        reduction = 'umap.cca',
        group.by = 'highlight',
        cols = highlight_colors,
        order = highlight_levels  # Plot highlighted cells on top
    )

    tsne_plot <- DimPlot(
        joined_integrated_data,
        reduction = 'tsne.cca',
        group.by = 'highlight',
        cols = highlight_colors,
        order = highlight_levels  # Plot highlighted cells on top
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


print("umap plots done")

saveRDS(joined_integrated_data, file = "joined_integrated_seurat_object.rds")

outdir <- "diffexpress_one_vs_all"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

joined_integrated_data$singleR.labels_fine <- sanitize_label_for_AggregateExpression(joined_integrated_data$singleR.labels_fine)
print(unique(joined_integrated_data$singleR.labels_fine))

bulk <- AggregateExpression(
    joined_integrated_data,
    return.seurat = TRUE,
    slot = "counts",
    assays = "RNA",
    group.by = c("singleR.labels_fine", "sample_id", "diagnosis")
)

Idents(bulk) <- "singleR.labels_fine"
print(unique(Idents(bulk)))
print(length(Idents(bulk)))
print(table(bulk$singleR.labels_fine))
write.csv(table(bulk$singleR.labels_fine), "aggregate_cell_counts.csv")


# add a pseudocount of 1 to every gene as per https://www.biostars.org/p/440379/
counts_matrix <- GetAssayData(bulk, layer = "counts")
bulk <- SetAssayData(bulk, layer = "counts", new.data = counts_matrix + 1)

for (label in unique(bulk$singleR.labels_fine)){
    n_cells <- sum(bulk$singleR.labels_fine == label)
    if (n_cells < args$min_cells){
        print(paste0(label, " had less than ", args$min_cells," cells"))
    }else{
        print(label)
        markers <- FindMarkers(
            bulk,
            ident.1 = label,
            ident.2 = NULL,
            slot = "counts",
            test.use = "DESeq2"
        )
        file_path <- file.path(outdir, paste0(sanitize_label_for_filename(label), "_vs_all_diffexpress.csv"))
        write.csv(markers, file_path)

        volcano_plot <- EnhancedVolcano(
            markers,
            lab = rownames(markers),
            x = "avg_log2FC",
            y = "p_val_adj",
            title = cell_type,
            pCutoff = 0.05,
            FCcutoff = 1
        )
        plot_path <- file.path(outdir, paste0(sanitize_label_for_filename(label), "_vs_all_diffexpress.pdf"))
        ggsave(plot_path, volcano_plot, width = 10, height = 8)
    }
}

#outdir <- "diffexpress_fibroblasts_only_one_vs_all"
#dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

#fibroblast_cells$singleR.labels_fine <- sanitize_label_for_AggregateExpression(fibroblast_cells$singleR.labels_fine)
#print(unique(fibroblast_cells$singleR.labels_fine))
#
#bulk <- AggregateExpression(
#    fibroblast_cells,
#    return.seurat = TRUE,
#    slot = "counts",
#    assays = "RNA",
#    group.by = c("singleR.labels_fine", "sample_id", "diagnosis")
#)
#
#Idents(bulk) <- "singleR.labels_fine"
#print(unique(Idents(bulk)))
#print(length(Idents(bulk)))
#
#
## add a pseudocount of 1 to every gene as per https://www.biostars.org/p/440379/
#counts_matrix <- GetAssayData(bulk, layer = "counts")
#bulk <- SetAssayData(bulk, layer = "counts", new.data = counts_matrix + 1)
#
#for (label in unique(bulk$singleR.labels_fine)){
#    n_cells <- sum(bulk$singleR.labels_fine == label)
#    if (n_cells < args$min_cells){
#        print(paste0(label, " had less than ", args$min_cells," cells"))
#    }else{
#        print(label)
#        markers <- FindMarkers(
#            bulk,
#            ident.1 = label,
#            ident.2 = NULL,
#            slot = "counts",
#            test.use = "DESeq2"
#        )
#        file_path <- file.path(outdir, paste0(sanitize_label_for_filename(label), "_vs_all_diffexpress_fibroblasts_only.csv"))
#        write.csv(markers, file_path)
#    }
#}
