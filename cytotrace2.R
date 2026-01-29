library(CytoTRACE2)
library(argparse)
library(Seurat)
library(ggplot2)
library(patchwork)
library(readr)
library(tibble)

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

run_DE <- function(seurat_object, potency, outdir, min_cells){
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

    bulk <- AggregateExpression(
        seurat_object,
        return.seurat = TRUE,
        slot = "counts",
        assays = "RNA",
        group.by = c("singleR.labels_fine", "sample_id", "diagnosis", "CytoTRACE2_Potency")
    )

    print("here1")
    Idents(bulk) <- "CytoTRACE2_Potency"
    print(unique(Idents(bulk)))
    print(length(Idents(bulk)))
    print(table(bulk$CytoTRACE2_Potency))

    print("here2")
    file_name <- file.path(outdir, "aggregate_cell_counts.csv")
    write.csv(table(bulk$CytoTRACE2_Potency), file_name)


    print("here3")
    # add a pseudocount of 1 to every gene as per https://www.biostars.org/p/440379/
    counts_matrix <- LayerData(bulk, layer = "counts")
    LayerData(bulk, layer = "counts") <- counts_matrix + 1

    print("here4")
    n_cells <- sum(bulk$CytoTRACE2_Potency == potency)
    if (n_cells < min_cells){
        print(paste0(CytoTRACE2_Potency, " had less than ", min_cells," cells"))
        print(unique(bulk$CytoTRACE2_Potency))
    }else{
        print(potency)
        markers <- FindMarkers(
            bulk,
            ident.1 = potency,
            ident.2 = NULL,
            slot = "counts",
            test.use = "DESeq2"
        )
        file_path <- file.path(outdir, paste0(sanitize_label_for_filename(potency), "_vs_all_fibroblasts_diffexpress.csv"))
        write.csv(markers, file_path)
    }
}

plot_data <- function(cytotrace2_result, annotation, expression_data, batch){
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
}

parser <- ArgumentParser(description = 'Integrate multiple scRNA-seq datasets into one Seurat object')
parser$add_argument("--rds", type="character", help="File path to joined_integrated_seurat_object.rds")
args <- parser$parse_args()

data <- readRDS(args$rds)

if (length(Layers(data[["RNA"]], search = "counts")) > 1) {
    cat("Joining counts layers...\n")
    data <- JoinLayers(data, assay = "RNA", layers = "counts")
}

#fibroblast_cells <- subset(
#    data,
#    singleR.labels_fine %in% fibroblast_celltypes
#)

all_potency <- data.frame()

batches <- unique(data$sample_id)
for (batch in batches) {
    cat("Processing batch:", batch, "\n")
    
    # Subset to this batch
    batch_data <-  subset(
        data,
        sample_id == batch
    )

    print(head(batch_data))

    # Extract expression data
    expression_data <- as.data.frame(GetAssayData(batch_data, assay = "RNA", layer = "counts"))
    
    # Run CytoTRACE2
    cytotrace2_result <- cytotrace2(expression_data, species = "human")
    
    potency_batch <- data.frame(
        CytoTRACE2_Potency = cytotrace2_result$CytoTRACE2_Potency,
        row.names = rownames(cytotrace2_result)
    )
    print("binding")
    all_potency <- rbind(all_potency, potency_batch)
    print("finished binding")
    
    # Prepare annotation
    annotation <- data.frame(Phenotype = batch_data$singleR.labels_fine)
    
    # Write stats
    cytotrace2_result_with_rownames <- merge(cytotrace2_result, annotation, by = 0, all = TRUE)
    names(cytotrace2_result_with_rownames)[1] <- "barcode"
    write_csv(cytotrace2_result_with_rownames, paste0("cytotrace_stats_", batch, "_fine.csv"))

#    plot_data(cytotrace2_result, annotation, expression_data, paste0(batch, "_fine"))
    
    cat("Saved: cytotrace_plots_", batch, ".pdf\n")
}

head(rownames(data@meta.data))
head(rownames(all_potency))

# Check for overlap
print(sum(rownames(data@meta.data) %in% rownames(all_potency)))


data$CytoTRACE2_Potency <- all_potency$CytoTRACE2_Potency[match(rownames(data@meta.data), rownames(all_potency))]

saveRDS(fibroblast_cells, "cytotrace_seurat.rds")

#fibroblast_cells <- readRDS("cytotrace_seurat.rds")
#
#
#for (celltype in unique(fibroblast_cells$singleR.labels_fine)){
#    cells <- subset(
#        fibroblast_cells,
#        singleR.labels_fine == celltype
#    )
#
#    plot <- DimPlot(
#      cells,
#      reduction = "umap.cca",
#      group.by = "CytoTRACE2_Potency",
#      label.size = 2
#    )
#    ggsave(paste0("umap_CytoTRACE2_Potency_", celltype, ".pdf"), plot, width = 16, height = 8, dpi = 300)
#
#    plot <- DimPlot(
#      cells,
#      reduction = "tsne.cca",
#      group.by = "CytoTRACE2_Potency",
#      label.size = 2
#    )
#    ggsave(paste0("tsne_CytoTRACE2_Potency_", celltype, ".pdf"), plot, width = 16, height = 8, dpi = 300)
#}
#
#
##for (potency in unique(fibroblast_cells$CytoTRACE2_Potency)){
##    run_DE(fibroblast_cells, potency, paste0(potency, "_vs_all_others_just_fibroblasts"), 8)
##}
