options(warn = 1)
library(CytoTRACE2)
library(argparse)
library(Seurat)
library(ggplot2)
library(patchwork)
library(readr)
library(tibble)

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
parser$add_argument("--outdir", type="character", help="where to write results")
parser$add_argument("--ncores", type="integer", help="number of cores to use")
args <- parser$parse_args()

dir.create(args$outdir, showWarnings = FALSE, recursive = TRUE)
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
    cat("Batch:", batch, "\n")
    cat("Dims:", nrow(expression_data), "genes x", ncol(expression_data), "cells\n")
    cat("Class:", class(expression_data), "\n")
    cat("Type of first element:", class(expression_data[1,1]), "\n")
    cat("Any NAs:", any(is.na(expression_data)), "\n")
    cat("Any Inf:", any(is.infinite(as.matrix(expression_data[1:10, 1:10]))), "\n")
    cat("Sum zero columns:", sum(colSums(expression_data) == 0), "\n")
    cat("Sum zero rows:", sum(rowSums(expression_data) == 0), "\n")
    
    # Run CytoTRACE2
#    cytotrace2_result <- cytotrace2(expression_data, species = "human", ncores = args$ncores)
    cytotrace2_result <- tryCatch({
        cytotrace2(expression_data, species = "human", ncores = args$ncores)
    }, error = function(e) {
        cat("Failed on batch:", batch, "—", conditionMessage(e), "\n")
        return(NULL)
    })
    if (is.null(cytotrace2_result)) next
    
    potency_batch <- data.frame(
        CytoTRACE2_Potency = cytotrace2_result$CytoTRACE2_Potency,
        CytoTRACE2_Score = cytotrace2_result$CytoTRACE2_Score,
        row.names = rownames(cytotrace2_result)
    )
    print("binding")
    all_potency <- rbind(all_potency, potency_batch)
    print("finished binding")
    
    # absolute scaling
    batch_data@assays$RNA@layers$abs_scale <- as.matrix(log2(LayerData(batch_data, assay = "RNA", layer = "counts") + 1))
    batch_data@meta.data$scaledcounts <- colSums(batch_data@assays$RNA@layers$abs_scale)

    batch_data@meta.data <- rownames_to_column(batch_data@meta.data, var = "cell_barcode")

    # Prepare annotation
    meta <- batch_data@meta.data[, c("cell_barcode", "singleR.labels_fine", "scaledcounts", "diagnosis", "sample_id", "dataset")]
    cytotrace2_result_with_rownames <- merge(cytotrace2_result, meta, by.x = 0, by.y = "cell_barcode", all = TRUE)
    names(cytotrace2_result_with_rownames)[1] <- "barcode"
    write_csv(cytotrace2_result_with_rownames, file.path(args$outdir, paste0("cytotrace_stats_", batch, "_fine.csv")))
#    plot_data(cytotrace2_result, annotation, expression_data, paste0(batch, "_fine"))
    
    cat("Saved: cytotrace_plots_", batch, ".pdf\n")
    # memory management
    rm(batch_data, expression_data, cytotrace2_result, potency_batch, meta, cytotrace2_result_with_rownames)
    gc()
}

head(rownames(data@meta.data))
head(rownames(all_potency))

# Check for overlap
print(sum(rownames(data@meta.data) %in% rownames(all_potency)))

print(head(data))
print(head(all_potency))

data <- AddMetaData(data, all_potency)

saveRDS(data, "cytotrace_seurat.rds")

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
#    plot <- Dimdata#      reduction = "tsne.cca",
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
