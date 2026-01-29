library(argparse)
library(Seurat)
library(ggplot2)
library(SingleR)
library(celldex)
library(tidyverse)
library(pheatmap)
library(EnhancedVolcano)


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

# hard code for now, do properly later
diagnosis = c(
    "SC1_cellranger" = "control",
    "SC4_cellranger" = "control",
    "SC18_cellranger" = "control",
    "SC32_cellranger" = "control",
    "SC33_cellranger" = "control",
    "SC34_cellranger" = "control",
    "SC50_cellranger" = "control",
    "SC68_cellranger" = "control",
    "SC124_cellranger" = "control",
    "SC125_cellranger" = "control",
    "SC2_cellranger" = "SSC",
    "SC5_cellranger" = "SSC",
    "SC19_cellranger" = "SSC",
    "SC49_cellranger" = "SSC",
    "SC60_cellranger" = "SSC",
    "SC69_cellranger" = "SSC",
    "SC70_cellranger" = "SSC",
    "SC86_cellranger" = "SSC",
    "SC119_cellranger" = "SSC",
    "SC185_cellranger" = "SSC",
    "SC188_cellranger" = "SSC",
    "SC189_cellranger" = "SSC"
)

sanitize_label_for_filename <- function(label) {
    gsub("[/\\:*?\"<>|]", "_", label)
}

sanitize_label_for_AggregateExpression <- function(label) {
    gsub("[/\\:*?\"<>|_]", "-", label)
}

run_DE <- function(seurat_object, cell_type, identity_of_interest, outdir, file_suff, min_cells){
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

    bulk <- AggregateExpression(
        seurat_object,
        return.seurat = TRUE,
        slot = "counts",
        assays = "RNA",
        group.by = c("singleR.labels_fine", "sample_id", "diagnosis", "is_fib")
    )

    Idents(bulk) <- identity_of_interest
    print(head(bulk))

    file_name <- file.path(outdir, paste0("aggregate_cell_counts_", file_suff, ".csv"))
    write.csv(table(bulk[[identity_of_interest]]), file_name)


    # add a pseudocount of 1 to every gene as per https://www.biostars.org/p/440379/
    counts_matrix <- GetAssayData(bulk, layer = "counts")
    bulk <- SetAssayData(bulk, layer = "counts", new.data = counts_matrix + 1)

    n_cells <- sum(bulk[[identity_of_interest, drop = TRUE]] == cell_type)
    print(n_cells)
    if (n_cells < min_cells){
        print(paste0(cell_type, " had less than ", min_cells," cells"))
        print(table(bulk[[identity_of_interest]]))
    }else{
        print(cell_type)
        markers <- FindMarkers(
            bulk,
            ident.1 = cell_type,
            ident.2 = NULL,
            slot = "counts",
            test.use = "DESeq2"
        )
        file_path <- file.path(outdir, paste0(sanitize_label_for_filename(cell_type), "_", file_suff, ".csv"))
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
        plot_path <- file.path(outdir, paste0(sanitize_label_for_filename(cell_type), "_", file_suff, "_volcano.pdf"))
        ggsave(plot_path, volcano_plot, width = 10, height = 8)
    }
}

# Create an ArgumentParser object
parser <- ArgumentParser(description = 'Integrate multiple scRNA-seq datasets into one Seurat object')
parser$add_argument("joined_integrated_seurat_object", help="File paths to joined_integrated_seurat_object.rds")
parser$add_argument("--cell_type", help="cell type to use", required = TRUE)
parser$add_argument("--cell_type_identity", help="identity to use", required = TRUE)
parser$add_argument("--outdir", help="directory to write to", required = TRUE)
parser$add_argument("--file_suff", help="file suffix", required = TRUE)
parser$add_argument("--min_cells", type="integer", help="Min cells per celltype for diff express", default = 8)
parser$add_argument("--subset", help="cell type to use, leaving this null subsets to all fibroblast cell types", default = NULL)
parser$add_argument("--subset_identity", help="identity to use", default = NULL)

# Parse the command-line arguments
args <- parser$parse_args()

print(args)

joined_integrated_data <- readRDS(args$joined_integrated_seurat_object)
joined_integrated_data$diagnosis <- trimws(joined_integrated_data$diagnosis)
joined_integrated_data$singleR.labels_fine <- sanitize_label_for_AggregateExpression(joined_integrated_data$singleR.labels_fine)

sanitized_fib_labels <- sanitize_label_for_AggregateExpression(fibroblast_celltypes)
joined_integrated_data$is_fib <- ifelse(
    joined_integrated_data$singleR.labels_fine %in% sanitized_fib_labels,
    "fibroblast",
    "other"
)

cells <- joined_integrated_data
#if (is.null(args$subset_identity)){
#    sanitized_fib_labels <- sanitize_label_for_AggregateExpression(fibroblast_celltypes)
#    cells <- subset(
#        joined_integrated_data,
#        singleR.labels_fine %in% sanitized_fib_labels
#    )
#}else{
#    cells <- joined_integrated_data[, joined_integrated_data[[args$subset_identity]] == args$subset]
#}

print(head(cells))
print(unique(cells$singleR.labels_fine))

run_DE(cells, args$cell_type, args$cell_type_identity, args$outdir, args$file_suff, args$min_cells)
