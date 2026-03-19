library(argparse)
library(Seurat)
library(ggplot2)
library(SingleR)
library(celldex)
library(tidyverse)
library(pheatmap)
library(EnhancedVolcano)


sanitize_label_for_filename <- function(label) {
    gsub("[/\\:*?\"<>|]", "_", label)
}

sanitize_label_for_AggregateExpression <- function(label) {
    gsub("[/\\:*?\"<>|_]", "-", label)
}

run_DE <- function(seurat_object, cell_type, identity_of_interest, outdir, file_suff, min_cells){
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
    print(file_suff)

    bulk <- AggregateExpression(
        seurat_object,
        return.seurat = TRUE,
        slot = "counts",
        assays = "RNA",
        group.by = c("diagnosis", "sample_id")
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
            ident.2 = "normal-skin",
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
parser$add_argument("--outdir", help="directory to write to", required = TRUE)
parser$add_argument("--min_cells", type="integer", help="Min cells per celltype for diff express", default = 5)
parser$add_argument("--cell_type", help="cell type to use", required = TRUE)

# Parse the command-line arguments
args <- parser$parse_args()

print(args)

joined_integrated_data <- readRDS(args$joined_integrated_seurat_object)

wanted_groups = c("fibroblasts", "pericyte")
just_subset <- subset(joined_integrated_data, cell_group %in% wanted_groups)

rm(joined_integrated_seurat_object)
gc()

just_subset$diagnosis <- trimws(just_subset$diagnosis)
just_subset$singleR.labels_fine <- sanitize_label_for_AggregateExpression(just_subset$singleR.labels_fine)

#just_subset_diag <- subset(just_subset, diagnosis == args$cell_type)
for (label_fine in unique(just_subset$singleR.labels_fine)){
    just_subset_label <- subset(just_subset, singleR.labels_fine == label_fine)
    run_DE(just_subset_label, args$cell_type, "diagnosis", file.path(args$outdir, args$cell_type), paste0("_vs_normal_skin_", label_fine), args$min_cells)
}
