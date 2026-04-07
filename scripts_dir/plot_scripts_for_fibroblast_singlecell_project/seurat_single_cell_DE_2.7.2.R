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

    bulk <- AggregateExpression(
        seurat_object,
        return.seurat = TRUE,
        slot = "counts",
        assays = "RNA",
        group.by = c("CytoTRACE2_Potency", "sample_id")
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
parser$add_argument("--outdir", help="directory to write to", required = TRUE)
parser$add_argument("--min_cells", type="integer", help="Min cells per celltype for diff express", default = 5)
parser$add_argument("--cell_type", help="cell type to use", required = TRUE)

# Parse the command-line arguments
args <- parser$parse_args()

print(args)

joined_integrated_data <- readRDS("rds_files/cytotrace_seurat.rds")

wanted_groups = c("fibroblasts", "pericyte")
just_subset <- subset(joined_integrated_data, cell_group %in% wanted_groups)
just_subset <- subset(just_subset, !is.na(CytoTRACE2_Potency))

for (label_fine in unique(just_subset$CytoTRACE2_Potency)){
    run_DE(just_subset, label_fine, "CytoTRACE2_Potency", args$outdir, "vs_all", args$min_cells)
}


rm(joined_integrated_seurat_object)
gc()

just_subset$diagnosis <- trimws(just_subset$diagnosis)

just_subset_diag <- subset(just_subset, diagnosis == args$cell_type)
for (label_fine in unique(just_subset$CytoTRACE2_Potency)){
    run_DE(just_subset_diag, label_fine, "CytoTRACE2_Potency", file.path(args$outdir, args$cell_type), "vs_all", args$min_cells)
}
