library(argparse)
library(Seurat)
library(ggplot2)
library(tidyverse)
library(pheatmap)



# Create an ArgumentParser object
parser <- ArgumentParser(description = 'Integrate multiple scRNA-seq datasets into one Seurat object')
parser$add_argument("joined_integrated_seurat_object_with_cytotrace_labels", help="File paths to joined_integrated_seurat_object.rds")
parser$add_argument("--outdir", help="directory to write to", required = TRUE)
parser$add_argument("--file_suff", help="file suffix", required = TRUE)

# Parse the command-line arguments
args <- parser$parse_args()

print(args)

seurat_object <- readRDS(args$joined_integrated_seurat_object_with_cytotrace_labels)
seurat_object$cell_group[seurat_object$cell_group == "melanocyes"] <- "melanocytes"
seurat_object$diagnosis <- trimws(seurat_object$diagnosis)
print(head(seurat_object))

dir.create(args$outdir, showWarnings = FALSE, recursive = TRUE)
date_and_time <- format(now(), "%Y%m%d_%H%M%S")

just_subset <- subset(seurat_object, cell_group != 'n/a')

filenm <- file.path(
    args$outdir,
    paste0("all_datasets_cell_group_", args$file_suff, "_", date_and_time, ".csv")
)

write.csv(table(just_subset$cell_group, just_subset$diagnosis), filenm)

filenm <- file.path(
    args$outdir,
    paste0("all_datasets_singleR.labels_fine_", args$file_suff, "_", date_and_time, ".csv")
)

write.csv(table(just_subset$singleR.labels_fine, just_subset$diagnosis), filenm)
