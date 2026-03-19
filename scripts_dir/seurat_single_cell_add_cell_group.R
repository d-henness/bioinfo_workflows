
library(argparse)
library(Seurat)
library(ggplot2)
library(tidyverse)
library(pheatmap)

parser <- ArgumentParser(description = 'Integrate multiple scRNA-seq datasets into one Seurat object')
parser$add_argument("joined_integrated_seurat_object_with_cytotrace_labels", help="File paths to joined_integrated_seurat_object.rds")
parser$add_argument("--group_file", help="file to group cell types", required = TRUE)
args <- parser$parse_args()

seurat_object <- readRDS(args$joined_integrated_seurat_object_with_cytotrace_labels)
group_data <- read.csv(args$group_file, sep = '\t')

lookup <- setNames(group_data$new_name_dotplot, group_data$singleR.label_fine)

seurat_object$cell_group <- ifelse(
    seurat_object$singleR.labels_fine %in% names(lookup),
    lookup[as.character(seurat_object$singleR.labels_fine)],
    "n/a"
)
print(head(seurat_object))

saveRDS(seurat_object, file = "joined_integrated_seurat_object_with_group_data.rds")
