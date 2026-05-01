library(argparse)
library(Seurat)
library(ggplot2)
library(tidyverse)
library(pheatmap)
library(purrr)
library(patchwork)

mapping <- read.csv("cell_type_mapping.csv")

sanitize_label_for_filename <- function(label) {
    gsub("[/\\:*?\"<>|]", "_", label)
}

mangle_label <- function(label) {
  gsub("[^a-zA-Z0-9]", "-", label)
}

fix_colnames <- function(x) {
  mapping <- c(
    "Hair-follicle"   = "Hair_follicle"
  )
  ifelse(x %in% names(mapping), mapping[x], x)
}

write_ave_express_by_percent <- function(df, identity, filenm1){
    ave_express <- AverageExpression(df, group.by = identity)

    if (identity != 'diagnosis'){
        original_labels <- mapping[[identity]]

        mangled <- sapply(as.character(original_labels), mangle_label)
        name_map <- setNames(original_labels, mangled)

        colnames(ave_express$RNA) <- fix_colnames(colnames(ave_express$RNA))
    }

    df_out <- as.data.frame(ave_express$RNA, check.names = FALSE)
    df_out <- df_out[, sort(colnames(df_out))]
    df_out <- log(df_out + 1)
    df_out <- df_out[c("MUC1", "MUC16"), ]

    write.table(df_out, filenm1, sep = ",", col.names = NA)
}

parser <- ArgumentParser(description = 'Integrate multiple scRNA-seq datasets into one Seurat cbject')
parser$add_argument("--outdir", help="directory to write to", required = TRUE)
args <- parser$parse_args()
print(args)
dir.create(args$outdir, showWarnings = FALSE, recursive = TRUE)

joined_integrated_data <- readRDS("rds_files/joined_integrated_seurat_object_with_group_data.rds")
wanted_groups = c("fibroblasts", "pericyte")

just_subset <- subset(joined_integrated_data, cell_group %in% wanted_groups)
just_subset$m_labels <- mapping$m_labels[match(just_subset$singleR.labels_fine, mapping$singleR.labels_fine)]
just_subset$new_cat_labels <- mapping$new_cat_labels[match(just_subset$singleR.labels_fine, mapping$singleR.labels_fine)]

filenm1 <- file.path(args$outdir, "average_expression_new_cat_labels_all.csv")
write_ave_express_by_percent(just_subset, "new_cat_labels", filenm1)
filenm1 <- file.path(args$outdir, "average_expression_m_labels_all.csv")
write_ave_express_by_percent(just_subset, "m_labels", filenm1)
filenm1 <- file.path(args$outdir, "average_expression_diagnosis_all.csv")
write_ave_express_by_percent(just_subset, "diagnosis", filenm1)

for (diag in unique(just_subset$diagnosis)){
    just_diag <- subset(just_subset, diagnosis == diag)
    filenm1 <- file.path(args$outdir, paste0("average_expression_new_cat_labels_", diag, ".csv"))
    write_ave_express_by_percent(just_diag, "new_cat_labels", filenm1)
    filenm1 <- file.path(args$outdir, paste0("average_expression_m_labels_", diag, ".csv"))
    write_ave_express_by_percent(just_diag, "m_labels", filenm1)
}
