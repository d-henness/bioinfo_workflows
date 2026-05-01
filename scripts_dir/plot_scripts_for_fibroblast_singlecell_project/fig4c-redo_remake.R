library(argparse)
library(Seurat)
library(ggplot2)
library(tidyverse)
library(pheatmap)
library(purrr)
library(patchwork)
library(BBmisc)
library(ggrepel)
library(ggpubr)

mapping <- read.csv("cell_type_mapping.csv")

gene_list <- c(
    "ARID1A",
    "ARID2",
    "ARID5B",
    "KMT2C",
    "SETD2",
    "EZH2",
    "CHD1",
    "EP300",
    "BRD4",
    "HEXIM1",
    "CUL3",
    "DICER1"
)

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

# Log2 transform
just_subset@assays$RNA@layers$data <- as.matrix(log2(just_subset@assays$RNA@layers$counts + 1))

make_gene_plot <- function(df, gene, label_col_title) {
  ggplot(df, aes(x = RNA_Counts, y = .data[[gene]])) +
    geom_point(aes(color = Diagnosis), size = 3) +
    geom_smooth(color = "black", fill = "grey90", linetype = "dashed", method = "lm", se = TRUE) +
    stat_cor(method = "spearman", label.x = 0.5, label.y = 9.5, size = 4) +
    geom_text_repel(
      aes(label = ifelse(RNA_Counts > 6 & .data[[gene]] > 5.5, as.character(Celltype), ""), color = Diagnosis),
      size = 4, point.padding = 0.1, box.padding = 0.5, segment.alpha = 1
    ) +
    theme_classic() +
    theme(
      aspect.ratio = 1,
      axis.text = element_text(size = 12, colour = "black"),
      axis.title = element_text(size = 12)
    ) +
    labs(y = gene, x = "Log2 Transcripts", title = label_col_title)
}

prep_gene_df <- function(gene, label_col) {
  df <- data.frame(
    gene_expr = as.numeric(as.matrix(just_subset@assays$RNA@layers$data)[rownames(just_subset@assays$RNA) == gene, ]),
    RNA_Counts = log2(just_subset@meta.data$nCount_RNA),
    Celltype = label_col,
    Diagnosis = just_subset@meta.data$diagnosis
  )
  colnames(df)[1] <- gene

  df[, "CelltypeDiag"] <- paste(df$Celltype, "/", df$Diagnosis, sep = "")
  df$Celltype <- NULL
  df$Diagnosis <- NULL

  df <- aggregate(. ~ CelltypeDiag, df, mean)
  df <- BBmisc::normalize(df, method = "range", range = c(0, 10), margin = 2L)
  rownames(df) <- df$CelltypeDiag
  df$CelltypeDiag <- NULL

  df[, "celltypediag"] <- rownames(df)
  df[, "Diagnosis"] <- sapply(X = strsplit(df$celltypediag, split = "/"), FUN = "[", 2)
  df[, "Celltype"] <- sapply(X = strsplit(df$celltypediag, split = "/"), FUN = "[", 1)
  df$celltypediag <- paste(df$Celltype, " (", df$Diagnosis, ")", sep = "")
  rownames(df) <- df$celltypediag
  df
}

for (gene in gene_list) {
  df1 <- prep_gene_df(gene, just_subset@meta.data$m_labels)
  df2 <- prep_gene_df(gene, just_subset@meta.data$new_cat_labels)

  p1 <- make_gene_plot(df1, gene, "m_labels")
  p2 <- make_gene_plot(df2, gene, "new_cat_labels")

  combined <- p1 + p2
  ggsave(file.path(args$outdir, paste0(gene, ".pdf")), combined, width = 12, height = 6)
}
