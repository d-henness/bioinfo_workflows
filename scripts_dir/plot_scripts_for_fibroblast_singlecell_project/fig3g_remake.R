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

# Build comparison dataframe
GENEcomparison <- data.frame(
  Chd1 = as.numeric(as.matrix(just_subset@assays$RNA@layers$data)[rownames(just_subset@assays$RNA) == "CHD1", ]),
  RNA_Counts = log2(just_subset@meta.data$nCount_RNA)
)

# --- Plot 1: m_labels ---
GENEcomp1 <- GENEcomparison
GENEcomp1[, "Celltype"] <- just_subset@meta.data$m_labels
GENEcomp1[, "Diagnosis"] <- just_subset@meta.data$diagnosis

GENEcomp1[, "CelltypeDiag"] <- paste(GENEcomp1$Celltype, "/", GENEcomp1$Diagnosis, sep = "")
GENEcomp1$Celltype <- NULL
GENEcomp1$Diagnosis <- NULL

GENEcomp1 <- aggregate(. ~ CelltypeDiag, GENEcomp1, mean)
GENEcomp1 <- BBmisc::normalize(GENEcomp1, method = "range", range = c(0, 10), margin = 2L)
rownames(GENEcomp1) <- GENEcomp1$CelltypeDiag
GENEcomp1$CelltypeDiag <- NULL

GENEcomp1[, "celltypediag"] <- rownames(GENEcomp1)
GENEcomp1[, "Diagnosis"] <- sapply(X = strsplit(GENEcomp1$celltypediag, split = "/"), FUN = "[", 2)
GENEcomp1[, "Celltype"] <- sapply(X = strsplit(GENEcomp1$celltypediag, split = "/"), FUN = "[", 1)
GENEcomp1$celltypediag <- paste(GENEcomp1$Celltype, " (", GENEcomp1$Diagnosis, ")", sep = "")
rownames(GENEcomp1) <- GENEcomp1$celltypediag

p1 <- ggplot(GENEcomp1, aes(x = RNA_Counts, y = Chd1)) +
  geom_point(aes(color = Diagnosis), size = 3) +
  geom_smooth(color = "black", fill = "grey90", linetype = "dashed", method = "lm", se = TRUE) +
  stat_cor(method = "spearman", label.x = 0.5, label.y = 9.5, size = 4) +
  geom_text_repel(
    aes(label = ifelse(RNA_Counts > 6 & Chd1 > 5.5, as.character(Celltype), ""), color = Diagnosis),
    size = 4, point.padding = 0.1, box.padding = 0.5, segment.alpha = 1
  ) +
  theme_classic() +
  theme(
    aspect.ratio = 1,
    axis.text = element_text(size = 12, colour = "black"),
    axis.title = element_text(size = 12)
  ) +
  labs(y = "Chd1", x = "Log2 Transcripts", title = "m_labels")

# --- Plot 2: new_cat_labels ---
GENEcomp2 <- GENEcomparison
GENEcomp2[, "Celltype"] <- just_subset@meta.data$new_cat_labels
GENEcomp2[, "Diagnosis"] <- just_subset@meta.data$diagnosis

GENEcomp2[, "CelltypeDiag"] <- paste(GENEcomp2$Celltype, "/", GENEcomp2$Diagnosis, sep = "")
GENEcomp2$Celltype <- NULL
GENEcomp2$Diagnosis <- NULL

GENEcomp2 <- aggregate(. ~ CelltypeDiag, GENEcomp2, mean)
GENEcomp2 <- BBmisc::normalize(GENEcomp2, method = "range", range = c(0, 10), margin = 2L)
rownames(GENEcomp2) <- GENEcomp2$CelltypeDiag
GENEcomp2$CelltypeDiag <- NULL

GENEcomp2[, "celltypediag"] <- rownames(GENEcomp2)
GENEcomp2[, "Diagnosis"] <- sapply(X = strsplit(GENEcomp2$celltypediag, split = "/"), FUN = "[", 2)
GENEcomp2[, "Celltype"] <- sapply(X = strsplit(GENEcomp2$celltypediag, split = "/"), FUN = "[", 1)
GENEcomp2$celltypediag <- paste(GENEcomp2$Celltype, " (", GENEcomp2$Diagnosis, ")", sep = "")
rownames(GENEcomp2) <- GENEcomp2$celltypediag

p2 <- ggplot(GENEcomp2, aes(x = RNA_Counts, y = Chd1)) +
  geom_point(aes(color = Diagnosis), size = 3) +
  geom_smooth(color = "black", fill = "grey90", linetype = "dashed", method = "lm", se = TRUE) +
  stat_cor(method = "spearman", label.x = 0.5, label.y = 9.5, size = 4) +
  geom_text_repel(
    aes(label = ifelse(RNA_Counts > 6 & Chd1 > 5.5, as.character(Celltype), ""), color = Diagnosis),
    size = 4, point.padding = 0.1, box.padding = 0.5, segment.alpha = 1
  ) +
  theme_classic() +
  theme(
    aspect.ratio = 1,
    axis.text = element_text(size = 12, colour = "black"),
    axis.title = element_text(size = 12)
  ) +
  labs(y = "Chd1", x = "Log2 Transcripts", title = "new_cat_labels")

# Combine and save
combined <- p1 + p2
ggsave(file.path(args$outdir, "panel_G_chd1.pdf"), combined, width = 12, height = 6)
