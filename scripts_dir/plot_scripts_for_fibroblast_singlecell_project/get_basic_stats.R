library(argparse)
library(Seurat)
library(ggplot2)
library(SingleR)
library(celldex)
library(tidyverse)
library(pheatmap)
library(EnhancedVolcano)

options(warn = 1)

mapping <- read.csv("cell_type_mapping.csv")


joined_integrated_data <- readRDS("rds_files/joined_integrated_seurat_object_with_group_data.rds")
wanted_groups = c("fibroblasts", "pericyte")

just_subset <- subset(joined_integrated_data, cell_group %in% wanted_groups)
just_subset$m_labels <- mapping$m_labels[match(just_subset$singleR.labels_fine, mapping$singleR.labels_fine)]
just_subset$new_cat_labels <- mapping$new_cat_labels[match(just_subset$singleR.labels_fine, mapping$singleR.labels_fine)]

cells_per_sample <- table(joined_integrated_data$dataset, joined_integrated_data$diagnosis)
print(cells_per_sample)

cells_per_sample <- table(just_subset$dataset, just_subset$diagnosis)
print(cells_per_sample)

cells_per_sample <- table(just_subset$dataset, just_subset$new_cat_labels)
print(cells_per_sample)

cells_per_sample <- table(just_subset$diagnosis, just_subset$new_cat_labels)
print(cells_per_sample)
