library(argparse)
library(Seurat)
library(ggplot2)
library(tidyverse)
library(pheatmap)
library(purrr)
library(patchwork)

mapping <- read.csv("cell_type_mapping.csv")

gene_set <- c(
    "ARID5B",
    "EP300",
    "KDM5C",
    "KMT2C",
    "NCOR1",
    "NSD1",
    "SETD2",
    "ARID1A",
    "ARID2",
    "CHD1",
    "CHD8",
    "H3F3C",
    "CTCF",
    "BRCA1",
    "POLQ",
    "SMC1A",
    "CUL3",
    "USP9X",
    "MGA",
    "DHX9",
    "DICER1",
    "XPO1",
    "BRD4",
    "SUPT6H",
    "HEXIM1",
    "LARP7",
    "EZH2",
    "DOT1L",
    "SUPT16H",
    "CDK9",
    "HEXIM1"
)

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

get_ave_express <- function(df, identity){
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
    df_out <- df_out[gene_set, ]

    df_out
}

parser <- ArgumentParser(description = 'Integrate multiple scRNA-seq datasets into one Seurat cbject')
parser$add_argument("--outdir", help="directory to write to", required = TRUE)
args <- parser$parse_args()
print(args)
dir.create(args$outdir, showWarnings = FALSE, recursive = TRUE)

joined_integrated_data <- readRDS("rds_files/cytotrace_seurat.rds")
wanted_groups = c("fibroblasts", "pericyte")

just_subset <- subset(joined_integrated_data, cell_group %in% wanted_groups)
just_subset <- subset(just_subset, !is.na(CytoTRACE2_Score))

just_subset$m_labels <- mapping$m_labels[match(just_subset$singleR.labels_fine, mapping$singleR.labels_fine)]
just_subset$new_cat_labels <- mapping$new_cat_labels[match(just_subset$singleR.labels_fine, mapping$singleR.labels_fine)]

for (diag in unique(just_subset$diagnosis)){
    just_diag <- subset(just_subset, diagnosis == diag)
    filenm1 <- file.path(args$outdir, paste0("average_expression_m_labels_", diag, ".csv"))
    ave_express <- get_ave_express(just_diag, "m_labels")

    just_diag@assays$RNA@layers$abs_scale <- as.matrix(log2(LayerData(just_diag, assay = "RNA", layer = "counts") + 1))
    just_diag@meta.data$scaledcounts <- colSums(just_diag@assays$RNA@layers$abs_scale)

    filenm_genes <- file.path(args$outdir, paste0("heatmap_ave_express_", sanitize_label_for_filename(diag), ".pdf"))
    pdf(filenm_genes)
    pheatmap(
        ave_express,
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        scale        = "none",
        main         = paste0("Average Expression (", diag, ")"),
        fontsize_row = 8,
        fontsize_col = 8
    )
    dev.off()

    meta_avgs <- just_diag@meta.data %>%
        group_by(m_labels) %>%
        summarise(
            CytoTRACE2_Score = mean(CytoTRACE2_Score, na.rm = TRUE),
            scaledcounts     = mean(scaledcounts, na.rm = TRUE),
            .groups = "drop"
        )

    print(meta_avgs)

    cytotrace_row  <- setNames(meta_avgs$CytoTRACE2_Score, meta_avgs$m_labels)
    scaledcounts_row <- setNames(meta_avgs$scaledcounts,   meta_avgs$m_labels)

    ave_express["CytoTRACE2_Score", ] <- cytotrace_row[colnames(ave_express)]
    ave_express["scaledcounts", ]     <- scaledcounts_row[colnames(ave_express)]
    filenm <- file.path(args$outdir, paste0("table_ave_express_and_meta_", sanitize_label_for_filename(diag), ".csv"))
    write.table(ave_express, filenm1, sep = ",", col.names = NA)
}
