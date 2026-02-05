library(argparse)
library(Seurat)
library(ggplot2)
library(tidyverse)
library(pheatmap)

gene_set1 <-c(
    "APCDD1",
    "COL18A1",
    "COL23A1",
    "COL13A1",
    "COMP",
    "NKD2",
    "RSPO1",
    "AXIN2",
    "WIF1",
    "SFRP2",
    "CD34",
    "PI16",
    "DPP4",
    "MFAP5",
    "PCOLCE2",
    "CTHRC1",
    "SLPI",
    "CD70",
    "LGR5",
    "CXCL12",
    "APOE",
    "EFEMP1",
    "APOC1",
    "C7",
    "PLA2G2A",
    "PPARG",
    "MYOC",
    "GDF10",
    "CCL19",
    "CD74",
    "CH25H",
    "TNFSF13B",
    "IL33",
    "IRF8",
    "IL15",
    "VCAM1",
    "HLA-DRA",
    "HLA-DRB1",
    "ASPN",
    "COL11A1",
    "MEF2C",
    "DPEP1",
    "MYL4",
    "TNN",
    "COCH",
    "CRABP1",
    "COL24A1",
    "RSPO4",
    "SLITRK6",
    "NRG3",
    "MKX",
    "TNMD",
    "CORIN",
    "BMP7",
    "WNT5A",
    "LEF1",
    "HHIP",
    "RSPO3",
    "INHBA",
    "PTCH1",
    "SCN7A",
    "FMO2",
    "FGFBP2",
    "OLFML2A",
    "PEAR1",
    "RAMP1",
    "RELN",
    "PLEKHA6",
    "IGFBP2",
    "SFRP1",
    "EBF2",
    "NGFR",
    "SFRP4",
    "ITGA6",
    "CDH19",
    "CLDN1"
)

gene_set2 <- c(
    "CCL19",
    "CD74",
    "CH25H",
    "ADAMDEC1",
    "CXCL9",
    "COL18A1",
    "COL23A1",
    "COL13A1",
    "NKD2",
    "CRABP1",
    "TNFRSF21",
    "CXCL1",
    "WNT5A",
    "COL7A1",
    "CXCL8",
    "CXCL5",
    "CXCL6",
    "CXCL13",
    "IL11",
    "IL24",
    "MMP1",
    "MMP3",
    "IL7R",
    "INHBA",
    "CSF3",
    "TDO2",
    "WWC1",
    "CHI3L1",
    "CCL5",
    "CCL11",
    "FAM167A",
    "HIF1A",
    "WNT2",
    "COL10A1",
    "LAMP5",
    "NRG1",
    "OGN",
    "ACTA2",
    "TAGLN",
    "ASPN",
    "COMP",
    "COL11A1",
    "KIF26B",
    "ZNF469",
    "RUNX2",
    "SULF1",
    "ADAM12",
    "ADAM19",
    "COL8A1",
    "LRRC15",
    "CREB3L1",
    "CTHRC1",
    "CCN4",
    "FABP5",
    "CDH2",
    "C1QTNF3",
    "CADM1",
    "LRRC17",
    "PIEZO2",
    "SFRP4",
    "KCNMA1",
    "ACAN",
    "SCX",
    "THBS4",
    "ITGA10",
    "FGF18",
    "PRG4",
    "CRTAC1"
)


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

sanitize_label_for_filename <- function(label) {
    gsub("[/\\:*?\"<>|]", "_", label)
}

sanitize_label_for_AggregateExpression <- function(label) {
    gsub("[/\\:*?\"<>|_]", "-", label)
}

make_plot <- function(seurat_object, gene_set, separator_genes, identity, filenm){
    line_positions <- which(gene_set %in% separator_genes) + 0.5

    plot <- DotPlot(seurat_object, gene_set, group.by = identity, cluster.idents = TRUE, dot.scale = 2) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6)) +
        theme(axis.text.y = element_text(hjust = 1, vjust = 0.5, size = 6)) +
        scale_color_gradient(low = "black", high = "red") +
        geom_vline(xintercept = line_positions, color = "black", linewidth = 0.3)

    ggsave(filenm, plot, width = 10, height = 5)
}

# Create an ArgumentParser object
parser <- ArgumentParser(description = 'Integrate multiple scRNA-seq datasets into one Seurat object')
parser$add_argument("joined_integrated_seurat_object_with_cytotrace_labels", help="File paths to joined_integrated_seurat_object.rds")
parser$add_argument("--cell_type_identity", help="identity to use", required = TRUE)
parser$add_argument("--outdir", help="directory to write to", required = TRUE)
parser$add_argument("--file_suff", help="file suffix", required = TRUE)

# Parse the command-line arguments
args <- parser$parse_args()

print(args)

seurat_object <- readRDS(args$joined_integrated_seurat_object_with_cytotrace_labels)
seurat_object$diagnosis <- trimws(seurat_object$diagnosis)

seurat_object_just_SSC <- subset(seurat_object, diagnosis == 'SSC')
seurat_object_just_control <- subset(seurat_object, diagnosis == 'control')

print(Assays(seurat_object_just_SSC))


dir.create(args$outdir, showWarnings = FALSE, recursive = TRUE)
date_and_time <- format(now(), "%Y%m%d_%H%M%S")

# Define genes after which to add separator lines
separator_genes_plottype_1 <- c("SFRP2", "LGR5", "GDF10", "HLA-DRB1",
"COL11A1", "MYL4", "TNMD", "PTCH1", "PEAR1", "SFRP1")

line_positions_plottype_1 <- which(gene_set1 %in% separator_genes_plottype_1) + 0.5
filenm <- file.path(
    args$outdir,
    paste0("just_SSC_singleR_labels_gene_set1_", args$file_suff, "_", date_and_time, ".pdf")
)
make_plot(seurat_object_just_SSC, gene_set1, separator_genes_plottype_1, "singleR.labels_fine", filenm)

filenm <- file.path(
    args$outdir,
    paste0("just_control_singleR_labels_gene_set1_", args$file_suff, "_", date_and_time, ".pdf")
)
make_plot(seurat_object_just_control, gene_set1, separator_genes_plottype_1, "singleR.labels_fine", filenm)

filenm <- file.path(
    args$outdir,
    paste0("just_SSC_CytoTRACE2_Potency_gene_set1_", args$file_suff, "_", date_and_time, ".pdf")
)
make_plot(seurat_object_just_SSC, gene_set1, separator_genes_plottype_1, "CytoTRACE2_Potency", filenm)

filenm <- file.path(
    args$outdir,
    paste0("just_control_CytoTRACE2_Potency_gene_set1_", args$file_suff, "_", date_and_time, ".pdf")
)
make_plot(seurat_object_just_control, gene_set1, separator_genes_plottype_1, "CytoTRACE2_Potency", filenm)

separator_genes_plottype_2 <- c(
    "CXCL9",
    "CXCL1",
    "COL7A1",
    "HIF1A",
    "LAMP5",
    "OGN",
    "CDH2",
    "KCNMA1",
    "CRTAC1"
)

line_positions_plottype_1 <- which(gene_set2 %in% separator_genes_plottype_2) + 0.5
filenm <- file.path(
    args$outdir,
    paste0("just_SSC_singleR_labels_gene_set2_", args$file_suff, "_", date_and_time, ".pdf")
)
make_plot(seurat_object_just_SSC, gene_set2, separator_genes_plottype_2, "singleR.labels_fine", filenm)

filenm <- file.path(
    args$outdir,
    paste0("just_control_singleR_labels_gene_set2_", args$file_suff, "_", date_and_time, ".pdf")
)
make_plot(seurat_object_just_control, gene_set2, separator_genes_plottype_2, "singleR.labels_fine", filenm)

filenm <- file.path(
    args$outdir,
    paste0("just_SSC_CytoTRACE2_Potency_gene_set2_", args$file_suff, "_", date_and_time, ".pdf")
)
make_plot(seurat_object_just_SSC, gene_set2, separator_genes_plottype_2, "CytoTRACE2_Potency", filenm)

filenm <- file.path(
    args$outdir,
    paste0("just_control_CytoTRACE2_Potency_gene_set2_", args$file_suff, "_", date_and_time, ".pdf")
)
make_plot(seurat_object_just_control, gene_set2, separator_genes_plottype_2, "CytoTRACE2_Potency", filenm)
