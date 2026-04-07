library(argparse)
library(Seurat)
library(ggplot2)
library(tidyverse)
library(pheatmap)


gene_set <- c(
    "COL23A1",
    "COL13A1",
    "COL7A1",
    "COL10A1",
    "COL11A1",
    "COL8A1",
    "COMP",
    "OGN",
    "ACAN",
    "PRG4",
    "CRTAC1",
    "THBS4",
    "ZNF469",
    "ITGA10",
    "P4HA3",
    "HAPLN1",
    "LUM",
    "DCN",
    "FBLN1",
    "FBLN2",
    "CILP",
    "OMD",
    "ECM1",
    "SPON2",
    "CHPF",
    "LOX",
    "SERPINE2",
    "SERPINF1",
    "TIMP1",
    "MMP2",
    "ADAMTS6",
    "ADAMTSL1",
    "PODN",
    "PODNL1",
    "ABI3BP",
    "KERA",
    "ELN",
    "MATN4",
    "ACTA2",
    "WNT2",
    "WNT5A",
    "INHBA",
    "IL11",
    "ADAM12",
    "ADAM19",
    "CREB3L1",
    "CTHRC1",
    "CCN4",
    "SFRP4",
    "LRRC15",
    "RUNX2",
    "FGF18",
    "NKD2",
    "SCX",
    "KIF26B",
    "NOX4",
    "FAP",
    "PDPN",
    "CCN3",
    "IL13RA2",
    "FGF9",
    "GDF11",
    "SOSTDC1",
    "CHRDL1",
    "BMP5",
    "LRG1",
    "DHH",
    "WNT6",
    "RARRES1",
    "CADM1",
    "KCNMA1",
    "WWC1",
    "NRG1",
    "CRABP1",
    "MYH10",
    "FLNA",
    "CAP2",
    "CFL2",
    "DNMBP",
    "PLEKHG3",
    "ARHGEF26",
    "ARAP3",
    "ARHGEF15",
    "RUFY1",
    "MICAL3",
    "PAK1",
    "RGS3",
    "RALGDS",
    "HMMR",
    "RCSD1",
    "CCL19",
    "CXCL1",
    "CXCL6",
    "CXCL8",
    "CXCL13",
    "CSF3",
    "IL24",
    "TNFRSF21",
    "CHI3L1",
    "TDO2",
    "IL7R",
    "CXCL2",
    "CXCL3",
    "CXCL14",
    "CCL20",
    "IL32",
    "PTX3",
    "TNFAIP6",
    "ACKR2",
    "CPVL",
    "MSR1",
    "TAB1",
    "NFATC3",
    "NDRG2",
    "DDIT4",
    "SBNO2",
    "HIF1A",
    "FABP5",
    "C1QTNF3",
    "AOX1",
    "HSD11B1",
    "SIRT6",
    "SREBF1",
    "ST3GAL2",
    "SLC4A2",
    "GPX3",
    "RBP4",
    "FABP4",
    "CD36",
    "APOD",
    "C3",
    "CFB",
    "CFH",
    "CLU",
    "AOC3",
    "COLEC11",
    "ACSL1",
    "MGST1",
    "STEAP4",
    "INPP4B",
    "MKI67",
    "CCNA2",
    "CCNB1",
    "CCNB2",
    "CDK1",
    "PLK1",
    "AURKB",
    "FOXM1",
    "TOP2A",
    "BIRC5",
    "CDC20",
    "BUB1",
    "BUB1B",
    "TTK",
    "CENPF",
    "CEP55",
    "CKAP2",
    "CKAP2L",
    "NUF2",
    "TPX2",
    "PRC1",
    "RRM2",
    "SGO1",
    "SKA1",
    "GTSE1",
    "DLGAP5",
    "KIF11",
    "KIF14",
    "KIF15",
    "KIF18A",
    "KIF18B",
    "KIF20A",
    "KIF23",
    "KIF4A",
    "KIFC1",
    "MELK",
    "NEK2",
    "UBE2C",
    "UBE2T",
    "ANLN",
    "ASPM",
    "HJURP",
    "NCAPG",
    "TACC3",
    "SPC25"
)



#cell_order <- c(
#    "B-cells",
#    "mast cells",
#    "neutrophils",
#    "dendritic cells",
#    "macrophage",
#    "T-cells",
#    "neurons",
#    "smooth muscle",
#    "pericyte",
#    "fibroblasts",
#    "endothelial",
#    "eccrine gland",
#    "melanocytes",
#    "keratinocytes"
#)


sanitize_label_for_filename <- function(label) {
    gsub("[/\\:*?\"<>|]", "_", label)
}

sanitize_label_for_AggregateExpression <- function(label) {
    gsub("[/\\:*?\"<>|_]just_control_CytoTRACE2_Potency_gene_set2_", "-", label)
}

make_plot <- function(seurat_object, gene_set, separator_genes, identity, filenm){
    #line_positions <- which(gene_set %in% separator_genes) + 0.5


    plot <- DotPlot(seurat_object, gene_set, group.by = identity, dot.scale = 2, scale = FALSE, scale.min = 0, scale.max = 60) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6)) +
        theme(axis.text.y = element_text(hjust = 1, vjust = 0.5, size = 6)) +
        scale_color_gradient(low = "black", high = "red") #+
    #   uncomment if you want the lines back
    #    geom_vline(xintercept = line_positions, color = "black", linewidth = 0.3)

    print(5 * (nrow(unique(seurat_object[[identity]])) / 10))
    ggsave(filenm, plot, width = 2 * length(gene_set) / 10, height = 2 * (nrow(unique(seurat_object[[identity]])) / 10), limitsize = FALSE)
    write.csv(plot$data, paste0(tools::file_path_sans_ext(filenm), ".csv"))
}

# Create an ArgumentParser object
parser <- ArgumentParser(description = 'Integrate multiple scRNA-seq datasets into one Seurat object')
parser$add_argument("joined_integrated_seurat_object_with_cytotrace_labels", help="File paths to joined_integrated_seurat_object.rds")
#parser$add_argument("--cell_type_identity", help="identity to use", required = TRUE)
#parser$add_argument("--group_file", help="file to group cell types", required = TRUE)
parser$add_argument("--outdir", help="directory to write to", required = TRUE)
parser$add_argument("--file_suff", help="file suffix", required = TRUE)

# Parse the command-line arguments
args <- parser$parse_args()

print(args)

seurat_object <- readRDS(args$joined_integrated_seurat_object_with_cytotrace_labels)
seurat_object$diagnosis <- trimws(seurat_object$diagnosis)
print(head(seurat_object))

dir.create(args$outdir, showWarnings = FALSE, recursive = TRUE)
date_and_time <- format(now(), "%Y%m%d_%H%M%S")

wanted_groups = c("fibroblasts", "pericyte")
just_subset <- subset(seurat_object, cell_group %in% wanted_groups)
missing <- gene_set[!gene_set %in% rownames(just_subset)]
print("here")
print(missing)
print("here")

just_subset$singleR.labels_fine_plus_diagnosis <- factor(
  paste(just_subset$singleR.labels_fine, just_subset$diagnosis, sep = "_"),
  levels = sort(unique(paste(just_subset$singleR.labels_fine, just_subset$diagnosis, sep = "_")))
)

filenm <- file.path(
    args$outdir,
    paste0("all_datasets_cell_group_", args$file_suff, "_", date_and_time, ".pdf")
)
make_plot(just_subset, gene_set, NULL, "singleR.labels_fine_plus_diagnosis", filenm)
