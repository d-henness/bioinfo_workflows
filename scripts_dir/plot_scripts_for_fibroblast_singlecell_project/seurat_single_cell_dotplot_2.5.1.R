library(argparse)
library(Seurat)
library(ggplot2)
library(tidyverse)
library(pheatmap)


gene_set <- c(
    "APCDD1",
    "COL23A1",
    "NKD2",
    "WIF1",
    "SFRP2",
    "PI16",
    "CD34",
    "DPP4",
    "MFAP5",
    "PPARG",
    "FABP4",
    "LPL",
    "CIDEA",
    "CCL19",
    "CD74",
    "HLA-DRA",
    "VCAM1",
    "ASPN",
    "COL11A1",
    "CILP2",
    "TNMD",
    "COCH",
    "SCN7A",
    "FMO2",
    "NGFR",
    "RAMP1",
    "CALCA",
    "ACTA2",
    "LRRC15",
    "CTHRC1",
    "LOX",
    "FAP",
    "ADAM12",
    "IL11",
    "CXCL8",
    "CXCL13",
    "IL7R",
    "PTX3",
    "CCL20",
    "MKI67",
    "TOP2A",
    "CDK1",
    "AURKB",
    "PLK1",
    "HIF1A",
    "DDIT4",
    "NDRG2",
    "HSPA6",
    "FBXO32",
    "ISG20"
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


    plot <- DotPlot(seurat_object, gene_set, group.by = identity, dot.scale = 2) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6)) +
        theme(axis.text.y = element_text(hjust = 1, vjust = 0.5, size = 6)) +
        scale_color_gradient(low = "black", high = "red") #+
    #   uncomment if you want the lines back
    #    geom_vline(xintercept = line_positions, color = "black", linewidth = 0.3)

    print(5 * (nrow(unique(seurat_object[[identity]])) / 10))
    ggsave(filenm, plot, width = 2 * length(gene_set) / 10, height = 2 * (nrow(unique(seurat_object[[identity]])) / 10), limitsize = FALSE)
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

cell_types <- c(
  "Fibroblasts:breast",
  "Fibroblasts:foreskin",
  "Tissue_stem_cells:BM_MSC",
  "MSC",
  "Tissue_stem_cells:iliac_MSC",
  "Chondrocytes:MSC-derived",
  "Tissue_stem_cells:BM_MSC:BMP2",
  "Tissue_stem_cells:BM_MSC:TGFb3",
  "Osteoblasts",
  "Osteoblasts:BMP2",
  "Smooth_muscle_cells:bronchial",
  "Smooth_muscle_cells:bronchial:vit_D",
  "Smooth_muscle_cells:vascular:IL-17",
  "iPS_cells:adipose_stem_cells",
  "iPS_cells:CRL2097_foreskin",
  "iPS_cells:fibroblasts",
  "iPS_cells:PDB_fibroblasts",
  "iPS_cells:skin_fibroblast",
  "Tissue_stem_cells:adipose-derived_MSC_AM3"
)

diagnoses <- c("normal_skin", "scar", "keloid", "ssc")

# expand.grid varies the second arg fastest by default,
# but we want diagnosis to cycle within each cell type:
ordered_levels <- rev(paste(
  rep(cell_types, each = length(diagnoses)),
  rep(diagnoses, times = length(cell_types)),
  sep = "_"
))

just_subset$singleR.labels_fine_plus_diagnosis <- factor(
  paste(just_subset$singleR.labels_fine, just_subset$diagnosis, sep = "_"),
  levels = ordered_levels
)

filenm <- file.path(
    args$outdir,
    paste0("all_datasets_cell_group_", args$file_suff, "_", date_and_time, ".pdf")
)
make_plot(just_subset, gene_set, NULL, "singleR.labels_fine_plus_diagnosis", filenm)
