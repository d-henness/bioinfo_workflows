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
    "PLIN1",
    "THRSP",
    "AWAT2",
    "CCL19",
    "CD74",
    "HLA-DRA",
    "VCAM1",
    "IL11",
    "CXCL8",
    "CXCL13",
    "IL7R",
    "PTX3",
    "CCL20",
    "CXCL1",
    "SELE",
    "SERPINB2",
    "MMP1",
    "ASPN",
    "COL11A1",
    "CILP2",
    "TNMD",
    "COCH",
    "HAPLN1",
    "HAS1",
    "ADAMTS18",
    "SENCR",
    "SCN7A",
    "FMO2",
    "NGFR",
    "RAMP1",
    "CALCA",
    "SOX10",
    "RELN",
    "PLP1",
    "GJA4",
    "MKI67",
    "TOP2A",
    "CDK1",
    "AURKB",
    "PLK1",
    "MELK",
    "TROAP",
    "TERT",
    "HIF1A",
    "DDIT4",
    "NDRG2",
    "HSPA6",
    "FBXO32",
    "ISG20",
    "TM7SF3",
    "TK2",
    "PLVAP",
    "SOX17",
    "ACTA2",
    "LRRC15",
    "CTHRC1",
    "LOX",
    "FAP",
    "ADAM12"
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
#    line_positions <- which(gene_set %in% separator_genes) + 0.5


    plot <- DotPlot(seurat_object, gene_set, group.by = identity, dot.scale = 2, scale = FALSE, scale.min = 0, scale.max = 60) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6)) +
        theme(axis.text.y = element_text(hjust = 1, vjust = 0.5, size = 6)) +
        scale_color_gradient(low = "black", high = "red") #+
    #   uncomment if you want the lines back
#        geom_vline(xintercept = line_positions, color = "black", linewidth = 0.3)

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
just_subset <- subset(just_subset, diagnosis == "normal_skin")
missing <- gene_set[!gene_set %in% rownames(just_subset)]
print("here")
print(missing)
print("here")

cell_types <- c(
    "Fibroblasts:breast",
    "iPS_cells:fibroblasts",
    "iPS_cells:PDB_fibroblasts",
    "iPS_cells:skin_fibroblast",
    #
    "Chondrocytes:MSC-derived",
    "Smooth_muscle_cells:bronchial",
    "Smooth_muscle_cells:bronchial:vit_D",
    #
    "Tissue_stem_cells:iliac_MSC",
    #
    "Fibroblasts:foreskin",
    "Smooth_muscle_cells:vascular:IL-17",
    #
    "Osteoblasts",
    "Osteoblasts:BMP2",
    "iPS_cells:CRL2097_foreskin",
    #
    "MSC",
    "Tissue_stem_cells:adipose-derived_MSC_AM3",
    #
    "Tissue_stem_cells:BM_MSC:BMP2",
    "Tissue_stem_cells:BM_MSC",
    "Tissue_stem_cells:BM_MSC:TGFb3",
    "iPS_cells:adipose_stem_cells"
)

diagnoses <- c("normal_skin")#, "scar", "keloid", "ssc")

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

separator_genes <- c(
    # F2 — Universal "PI16",
    "CD34",
    # AdipoMSC — Lipogenic "PPARG",
    "FABP4",
    # F3 — Immunomodulatory "CCL19",
    "CD74",
    # D1 — Inflammatory "IL11",
    "CXCL8",
    # F4 — Niche / Dermal Papilla "ASPN",
    "COL11A1",
    # F5 — Schwann / Neural Crest "SCN7A",
    "FMO2",
    # D3 — Proliferating "MKI67",
    "TOP2A",
    # Stress / UPR "HIF1A",
    "DDIT4",
    # MEndoT — Mesenchymal-to-Endothelial "PLVAP",
    "SOX17",
    # D2 — Contractile / Myofibroblast "ACTA2",
    "LRRC15"
)

filenm <- file.path(
    args$outdir,
    paste0("all_datasets_cell_group_", args$file_suff, "_", date_and_time, ".pdf")
)
make_plot(just_subset, gene_set, separator_genes, "singleR.labels_fine_plus_diagnosis", filenm)
