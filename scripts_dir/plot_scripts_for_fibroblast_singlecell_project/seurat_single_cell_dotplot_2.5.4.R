library(argparse)
library(Seurat)
library(ggplot2)
library(tidyverse)
library(pheatmap)


gene_set <- c( 
    "SFRP2",
    "APCDD1",
    "COL23A1",
    "NKD2",
    "WIF1",
    "PI16",
    "CD34",
    "MFAP5",
    "PPARG",
    "FABP4",
    "LPL",
    "NDRG2",
    "CXCL1",
    "CCL19",
    "CD74",
    "HLA-DRA",
    "VCAM1",
    "CXCL8",
    "PTX3",
    "ASPN",
    "COL11A1",
    "TNMD",
    "COCH",
    "RAMP1",
    "PLP1",
    "MKI67",
    "TOP2A",
    "CDK1",
    "AURKB",
    "MELK",
    "TROAP",
    "DDIT4",
    "HIF1A",
    "ACTA2",
    "LOX"
)

cell_type_mapping <- c(
    "Fibroblasts:breast" = "Superficial APCDD1",
    "iPS_cells:fibroblasts" = "Superficial APCDD1",
    "iPS_cells:PDB_fibroblasts" = "Superficial APCDD1",
    "iPS_cells:skin_fibroblast" = "Superficial APCDD1",
    "iPS_cells:adipose_stem_cells" = "Superficial APCDD1",
    "Chondrocytes:MSC-derived" = "Reticular_PI16",
    "Smooth_muscle_cells:bronchial" = "Reticular_PI16",
    "Smooth_muscle_cells:bronchial:vit_D" = "Reticular_PI16",
    "Tissue_stem_cells:iliac_MSC" = "Adipogenic_FABP4",
    "Fibroblasts:foreskin" = "Perivascular",
    "Smooth_muscle_cells:vascular:IL-17" = "Perivascular",
    "Osteoblasts" = "Hair_follicle_like",
    "Osteoblasts:BMP2" = "Hair_follicle_like",
    "iPS_cells:CRL2097_foreskin" = "Schwann-like",
    "MSC" = "Mitotic",
    "Tissue_stem_cells:adipose-derived_MSC_AM3" = "Mitotic",
    "Tissue_stem_cells:BM_MSC:BMP2" = "Stressed/Contractile",
    "Tissue_stem_cells:BM_MSC" = "Stressed/Contractile",
    "Tissue_stem_cells:BM_MSC:TGFb3" = "Stressed/Contractile"
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


    plot <- DotPlot(seurat_object, gene_set, group.by = identity, col.max = 4, dot.scale = 8, scale = FALSE, scale.min = 0, scale.max = 100) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 14),
        axis.text.y = element_text(hjust = 1, vjust = 0.5, size = 14), axis.title.x = element_blank(), axis.title.y = element_blank()) +
        scale_color_gradient(low = "#FF9999", high = "darkred", limits = c(0, 4)) #+
    #   uncomment if you want the lines back
#        geom_vline(xintercept = line_positions, color = "black", linewidth = 0.3)

    print(max(3.5, 2 * (nrow(unique(seurat_object[[identity]])) / 10)))
    ggsave(paste0(filenm, ".pdf"), plot, width = 3.5 * length(gene_set) / 10, height = max(3.5, 2 * (nrow(unique(seurat_object[[identity]])) / 10)), limitsize = FALSE)
    write.csv(plot$data, paste0(filenm, ".csv"))
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
for (diag in unique(seurat_object$diagnosis)){
    just_diag <- subset(just_subset, diagnosis == diag)

    print("here")
    just_diag$new_cat_labels <- unname(cell_type_mapping[just_diag$singleR.labels_fine])

    missing <- gene_set[!gene_set %in% rownames(just_diag)]
    print("here")
    print(missing)
    print("here")

    cell_types <- c(
        "Superficial APCDD1",
        "Reticular_PI16",
        "Adipogenic_FABP4",
        "Perivascular",
        "Hair_follicle_like",
        "Schwann-like",
        "Mitotic",
        "Stressed/Contractile"
    )

    just_diag$new_cat_labels <- factor(just_diag$new_cat_labels, levels = rev(cell_types))

    filenm <- file.path(
        args$outdir,
        paste0("all_datasets_cell_group_", diag, "_", args$file_suff, "_", date_and_time)
    )
    make_plot(just_diag, gene_set, separator_genes, "new_cat_labels", filenm)
}
