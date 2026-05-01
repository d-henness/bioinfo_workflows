library(argparse)
library(Seurat)
library(ggplot2)
library(tidyverse)
library(pheatmap)


mapping <- read.csv("cell_type_mapping.csv")


all_gene_sets <- list(
    gene_set = c(
        "WIF1",
        "SFRP2",
        "GREM2",
        "APCDD1",
        "LEPR",
        "F13A1",
        "MAMDC2",
        "MXRA5",
        "DPP4",
        "GAS1",
        "ELN",
        "PI16",
        "OGN",
        "IGF1",
        "CTHRC1",
        "ABCA10",
        "ASPN",
        "ACTA2",
        "TAGLN",
        "MYH11",
        "MYL9",
        "RGS5",
        "ANGPT2",
        "FABP4",
        "CD36",
        "IL6",
        "CXCL1",
        "CXCL2",
        "CXCL8",
        "CD74",
        "COL11A1",
        "DPEP1",
        "KRT5",
        "DSP",
        "NDRG2",
        "TGFB1",
        "MKI67",
        "TOP2A",
        "CDK1",
        "TK1",
        "CENPK",
        "DIAPH3",
        "PCLAF"
    )
)

sep_genes <- c(
    "ELN",
    "ASPN",
    "ANGPT2",
    "CD36",
    "CD74",
    "TGFB1"
)

sep_cells <- c(
    "M17",
    "M7",
    "M14",
    "M8",
    "M15",
    "M16"
)

cell_types_level_order <- c(
    "M3",
    "M5",
    "M11",
    "M12",
    "M17",
    "M2",
    "M6",
    "M7",
    "M1",
    "M4",
    "M10",
    "M14",
    "M8",
    "M13",
    "M15",
    "M9",
    "M16",
    "M18",
    "M19"
)

sanitize_label_for_filename <- function(label) {
    gsub("[/\\:*?\"<>|]", "_", label)
}

sanitize_label_for_AggregateExpression <- function(label) {
    gsub("[/\\:*?\"<>|_]just_control_CytoTRACE2_Potency_gene_set2_", "-", label)
}

make_plot <- function(seurat_object, gene_set, separator_genes, separator_cells, identity, filenm){
    vline_positions <- if (!is.null(separator_genes)) which(gene_set %in% separator_genes) + 0.5 else c()
    identity_levels <- levels(seurat_object[[identity]][[identity]])
    hline_positions <- if (!is.null(separator_cells)) which(identity_levels %in% separator_cells) + 0.5 else c()

    plot <- DotPlot(seurat_object, gene_set, group.by = identity, col.max = 4, dot.scale = 8, scale = FALSE, scale.min = 0, scale.max = 100) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 14),
        axis.text.y = element_text(hjust = 1, vjust = 0.5, size = 14), axis.title.x = element_blank(), axis.title.y = element_blank()) +
        scale_color_gradient(low = "#FF9999", high = "darkred", limits = c(0, 4), oob = scales::squish) +
        geom_vline(xintercept = vline_positions, color = "black", linewidth = 0.3) +
        geom_hline(yintercept = hline_positions, color = "black", linewidth = 0.3)

    height <- max(3.5, 3 * (nrow(unique(seurat_object[[identity]])) / 10))
    print(height)
    ggsave(paste0(filenm, ".pdf"), plot, width = 3.5 * length(gene_set) / 10, height = height, limitsize = FALSE)
    write.csv(plot$data, paste0(filenm, ".csv"))
}

# Create an ArgumentParser object
parser <- ArgumentParser(description = 'Integrate multiple scRNA-seq datasets into one Seurat object')
parser$add_argument("--outdir", help="directory to write to", required = TRUE)
parser$add_argument("--file_suff", help="file suffix", required = TRUE)

# Parse the command-line arguments
args <- parser$parse_args()

print(args)

dir.create(args$outdir, showWarnings = FALSE, recursive = TRUE)
date_and_time <- format(now(), "%Y%m%d_%H%M%S")

joined_integrated_data <- readRDS("rds_files/joined_integrated_seurat_object_with_group_data.rds")
wanted_groups = c("fibroblasts", "pericyte")

just_subset <- subset(joined_integrated_data, cell_group %in% wanted_groups)
just_subset$m_labels <- mapping$m_labels[match(just_subset$singleR.labels_fine, mapping$singleR.labels_fine)]
just_subset$new_cat_labels <- mapping$new_cat_labels[match(just_subset$singleR.labels_fine, mapping$singleR.labels_fine)]

just_subset$m_labels <- factor(just_subset$m_labels, levels = rev(cell_types_level_order))

for (gs_name in names(all_gene_sets)){
    dir_name <- file.path(args$outdir, gs_name)
    dir.create(dir_name, showWarnings = FALSE, recursive = TRUE)

    filenm <- file.path(
        dir_name,
        paste0("all_datasets_cell_group_all_diag_", args$file_suff, "_", date_and_time)
    )
    gene_set <- all_gene_sets[[gs_name]]
    make_plot(just_subset, gene_set, sep_genes, sep_cells, "m_labels", filenm)
}

for (diag in unique(just_subset$diagnosis)){
    just_diag <- subset(just_subset, diagnosis == diag)

    for (gs_name in names(all_gene_sets)){
        dir_name <- file.path(args$outdir, gs_name)
        dir.create(dir_name, showWarnings = FALSE, recursive = TRUE)

        filenm <- file.path(
            dir_name,
            paste0("all_datasets_cell_group_", diag, "_", args$file_suff, "_", date_and_time)
        )
        gene_set <- all_gene_sets[[gs_name]]
        make_plot(just_diag, gene_set, sep_genes, sep_cells, "m_labels", filenm)
    }
}
