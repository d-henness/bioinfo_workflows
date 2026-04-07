library(argparse)
library(Seurat)
library(ggplot2)
library(tidyverse)
library(pheatmap)


mapping <- read.csv("cell_type_mapping.csv")


all_gene_sets <- list(
    gene_set = c(
        "CCN2",
        "COL14A1",
        "COL5A3",
        "BGN",
        "PCOLCE",
        "POSTN",
        "CTHRC1",
        "TNC",
        "COMP",
        "LOX",
        "TIMP1",
        "FAP",
        "THY1",
        "IL6",
        "TGFB1",
        "TGFBI",
        "CXCL1",
        "CXCL8",
        "SAA1",
        "PTX3",
        "CHI3L1",
        "VEGFB",
        "ADM",
        "ID4",
        "HTRA1",
        "ACTG2",
        "ACTA2",
        "TAGLN2",
        "TPM1",
        "PLOD2",
        "LOXL2",
        "NOX4",
        "ADAM12",
        "IL11",
        "LRRC15",
        "RUNX2",
        "KIF26B",
        "MKI67",
        "CCNA2",
        "CCNB1",
        "TOP2A",
        "BIRC5",
        "CENPF",
        "UBE2C",
        "AURKA",
        "TYMS",
        "CDK1",
        "PLK1",
        "FOXM1",
        "VEGFA",
        "ANGPT2",
        "RGS5",
        "NDUFA4L2",
        "CXCL9",
        "CXCL14",
        "PRSS23",
        "COL4A3",
        "COL4A4",
        "MMP11",
        "CEMIP",
        "SERPINE2",
        "CILP2",
        "LTBP2",
        "TNFRSF6B",
        "STMN2",
        "CDKN2A",
        "MDK",
        "SMOC2"
    )
)

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


for (diag in unique(just_subset$diagnosis)){
    just_diag <- subset(just_subset, diagnosis == diag)

    cell_types_level_order <- c(
        "Superficial",
        "Reticular",
        "Adipogenic",
        "Perivascular",
        "Hair_follicle",
        "Mitotic",
        "Stress/Contractile"
    )

    just_diag$new_cat_labels <- factor(just_diag$new_cat_labels, levels = rev(cell_types_level_order))

    for (gs_name in names(all_gene_sets)){
        dir_name <- file.path(args$outdir, gs_name)
        dir.create(dir_name, showWarnings = FALSE, recursive = TRUE)

        filenm <- file.path(
            dir_name,
            paste0("all_datasets_cell_group_", diag, "_", args$file_suff, "_", date_and_time)
        )
        gene_set <- all_gene_sets[[gs_name]]
        make_plot(just_diag, gene_set, NULL, "new_cat_labels", filenm)
    }
}
