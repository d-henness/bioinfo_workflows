library(argparse)
library(Seurat)
library(ggplot2)
library(tidyverse)
library(pheatmap)


mapping <- read.csv("cell_type_mapping.csv")


all_gene_sets <- list(
    # ECM & collagen (n=70) 
    ECM_and_collagen = c(
        "ADAM33",
        "ADAMTS2",
        "ADAMTS7",
        "ADAMTSL1",
        "ADAMTSL4",
        "ARSB",
        "BGN",
        "CCN3",
        "CCN4",
        "CCN5",
        "CHST12",
        "CHST14",
        "CILP",
        "CILP2",
        "COL10A1",
        "COL16A1",
        "COL24A1",
        "COL25A1",
        "COL5A1",
        "COL6A6",
        "COMP",
        "CTHRC1",
        "DCN",
        "ECM1",
        "ECM2",
        "EFEMP1",
        "EFEMP2",
        "ELN",
        "EMILIN1",
        "EMILIN2",
        "FBLN1",
        "FBLN2",
        "FBLN5",
        "FBLN7",
        "FBN1",
        "HAPLN1",
        "HSPG2",
        "HYAL2",
        "ITIH5",
        "LOX",
        "LOXL3",
        "LOXL4",
        "LTBP3",
        "MAMDC2",
        "MATN2",
        "MFAP3",
        "MFAP3L",
        "MMP11",
        "MMP23B",
        "MMP27",
        "MXRA5",
        "MXRA8",
        "NTN1",
        "OGN",
        "P3H1",
        "P3H3",
        "P3H4",
        "P4HA2",
        "P4HTM",
        "PAPLN",
        "PLOD3",
        "PODN",
        "SMOC1",
        "SMOC2",
        "SPON2",
        "SVEP1",
        "THBS2",
        "THBS3",
        "THBS4",
        "TNXB"
    ),

    #Wnt / BMP signalling (n=39)
    Wnt_BMP_signalling = c(
        "ACVRL1",
        "BMP5",
        "BMP6",
        "CHRD",
        "CHRDL1",
        "EVC2",
        "FOXF2",
        "FST",
        "FSTL1",
        "FSTL3",
        "FZD6",
        "GREM1",
        "GREM2",
        "INHBB",
        "LRP4",
        "NBL1",
        "NKD1",
        "PORCN",
        "PRDM5",
        "RSPO1",
        "RSPO3",
        "SFRP1",
        "SFRP2",
        "SFRP4",
        "SLCO1C1",
        "SMAD6",
        "SMAD9",
        "SOSTDC1",
        "SP5",
        "TBX1",
        "TBX5",
        "TGFBR1",
        "TGFBR3",
        "WIF1",
        "WNT10A",
        "WNT16",
        "WNT5B",
        "WNT6",
        "WNT9A"
    ),

    #Adipogenic / lipid metabolism (n=58)
    Adipogenic_lipid_metabolism = c(
        "ACAA1",
        "ACOX1",
        "ACOX2",
        "ACOX3",
        "ACSBG1",
        "ACSL1",
        "ANGPT1",
        "ANGPT4",
        "APMAP",
        "APOC1",
        "AQP7",
        "AWAT2",
        "CD36",
        "CERS2",
        "CHI3L1",
        "CIDEA",
        "CNR1",
        "DECR2",
        "DEGS2",
        "DGAT2",
        "DHCR24",
        "DHCR7",
        "EBF2",
        "ECHS1",
        "ELOVL3",
        "FABP3",
        "FABP4",
        "FADS1",
        "FADS2",
        "FASN",
        "FBP1",
        "GPD1",
        "HMGCLL1",
        "HSD17B13",
        "HSD17B14",
        "LCAT",
        "LPL",
        "LSS",
        "MLYCD",
        "MRAP",
        "MSMO1",
        "PCK1",
        "PECR",
        "PLIN1",
        "PLIN4",
        "PLIN5",
        "PPARG",
        "SCD",
        "SLCO1C1",
        "SPTLC1",
        "SPTLC3",
        "SRD5A1",
        "STEAP4",
        "SUCNR1",
        "THRSP",
        "TM7SF2",
        "TMEM97",
        "VIPR1"
    ),

    #Progenitor / stem cell identity (n=33)
    Progenitor_stem_cell_identity = c(
        "ALDH1A1",
        "ALDH1B1",
        "ALDH3A1",
        "ALDH3A2",
        "ALDH3B1",
        "ALDH9A1",
        "ANGPTL1",
        "ANGPTL2",
        "ANGPTL5",
        "ANGPTL7",
        "CRABP1",
        "DPP4",
        "FGF10",
        "FGF12",
        "FGF13",
        "FGF9",
        "GAS1",
        "GAS6",
        "GDF10",
        "HAND2",
        "IGF1",
        "MEIS1",
        "NTN1",
        "PCOLCE2",
        "PDPN",
        "PRDM6",
        "PRG4",
        "SEMA3B",
        "SEMA3C",
        "SEMA3G",
        "TWIST2",
        "VEGFC",
        "VEGFD"
    ),

    #ER stress / UPR (n=39)
    ER_stress_UPR = c(
        "CREB3",
        "CREB5",
        "CRELD2",
        "DERL2",
        "DNAJB11",
        "DNAJB12",
        "DNAJB2",
        "DNAJC13",
        "DNAJC16",
        "DNAJC25",
        "DNAJC4",
        "DNAJC5",
        "EIF4EBP1",
        "EIF4EBP3",
        "EMC1",
        "ERLIN1",
        "ERLIN2",
        "FBXO25",
        "FBXO4",
        "FKBP10",
        "FKBP7",
        "FKBP9",
        "GANAB",
        "HSPA12A",
        "HSPA12B",
        "HYOU1",
        "MAN1B1",
        "PDIA4",
        "PDIA5",
        "RNF121",
        "RNF135",
        "RNF170",
        "SDF4",
        "SYVN1",
        "TBC1D20",
        "TM7SF3",
        "TMEM208",
        "UGGT1",
        "VAPB"
    ),

    #Mitochondrial stress (n=49)
    Mitochondrial_stress = c(
        "BCS1L",
        "BNIP1",
        "BNIP3",
        "COQ2",
        "COQ5",
        "COQ6",
        "DELE1",
        "ERAL1",
        "ETFA",
        "FASTKD3",
        "GLT8D2",
        "HARS",
        "HARS2",
        "HSDL2",
        "LARS2",
        "MAVS",
        "MCAT",
        "MECR",
        "MRPL15",
        "MRPL17",
        "MRPL2",
        "MRPL38",
        "MRPL46",
        "MRPS14",
        "MRPS18A",
        "MRPS22",
        "MRPS30",
        "MRRF",
        "MTERF1",
        "MTX1",
        "MTX3",
        "NAGK",
        "NDUFA4L2",
        "NDUFAF1",
        "NDUFAF3",
        "NDUFAF5",
        "NSUN3",
        "PINK1",
        "SDHA",
        "SQOR",
        "SVBP",
        "TACO1",
        "TERF1",
        "TIMM22",
        "TIMM23",
        "TK2",
        "TOMM34",
        "TRMT2A",
        "TRMT2B"
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
