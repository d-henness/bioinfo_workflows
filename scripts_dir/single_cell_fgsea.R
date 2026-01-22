library(argparse)
library(Seurat)
library(ggplot2)
library(SingleR)
library(celldex)
library(tidyverse)
library(pheatmap)
library(msigdbr)
library(fgsea)


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

# hard code for now, do properly later
diagnosis = c(
    "SC1_cellranger" = "control",
    "SC4_cellranger" = "control",
    "SC18_cellranger" = "control",
    "SC32_cellranger" = "control",
    "SC33_cellranger" = "control",
    "SC34_cellranger" = "control",
    "SC50_cellranger" = "control",
    "SC68_cellranger" = "control",
    "SC124_cellranger" = "control",
    "SC125_cellranger" = "control",
    "SC2_cellranger" = "SSC",
    "SC5_cellranger" = "SSC",
    "SC19_cellranger" = "SSC",
    "SC49_cellranger" = "SSC",
    "SC60_cellranger" = "SSC",
    "SC69_cellranger" = "SSC",
    "SC70_cellranger" = "SSC",
    "SC86_cellranger" = "SSC",
    "SC119_cellranger" = "SSC",
    "SC185_cellranger" = "SSC",
    "SC188_cellranger" = "SSC",
    "SC189_cellranger" = "SSC"
)

target_sets <- c(
    "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
    "HALLMARK_DNA_REPAIR",
    "REACTOME_CELLULAR_SENESCENCE",
    "REACTOME_DNA_DAMAGE_TELOMERE_STRESS_INDUCED_SENESCENCE",
    "REACTOME_ONCOGENE_INDUCED_SENESCENCE",
    "REACTOME_OXIDATIVE_STRESS_INDUCED_SENESCENCE",
    "REACTOME_SENESCENCE_ASSOCIATED_SECRETORY_PHENOTYPE_SASP",
    "WP_SENESCENCEASSOCIATED_SECRETORY_PHENOTYPE_SASP",
    "GOBP_EXIT_FROM_MITOSIS",
    "GOBP_MICROTUBULE_CYTOSKELETON_ORGANIZATION_INVOLVED_IN_MITOSIS",
    "GOBP_POSITIVE_REGULATION_OF_EXIT_FROM_MITOSIS",
    "GOBP_REGULATION_OF_EXIT_FROM_MITOSIS",
    "HERNANDEZ_ABERRANT_MITOSIS_BY_DOCETACEL_2NM_DN",
    "HERNANDEZ_ABERRANT_MITOSIS_BY_DOCETACEL_2NM_UP",
    "HERNANDEZ_ABERRANT_MITOSIS_BY_DOCETACEL_4NM_DN",
    "HERNANDEZ_ABERRANT_MITOSIS_BY_DOCETACEL_4NM_UP",
    "REACTOME_APC_C_CDH1_MEDIATED_DEGRADATION_OF_CDC20_AND_OTHER_APC_C_CDH1_TARGETED_PROTEINS_IN_LATE_MITOSIS_EARLY_G1",
    "REACTOME_FBXL7_DOWN_REGULATES_AURKA_DURING_MITOTIC_ENTRY_AND_IN_EARLY_MITOSIS",
    "REICHERT_MITOSIS_LIN9_TARGETS"
)

FOXO1_gene_list = c(
    "FOXO1",
    "PIK3CA",
    "AKT1",
    "AKT2",
    "PDPK1",
    "PTEN",
    "MAPK8",
    "MAPK9",
    "MAPK14",
    "STK4",
    "PRKAA1",
    "EP300",
    "CREBBP",
    "SIRT1",
    "SKP2",
    "MDM2",
    "YWHAZ",
    "XPO1",
    "TNPO1",
    "METTL3",
    "FTO",
    "USP7",
    "CTNNB1",
    "KAT2B",
    "SMAD3",
    "SIN3A",
    "PCK1",
    "PDK4",
    "FBP1",
    "PNPLA2",
    "PPARG",
    "SOD2",
    "CAT",
    "GADD45A",
    "MAP1LC3B",
    "GABARAPL1",
    "FBXO32",
    "TRIM63",
    "CDKN1B",
    "CCNG2",
    "RBL2",
    "BCL2L11",
    "BBC3",
    "TNFSF10",
    "IL1B",
    "TLR4",
    "INSR",
    "TCF7",
    "IL7R",
    "CCR7",
    "SELL",
    "POMC",
    "MYC"
)

old_mitosis_list <- c(
    "AKAP8",
    "ANAPC10",
    "ANAPC11",
    "ANAPC4",
    "ANAPC5",
    "ANLN",
    "ATM",
    "AURKA",
    "BIRC5",
    "BRSK1",
    "BUB1",
    "BUB1B",
    "CCNA2",
    "CD28",
    "CDC16",
    "CDC23",
    "CDC25B",
    "CDC25C",
    "CDC27",
    "CDCA5",
    "CDK13",
    "CDKN2B",
    "CENPE",
    "CETN1",
    "CHFR",
    "CHMP1A",
    "CIT",
    "CLIP1",
    "DCTN2",
    "DCTN3",
    "DDX11",
    "EGF",
    "EPGN",
    "EREG",
    "ESPL1",
    "GML",
    "KIF11",
    "KIF15",
    "KIF22",
    "KIF25",
    "KIF2C",
    "KNTC1",
    "MAD2L1",
    "MAD2L2",
    "NBN",
    "NCAPH",
    "NDC80",
    "NEK2",
    "NEK6",
    "NOLC1",
    "NPM2",
    "NUMA1",
    "NUSAP1",
    "PAM",
    "PBRM1",
    "PCBP4",
    "PDS5B",
    "PIN1",
    "PKMYT1",
    "PLK1",
    "PML",
    "PPP5C",
    "PRMT5",
    "RAD17",
    "RAN",
    "RCC1",
    "RINT1",
    "SMC1A",
    "SMC3",
    "SMC4",
    "SUGT1",
    "TARDBP",
    "TGFA",
    "TGFB1",
    "TPX2",
    "TRIAP1",
    "TTK",
    "TTN",
    "UBE2C",
    "ZNRD2",
    "ZW10",
    "ZWINT"
)

sanitize_label_for_filename <- function(label) {
    gsub("[/\\:*?\"<>|]", "_", label)
}

sanitize_label_for_AggregateExpression <- function(label) {
    gsub("[/\\:*?\"<>|_]", "-", label)
}

run_DE <- function(seurat_object, cell_type, outdir, min_cells){
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

    bulk <- AggregateExpression(
        seurat_object,
        return.seurat = TRUE,
        slot = "counts",
        assays = "RNA",
        group.by = c("singleR.labels_fine", "sample_id", "diagnosis")
    )

    Idents(bulk) <- "diagnosis"
    print(unique(Idents(bulk)))
    print(length(Idents(bulk)))
    print(head(Idents(bulk)))
    print(table(bulk$diagnosis))

    file_name <- file.path(outdir, "aggregate_cell_counts.csv")
    write.csv(table(bulk$diagnosis), file_name)


    # add a pseudocount of 1 to every gene as per https://www.biostars.org/p/440379/
    counts_matrix <- GetAssayData(bulk, layer = "counts")
    bulk <- SetAssayData(bulk, layer = "counts", new.data = counts_matrix + 1)

    n_cells <- sum(bulk$diagnosis == cell_type)
    if (n_cells < min_cells){
        print(n_cells)
        print(paste0(cell_type, " had less than ", min_cells," cells"))
        print(unique(bulk$diagnosis))
    }else{
        print(cell_type)
        markers <- FindMarkers(
            bulk,
            ident.1 = cell_type,
            ident.2 = NULL,
            slot = "counts",
            test.use = "DESeq2"
        )
        file_path <- file.path(outdir, paste0(sanitize_label_for_filename(cell_type), "_vs_control_diffexpress.csv"))
        write.csv(markers, file_path)
    }
    return(markers)
}

run_fgsea <- function(markers, pathways, file_prefix, threads = 8, top_up = 20, top_down = 20){
    # following https://biostatsquid.com/fgsea-tutorial-gsea/
    markers <- na.omit(markers)

    print(head(markers))
    rankings <- sign(markers$avg_log2FC) * (-log10(markers$p_val_adj))
    names(rankings) <- rownames(markers)
    rankings <- sort(rankings, decreasing = TRUE)

    print(head(rankings))

    max_ranking <- max(rankings[is.finite(rankings)])
    print(paste0("max rank is: ", max_ranking))
    min_ranking <- min(rankings[is.finite(rankings)])
    print(paste0("min rank is: ", min_ranking))
    if (max_ranking == 0){
        message("No up regulated genes")
        quit(status = 1)
    }
    rankings <- replace(rankings, rankings > max_ranking, max_ranking * 10)
    rankings <- replace(rankings, rankings < min_ranking, min_ranking * 10)
    rankings <- sort(rankings, decreasing = TRUE) # sort genes by rankk

    print(head(rankings))

    GSEAres <- fgsea(pathways = pathways, # List of gene sets to check
                     stats = rankings,
                     scoreType = 'std', # in this case we have both pos and neg rankings. if only pos or neg, set to 'pos', 'neg'
                     minSize = 10,
                     maxSize = 500,
                     nproc = threads
    )

    print(GSEAres)

    print(sum(GSEAres[, padj < 0.01]))
    print(sum(GSEAres[, pval < 0.01]))
    print(min(GSEAres$padj))
    print(min(GSEAres$pval))

    number_of_top_pathways_up <- min(top_up, length(GSEAres$pathway))
    number_of_top_pathways_down <- min(top_down, length(GSEAres$pathway))

    topPathwaysUp <- GSEAres[ES > 0][head(order(padj), n = number_of_top_pathways_up), pathway]
    topPathwaysDown <- GSEAres[ES < 0][head(order(padj), n = number_of_top_pathways_down), pathway]
    topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
    print(topPathways)

    plot <- plotGseaTable(pathways[topPathways], stats = rankings, fgseaRes = GSEAres, gseaParam = 0.5)
    ggsave(paste0(file_prefix, ".pdf"), plot, width = 16, height = 8, dpi = 300)


    GSEAres <- GSEAres[order(padj)]
    GSEAres$leadingEdge <- sapply(GSEAres$leadingEdge, paste, collapse = ";")

    write.csv(GSEAres, paste0(file_prefix, ".csv"))
}

# Create an ArgumentParser object
parser <- ArgumentParser(description = 'Integrate multiple scRNA-seq datasets into one Seurat object')
parser$add_argument("joined_integrated_seurat_object", help="File paths to joined_integrated_seurat_object.rds")
parser$add_argument("--cell_type", help="cell type to use")
parser$add_argument("--min_cells", type="integer", help="Min cells per celltype for diff express", default = 8)

# Parse the command-line arguments
args <- parser$parse_args()

print(args$cell_type)

joined_integrated_data <- readRDS(args$joined_integrated_seurat_object)
msig <- msigdbr(species = "Homo sapiens")

# Filter to your sets and convert to named list
pathways <- msig %>%
    filter(gs_name %in% target_sets) %>%
    split(x = .$gene_symbol, f = .$gs_name)

# Check which ones were found
cat("Found", length(pathways), "of", length(target_sets), "gene sets\n")
setdiff(target_sets, names(pathways))  # Shows any missing sets

pathways[["CUSTOM_FOXO1_LIST"]] <- FOXO1_gene_list
pathways[["CUSTOM_OLD_MITOSIS_LIST"]] <- old_mitosis_list

cat("Found", length(pathways), "of", length(target_sets), "gene sets\n")
setdiff(target_sets, names(pathways))  # Shows any missing sets

just_fibroblast_cells <- subset(
    joined_integrated_data,
    singleR.labels_fine %in% fibroblast_celltypes
)


# run on all gene sets
all_pathways <- msig %>%
  split(x = .$gene_symbol, f = .$gs_name)

just_fibroblast_cells$singleR.labels_fine <- sanitize_label_for_AggregateExpression(just_fibroblast_cells$singleR.labels_fine)

just_cell_type <- subset(just_fibroblast_cells, singleR.labels_fine == args$cell_type)
fib_markers <- run_DE(just_cell_type, "SSC", paste0(args$cell_type, "_DE"), args$min_cells)
run_fgsea(fib_markers, pathways, paste0("targeted_gsea_results_just_", args$cell_type))
run_fgsea(fib_markers, all_pathways, paste0("all_gsea_results_just_", args$cell_type))
