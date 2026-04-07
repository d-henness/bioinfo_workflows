library(argparse)
library(Seurat)
library(ggplot2)
library(SingleR)
library(celldex)
library(tidyverse)
library(pheatmap)
library(msigdbr)
library(fgsea)
library(EnhancedVolcano)
library(lubridate)

mapping <- read.csv("cell_type_mapping.csv")

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

run_DE <- function(seurat_object, cell_type1, identity_of_interest, outdir, min_cells){
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

    bulk <- AggregateExpression(
        seurat_object,
        return.seurat = TRUE,
        slot = "counts",
        assays = "RNA",
        group.by = c(identity_of_interest, "sample_id")
    )

    Idents(bulk) <- identity_of_interest
    print(unique(Idents(bulk)))
    print(length(Idents(bulk)))
    print(head(Idents(bulk)))
    print(table(bulk[[identity_of_interest]]))

    file_name <- file.path(outdir, "aggregate_cell_counts.csv")
    write.csv(table(bulk[[identity_of_interest]]), file_name)

    # add a pseudocount of 1 to every gene as per https://www.biostars.org/p/440379/
    counts_matrix <- GetAssayData(bulk, layer = "counts")
    bulk <- SetAssayData(bulk, layer = "counts", new.data = counts_matrix + 1)

    n_cells <- sum(bulk[[identity_of_interest]] == cell_type1)
    markers <- NULL
    if (n_cells < min_cells){
        print(n_cells)
        print(paste0(cell_type1, " had less than ", min_cells," cells"))
        print(unique(bulk[[identity_of_interest]]))
    }else{
        print(cell_type1)
        markers <- FindMarkers(
            bulk,
            ident.1 = cell_type1,
            slot = "counts",
            test.use = "DESeq2"
        )
        file_path <- file.path(outdir, paste0(cell_type1, "_vs_all_diffexpress.csv"))
        write.csv(markers, file_path)

        volcano_plot <- EnhancedVolcano(
            markers,
            lab = rownames(markers),
            x = "avg_log2FC",
            y = "p_val_adj",
            title = cell_type1,
            pCutoff = 0.05,
            FCcutoff = 1
        )
        plot_path <- file.path(outdir, paste0(cell_type1, "_vs_all_volcano.pdf"))
        ggsave(plot_path, volcano_plot, width = 10, height = 8)
    }
    return(markers)
}

run_fgsea <- function(markers, pathways, outdir, file_prefix, threads = 8, top_up = 20, top_down = 20){
    dir.create(file.path(outdir, "plots"), showWarnings = FALSE, recursive = TRUE)

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
    
    if (min_ranking == 0 && max_ranking == 0) {
        print(paste0("No up or down ranked genes, skipping ", outdir))
        return(NULL)
    } else if (max_ranking <= 0) {
        score_type <- "neg"
    } else if (min_ranking >= 0) {
        score_type <- "pos"
    } else {
        score_type <- "std"
    }

    rankings <- replace(rankings, rankings > max_ranking, max_ranking * 10)
    rankings <- replace(rankings, rankings < min_ranking, min_ranking * 10)
    rankings <- sort(rankings, decreasing = TRUE) # sort genes by rankk

    print(head(rankings))

    GSEAres <- fgsea(pathways = pathways, # List of gene sets to check
                     stats = rankings,
                     scoreType = score_type,
                     minSize = 10,
                     maxSize = 500,
                     nproc = threads
    )

    print(GSEAres)

    print(sum(GSEAres[, padj < 0.01], na.rm = TRUE))
    print(sum(GSEAres[, pval < 0.01], na.rm = TRUE))
    print(min(GSEAres$padj, na.rm = TRUE))
    print(min(GSEAres$pval, na.rm = TRUE))

    number_of_top_pathways_up <- min(top_up, length(GSEAres$pathway))
    number_of_top_pathways_down <- min(top_down, length(GSEAres$pathway))

    topPathwaysUp <- GSEAres[ES > 0][head(order(padj), n = number_of_top_pathways_up), pathway]
    topPathwaysDown <- GSEAres[ES < 0][head(order(padj), n = number_of_top_pathways_down), pathway]
    topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
    print(topPathways)

    plot <- plotGseaTable(pathways[topPathways], stats = rankings, fgseaRes = GSEAres, gseaParam = 0.5)
    file_path <- file.path(outdir, paste0(file_prefix,".pdf"))
    print(file_path)
    ggsave(file_path, plot, width = 16, height = 8, dpi = 300)


    GSEAres <- GSEAres[order(padj)]
    GSEAres$leadingEdge <- sapply(GSEAres$leadingEdge, paste, collapse = ";")

    file_path <- file.path(outdir, paste0(file_prefix,".csv"))
    print(file_path)
    write.csv(GSEAres, file_path)

    sig_pathways <- GSEAres[padj < 0.05, pathway]

    for (pathway_name in sig_pathways) {
      plot <- plotEnrichment(pathways[[pathway_name]], rankings) +
        ggtitle(pathway_name)
      
      file_path <- file.path(outdir, "plots", paste0(pathway_name, "_enrichment.pdf"))
      ggsave(file_path, plot, width = 8, height = 6, dpi = 300)
    }

}

# Create an ArgumentParser object
parser <- ArgumentParser(description = 'Integrate multiple scRNA-seq datasets into one Seurat object')
parser$add_argument("--min_cells", type="integer", help="Min cells per celltype for diff express", default = 5)
parser$add_argument("--outdir", help="directory to write to", required = TRUE)
parser$add_argument("--file_pref", help="file prefix", required = TRUE)
parser$add_argument("--diagnosis", help="diagnosis", required = TRUE)
# Parse the command-line arguments
args <- parser$parse_args()

print(args)

joined_integrated_data <- readRDS("rds_files/joined_integrated_seurat_object_with_group_data.rds")
wanted_groups = c("fibroblasts", "pericyte")

just_subset <- subset(joined_integrated_data, cell_group %in% wanted_groups)
just_subset$m_labels <- mapping$m_labels[match(just_subset$singleR.labels_fine, mapping$singleR.labels_fine)]
just_subset$new_cat_labels <- mapping$new_cat_labels[match(just_subset$singleR.labels_fine, mapping$singleR.labels_fine)]

msig <- msigdbr(species = "Homo sapiens")

all_pathways <- msig %>%
    split(x = .$gene_symbol, f = .$gs_name)

pathways <- msig %>%
    filter(gs_name %in% target_sets) %>%
    split(x = .$gene_symbol, f = .$gs_name)

# Check which ones were found
cat("Found", length(pathways), "of", length(target_sets), "gene sets\n")
setdiff(target_sets, names(pathways))  # Shows any missing sets

pathways[["CUSTOM_FOXO1_LIST"]] <- FOXO1_gene_list
pathways[["CUSTOM_OLD_MITOSIS_LIST"]] <- old_mitosis_list

cat("Found", length(pathways), "of", length(target_sets), "gene sets\n")

just_subset$singleR.labels_fine <- sanitize_label_for_AggregateExpression(just_subset$singleR.labels_fine)

just_subset_diag <- subset(just_subset, diagnosis == args$diagnosis)
print(paste0("Will runs on these types: ", unique(just_subset_diag$singleR.labels_fine)))

for (cell_type in unique(just_subset_diag$singleR.labels_fine)){
    print(paste0("Will run on: ", cell_type))
    outdir <- file.path(args$outdir, args$diagnosis, cell_type, "DE")
    markers <- run_DE(just_subset_diag, cell_type, "singleR.labels_fine", outdir, args$min_cells)
    if (!is.null(markers)){
        outdir <- file.path(args$outdir, args$diagnosis, cell_type, "fgsea_target")
        run_fgsea(markers, pathways, outdir, paste0(args$file_pref, "_", cell_type, "_", args$diagnosis))
        outdir <- file.path(args$outdir, args$diagnosis, cell_type, "fgsea_all")
        run_fgsea(markers, all_pathways, outdir, paste0(args$file_pref, "_", cell_type, "_", args$diagnosis))
    }
}
print("Ran to the end")
