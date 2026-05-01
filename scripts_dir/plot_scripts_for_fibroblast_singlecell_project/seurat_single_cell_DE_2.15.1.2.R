library(argparse)
library(Seurat)
library(ggplot2)
library(SingleR)
library(celldex)
library(tidyverse)
library(pheatmap)
library(EnhancedVolcano)

options(warn = 1)

mapping <- read.csv("cell_type_mapping.csv")

gene_sets <- list(
    mitochondrial_genes = c(
        "MT-ATP6",
        "MT-CO2",
        "MT-CO1",
        "MT-ND2", 
        "MT-ND4",
        "MT-ND5",
        "MT-CYB",
        "MT-ATP8",
        "MT-CO3",
        "MT-ND3",
        "MT-ND1",
        "MT-ND4L",
        "MT-ND6"
    ),
    chromosome_X_genes = c(
        "PDHA1",
        "AIFM1",
        "APOOL",
        "IDH3G",
        "HCCS",
        "ABCB7",
        "SLC25A5",
        "PDK3",
        "NDUFB11",
        "MPC1L",
        "TIMM17B",
        "GK",
        "HSD17B10",
        "TIMM8A",
        "COX7B",
        "TRMT2B",
        "CA5B",
        "ALAS2",
        "FUNDC2",
        "APOO",
        "SLC25A14",
        "NDUFA1"
    ),
    somatic_chromosomes = c(
        "PRODH",
        "GLDC",
        "OGDHL",
        "ACADL",
        "SLC25A35",
        "LIPT2",
        "SLC25A15",
        "MRM1",
        "SLC25A10",
        "MARS2",
        "TOMM40L",
        "ME3",
        "PPIF",
        "ALDH4A1",
        "SLC25A30",
        "CMC2",
        "COQ10A",
        "C6orf136",
        "SDHAF4",
        "ABCB8",
        "SPRYD4",
        "IBA57",
        "CLPB",
        "SLC25A42",
        "RTN4IP1",
        "ACAD10",
        "SLC25A20",
        "ADCK5",
        "GATD3A",
        "BDH1",
        "CPT2",
        "GUF1",
        "MFN1",
        "FECH",
        "SLC25A19",
        "SUPV3L1",
        "AFG1L",
        "GATC",
        "GLS",
        "DHTKD1",
        "GATB",
        "FDXR",
        "SLC25A25",
        "ALAS1",
        "MRPL44",
        "CS",
        "LARS2",
        "ERAL1",
        "FXN",
        "HSD17B8",
        "MRPL30",
        "PCCB",
        "YARS2",
        "DLAT",
        "ATPAF2",
        "MTO1",
        "HSPD1",
        "LACTB",
        "COQ3",
        "L2HGDH",
        "NFS1",
        "COQ6",
        "TIMM21",
        "MRPL53",
        "PITRM1",
        "MTIF2",
        "ABHD11",
        "ADHFE1",
        "NDUFAF5",
        "MCAT",
        "ACSF3",
        "LETM1",
        "POLRMT",
        "TACO1",
        "PDK1",
        "PPTC7",
        "LIAS",
        "MRPL46",
        "MTERF2",
        "GCDH",
        "DNAJA3",
        "COX15",
        "ISCA1",
        "MRPS17",
        "GFM2",
        "NDUFAF4",
        "ACAD8",
        "ECHDC3",
        "IDH2",
        "LONP1",
        "MRPL17",
        "CRAT",
        "DBT"
    ), transcription_factors = c( "ARID5B", "EP300", "KDM5C", "KMT2C", "NCOR1", "NSD1", "SETD2", "ARID1A", "ARID2", "CHD1", "CHD8", "H3F3C", "CTCF",
        "BRCA1",
        "POLQ",
        "SMC1A",
        "CUL3",
        "USP9X",
        "MGA",
        "DHX9",
        "DICER1",
        "XPO1",
        "BRD4",
        "SUPT6H",
        "HEXIM1",
        "LARP7",
        "EZH2",
        "DOT1L",
        "SUPT16H",
        "CDK9",
        "HEXIM1"
    )
)

sanitize_label_for_filename <- function(label) {
    gsub("[/\\:*?\"<>|]", "_", label)
}

sanitize_label_for_AggregateExpression <- function(label) {
    gsub("[/\\:*?\"<>|_]", "-", label)
}

make_gene_set_plots <- function(markers, outdir, min_val = -2, max_val = 2){
    dir.create(file.path(outdir), showWarnings = FALSE, recursive = TRUE)
    color_palette = colorRampPalette(c("Darkblue", "white","red"))(100)
    breaks <- seq(min_val, max_val, length.out = 101)
    for (gene_set in names(gene_sets)){
        just_gene_set <- markers[intersect(gene_sets[[gene_set]], rownames(markers)), , drop = FALSE]

        plot <- pheatmap(
            just_gene_set,
            cluster_rows = TRUE,
            cluster_cols = TRUE,
            scale        = "none",
            main         = paste0("avg_log2FC (", gene_set, ")"),
            fontsize_row = 8,
            fontsize_col = 8,
            color        = color_palette,
            breaks       = breaks
        )

        filenm <- file.path(outdir, paste0("heatmap_diff_express_", gene_set, ".pdf"))
        print(filenm)
        height <- max(5, length(gene_sets[[gene_set]]) / 5)
        width <- length(colnames(markers)) / 5
        print(paste(filenm, height, width, length(gene_sets[[gene_set]]), length(colnames(markers))))
        ggsave(filenm, plot, width = width, height = height, limitsize = FALSE)
    }
}

run_DE <- function(seurat_object, cell_type, identity_of_interest, outdir, file_suff, min_cells){
    dir.create(file.path(outdir, "filtered"), showWarnings = FALSE, recursive = TRUE)

    cells_per_sample <- table(seurat_object$sample_id, seurat_object[[identity_of_interest, drop = TRUE]])
    file_name <- file.path(outdir, paste0("cell_counts_", cell_type, "_", file_suff, ".csv"))
    write.csv(cells_per_sample, file_name)

    #filter out bad samples
    good_samples <- rownames(cells_per_sample)[rowSums(cells_per_sample) >= 10]
    seurat_object <- subset(seurat_object, sample_id %in% good_samples)

    if (length(unique(seurat_object$sample_id)) < 2) {
        print(paste0("Skipping ", cell_type, ": only ", length(unique(seurat_object$sample_id)), " sample(s) after filtering — need >= 2 for pseudobulk DE"))
        return(NULL)
    }

    bulk <- AggregateExpression(
        seurat_object,
        return.seurat = TRUE,
        slot = "counts",
        assays = "RNA",
        group.by = c(identity_of_interest, "sample_id")
    )

    Idents(bulk) <- identity_of_interest
    print(head(bulk))

    file_name <- file.path(outdir, paste0("aggregate_cell_counts_", file_suff, ".csv"))
    write.csv(table(bulk[[identity_of_interest]]), file_name)

#    # add a pseudocount of 1 to every gene as per https://www.biostars.org/p/440379/
#    counts_matrix <- GetAssayData(bulk, layer = "counts")
#    bulk <- SetAssayData(bulk, layer = "counts", new.data = counts_matrix + 1)

    n_cells <- sum(bulk[[identity_of_interest, drop = TRUE]] == cell_type)
    print(n_cells)
    if (n_cells < min_cells){
        print(paste0(cell_type, " had less than ", min_cells," cells"))
        print(table(bulk[[identity_of_interest]]))
    }else{
        print(cell_type)
        markers <- FindMarkers(
            bulk,
            ident.1 = cell_type,
            ident.2 = sanitize_label_for_AggregateExpression("normal_skin"),
            slot = "counts",
            test.use = "DESeq2"
        )
        file_path <- file.path(outdir, paste0(sanitize_label_for_filename(cell_type), "_", file_suff, ".csv"))
        write.csv(markers, file_path)

        volcano_plot <- EnhancedVolcano(
            markers,
            lab = rownames(markers),
            x = "avg_log2FC",
            y = "p_val_adj",
            title = cell_type,
            pCutoff = 0.05,
            FCcutoff = 1
        )
        plot_path <- file.path(outdir, paste0(sanitize_label_for_filename(cell_type), "_", file_suff, "_volcano.pdf"))
        ggsave(plot_path, volcano_plot, width = 10, height = 8)

        markers_filtered <- markers[!is.na(markers$p_val_adj) & markers$p_val_adj < 0.05 & rownames(markers) %in% unlist(gene_sets), ]
        file_path <- file.path(outdir, "filtered", paste0(sanitize_label_for_filename(cell_type), "_filtered_", file_suff, ".csv"))
        write.csv(markers_filtered, file_path)
    }
    return(markers)
}

# Create an ArgumentParser object
parser <- ArgumentParser(description = 'Integrate multiple scRNA-seq datasets into one Seurat object')
parser$add_argument("--outdir", help="directory to write to", required = TRUE)
parser$add_argument("--min_cells", type="integer", help="Min cells per celltype for diff express", default = 5)
args <- parser$parse_args()

print(args)

joined_integrated_data <- readRDS("rds_files/joined_integrated_seurat_object_with_group_data.rds")
wanted_groups = c("fibroblasts", "pericyte")

just_subset <- subset(joined_integrated_data, cell_group %in% wanted_groups)
just_subset$m_labels <- mapping$m_labels[match(just_subset$singleR.labels_fine, mapping$singleR.labels_fine)]
just_subset$new_cat_labels <- mapping$new_cat_labels[match(just_subset$singleR.labels_fine, mapping$singleR.labels_fine)]

all_markers <- list()
for (new_cat_label in unique(just_subset$new_cat_labels)){
    just_subset_new <- subset(just_subset, new_cat_labels == new_cat_label)
    just_subset_new$new_cat_labels <- sanitize_label_for_AggregateExpression(just_subset_new$new_cat_labels)

    for (diag in unique(just_subset_new$diagnosis)){
        if (diag != "normal_skin") {
            diag_and_new_cat = paste0(diag, "_", new_cat_label)
            all_markers[[diag_and_new_cat]] <- run_DE(just_subset_new, diag, "diagnosis", file.path(args$outdir, new_cat_label), "vs_normal_skin", args$min_cells)
        }
    }
}

avail_genes <- Reduce(intersect, lapply(all_markers, rownames))
marker_matrix <- data.frame()
for (diag_and_new_cat in names(all_markers)) {
    diag    <- Filter(function(d) startsWith(diag_and_new_cat, d), unique(just_subset$diagnosis))[[1]]
    new_cat <- sub(paste0("^", diag, "_"), "", diag_and_new_cat)
    new_df <- data.frame(
        genes     = avail_genes,
        diag      = diag,
        new_cat   = new_cat,
        log2FC    = all_markers[[diag_and_new_cat]][avail_genes, "avg_log2FC"],
        p_val     = all_markers[[diag_and_new_cat]][avail_genes, "p_val"],
        p_val_adj = all_markers[[diag_and_new_cat]][avail_genes, "p_val_adj"]
    )
    marker_matrix <- rbind(marker_matrix, new_df)
}

filenm <- file.path(args$outdir, "all_fold_change_pval.csv")
write.csv(marker_matrix, filenm)

for (gene_set in names(gene_sets)){
    avail_genes    <- Reduce(intersect, lapply(all_markers, rownames))
    gene_set_genes <- intersect(gene_sets[[gene_set]], avail_genes)
    marker_matrix  <- data.frame()
    for (diag_and_new_cat in names(all_markers)) {
        diag    <- Filter(function(d) startsWith(diag_and_new_cat, d), unique(just_subset$diagnosis))[[1]]
        new_cat <- sub(paste0("^", diag, "_"), "", diag_and_new_cat)
        new_df  <- data.frame(
            genes     = gene_set_genes,
            diag      = diag,
            new_cat   = new_cat,
            log2FC    = all_markers[[diag_and_new_cat]][gene_set_genes, "avg_log2FC"],
            p_val     = all_markers[[diag_and_new_cat]][gene_set_genes, "p_val"],
            p_val_adj = all_markers[[diag_and_new_cat]][gene_set_genes, "p_val_adj"]
        )
        marker_matrix <- rbind(marker_matrix, new_df)
    }
    filenm <- file.path(args$outdir, paste0(gene_set, "_fold_change_pval.csv"))
    write.csv(marker_matrix, filenm)
}

avail_genes <- Reduce(intersect, lapply(all_markers, rownames))
marker_matrix <- data.frame(row.names = avail_genes)
for (diag_and_new_cat in names(all_markers)) {
    marker_matrix[[diag_and_new_cat]]  <- all_markers[[diag_and_new_cat]][avail_genes, "avg_log2FC"]   
}

for (new_cat_label in unique(just_subset$new_cat_labels)){
    new_cat_label_cols <- grep(paste0("^", new_cat_label, "_"), colnames(marker_matrix), value = TRUE)
    print(new_cat_label_cols)
    if (length(new_cat_label_cols) > 0) {
        new_cat_label_matrix <- marker_matrix[, new_cat_label_cols, drop = FALSE]
        make_gene_set_plots(new_cat_label_matrix, file.path(args$outdir, new_cat_label))
    }
}

make_gene_set_plots(marker_matrix, file.path(args$outdir, "all_combined"))
