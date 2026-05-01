library(argparse)
library(Seurat)
library(ggplot2)
library(SingleR)
library(celldex)
library(tidyverse)
library(pheatmap)
library(EnhancedVolcano)
library(grid)
library(gtable)
library(RColorBrewer)

options(warn = 1)

mapping <- read.csv("cell_type_mapping.csv")

gene_sets <- list(
    chromatin_modifiers = c(
        "HEXIM1",
        "DICER1",
        "ARID2",
        "ARID1A",
        "CHD1",
        "CUL3",
        "KMT2C",
        "SETD2",
        "EZH2",
        "BRD4",
        "ARID5B",
        "EP300"
  )
    
)

sanitize_label_for_filename <- function(label) {
    gsub("[/\\:*?\"<>|]", "_", label)
}

sanitize_label_for_AggregateExpression <- function(label) {
    gsub("[/\\:*?\"<>|_]", "-", label)
}

make_gene_set_plots <- function(markers,
                                outdir,
                                min_val = -2,
                                max_val = 2,
                                col_annotation = NULL,   # data.frame, rownames = colnames(markers)
                                ann_colors     = NULL,   # named list of color vectors
                                cutree_cols    = NA,     # e.g. 3 to color dendrogram branches
                                highlight_cols = NULL) { # named list: list(red = c("col1","col2"), ...)
 
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
 
  color_palette = colorRampPalette(c("Darkblue", "white","red"))(100)
  breaks <- seq(min_val, max_val, length.out = 101)
 
  add_boxes <- function(ph, mat, highlight_cols) {
    if (is.null(highlight_cols)) return(ph$gtable)
    col_order <- if (!is.null(ph$tree_col)) colnames(mat)[ph$tree_col$order] else colnames(mat)
    n_cols <- length(col_order)
    boxes <- grobTree()
    for (color in names(highlight_cols)) {
      positions <- which(col_order %in% highlight_cols[[color]])
      if (length(positions) == 0) next
      x_start <- (min(positions) - 1) / n_cols
      x_end   <- max(positions)       / n_cols
      boxes <- addGrob(boxes, rectGrob(
        x = unit(x_start, "npc"), y = unit(0, "npc"),
        width  = unit(x_end - x_start, "npc"),
        height = unit(1, "npc"),
        just = c("left", "bottom"),
        gp = gpar(col = color, fill = NA, lwd = 2)
      ))
    }
    matrix_pos <- which(ph$gtable$layout$name == "matrix")
    gt <- ph$gtable
    gt <- gtable_add_grob(gt, boxes,
                          t = gt$layout$t[matrix_pos],
                          l = gt$layout$l[matrix_pos],
                          b = gt$layout$b[matrix_pos],
                          r = gt$layout$r[matrix_pos],
                          name = "highlight_boxes")
    gt
  }
 
  draw_one <- function(mat, title, filenm) {
    ph <- pheatmap(
      mat,
      cluster_rows        = TRUE,
      cluster_cols        = TRUE,
      scale               = "none",
      main                = title,
      fontsize_row        = 8,
      fontsize_col        = 8,
      color               = color_palette,
      breaks              = breaks,
      annotation_col      = col_annotation,
      annotation_colors   = ann_colors,
      cutree_cols         = cutree_cols,
      angle_col           = 90,
      border_color        = "grey80",
      legend_breaks       = c(min_val, 0, max_val),
      silent              = TRUE
    )
    final <- add_boxes(ph, mat, highlight_cols)
    height <- max(5, nrow(mat) / 5)
    width  <- max(6, ncol(mat) / 3)
    ggsave(filenm, final, width = width, height = height, limitsize = FALSE)
  }
 
  for (gene_set in names(gene_sets)) {
    genes <- intersect(gene_sets[[gene_set]], rownames(markers))
    if (length(genes) == 0) next
    mat <- markers[genes, , drop = FALSE]
    filenm <- file.path(outdir, paste0("heatmap_diff_express_", gene_set, ".pdf"))
    draw_one(mat, paste0("avg_log2FC (", gene_set, ")"), filenm)
  }
 
  all_genes <- intersect(unlist(gene_sets), rownames(markers))
  if (length(all_genes) > 0) {
    draw_one(markers[all_genes, , drop = FALSE],
             "avg_log2FC (all gene sets)",
             file.path(outdir, "heatmap_diff_express_all_gene_sets.pdf"))
  }
}

run_DE <- function(seurat_object, cell_type, identity_of_interest, outdir, file_suff, min_cells){
    dir.create(file.path(outdir, "filtered"), showWarnings = FALSE, recursive = TRUE)

    cells_per_sample <- table(seurat_object$sample_id, seurat_object[[identity_of_interest, drop = TRUE]])
    file_name <- file.path(outdir, paste0("cell_counts_", cell_type, "_", file_suff, ".csv"))
    write.csv(cells_per_sample, file_name)

    print("Before filter")
    print(table(seurat_object$sample_id, seurat_object[[identity_of_interest, drop = TRUE]]))
    #filter out bad samples
    good_samples <- rownames(cells_per_sample)[rowSums(cells_per_sample) >= 10]
    seurat_object <- subset(seurat_object, sample_id %in% good_samples)
    print("After filter")
    print(table(seurat_object$sample_id, seurat_object[[identity_of_interest, drop = TRUE]]))

    if (length(unique(seurat_object$sample_id)) < min_cells) {
        print(paste0("Skipping ", cell_type, ": only ", length(unique(seurat_object$sample_id)), " sample(s) after filtering — need > 2 ", min_cells, " for pseudobulk DE"))
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
    n_control_cells <- sum(bulk[[identity_of_interest, drop = TRUE]] == "normal-skin")
    print(n_cells)
    markers <- NULL
    if (n_cells < min_cells || n_control_cells < 3){
        print(paste0(cell_type, " had less than ", min_cells," cells, or normal-skin had fewer than 3 pseudobulk samples (", n_control_cells, ")"))
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
parser$add_argument("--min_cells", type="integer", help="Min cells per celltype for diff express", default = 3)
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
            print(paste0("running: ", new_cat_label, " ", diag))
            diag_and_new_cat = paste0(diag, "_", new_cat_label)
            all_markers[[diag_and_new_cat]] <- run_DE(just_subset_new, diag, "diagnosis", file.path(args$outdir, new_cat_label), "vs_normal_skin", args$min_cells)
        }
    }
}

all_markers <- Filter(Negate(is.null), all_markers)

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

build_col_ann <- function(col_names, all_diagnoses) {
    diags    <- sapply(col_names, function(col) Filter(function(d) startsWith(col, d), all_diagnoses)[[1]])
    new_cats <- mapply(function(col, diag) sub(paste0("^", diag, "_"), "", col), col_names, diags)
    data.frame(
        Condition          = diags,
        `Fibroblast Class` = new_cats,
        row.names          = col_names,
        check.names        = FALSE
    )
}

build_ann_colors <- function(col_ann) {
    classes <- unique(col_ann[["Fibroblast Class"]])
    n <- length(classes)
    hues <- hsv(seq(0, 1, length.out = n + 1)[seq_len(n)], s = 0.75, v = 0.9)
    names(hues) <- classes
    list(`Fibroblast Class` = hues)
}

for (new_cat_label in unique(just_subset$new_cat_labels)){
    new_cat_label_cols <- grep(paste0("_", new_cat_label, "$"), colnames(marker_matrix), value = TRUE)
    print(new_cat_label_cols)
    if (length(new_cat_label_cols) > 0) {
        new_cat_label_matrix <- marker_matrix[, new_cat_label_cols, drop = FALSE]
        col_ann <- build_col_ann(colnames(new_cat_label_matrix), unique(just_subset$diagnosis))
        make_gene_set_plots(new_cat_label_matrix, file.path(args$outdir, new_cat_label),
                            col_annotation = col_ann, ann_colors = build_ann_colors(col_ann))
    }
}

col_ann <- build_col_ann(colnames(marker_matrix), unique(just_subset$diagnosis))
make_gene_set_plots(marker_matrix, file.path(args$outdir, "all_combined"),
                    col_annotation = col_ann, ann_colors = build_ann_colors(col_ann))
