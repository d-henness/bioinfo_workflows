library(argparse)
library(Seurat)
library(ggplot2)
library(SingleR)
library(celldex)
library(tidyverse)
library(pheatmap)


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

    Idents(bulk) <- "singleR.labels_fine"
    print(unique(Idents(bulk)))
    print(length(Idents(bulk)))
    print(table(bulk$singleR.labels_fine))

    file_name <- file.path(outdir, "aggregate_cell_counts.csv")
    write.csv(table(bulk$singleR.labels_fine), file_name)


    # add a pseudocount of 1 to every gene as per https://www.biostars.org/p/440379/
    counts_matrix <- GetAssayData(bulk, layer = "counts")
    bulk <- SetAssayData(bulk, layer = "counts", new.data = counts_matrix + 1)

    n_cells <- sum(bulk$singleR.labels_fine == cell_type)
    if (n_cells < min_cells){
        print(paste0(cell_type, " had less than ", min_cells," cells"))
        print(unique(bulk$singleR.labels_fine))
    }else{
        print(cell_type)
        markers <- FindMarkers(
            bulk,
            ident.1 = cell_type,
            ident.2 = NULL,
            slot = "counts",
            test.use = "DESeq2"
        )
        file_path <- file.path(outdir, paste0(sanitize_label_for_filename(cell_type), "_vs_all_diffexpress.csv"))
        write.csv(markers, file_path)
    }
}

# Create an ArgumentParser object
parser <- ArgumentParser(description = 'Integrate multiple scRNA-seq datasets into one Seurat object')
parser$add_argument("joined_integrated_seurat_object", help="File paths to joined_integrated_seurat_object.rds")
parser$add_argument("--cell_type", help="cell type to use")
parser$add_argument("--min_cells", type="integer", help="Min cells per celltype for diff express", default = 8)

# Parse the command-line arguments
args <- parser$parse_args()

joined_integrated_data <- readRDS(args$joined_integrated_seurat_object)

print(unique(sanitize_label_for_AggregateExpression(joined_integrated_data$singleR.labels_fine)))

just_fibroblast_cells <- subset(
    joined_integrated_data,
    singleR.labels_fine %in% fibroblast_celltypes
)

ssc_cells <- subset(
    joined_integrated_data,
    diagnosis == "SSC"
)

control_cells <- subset(
    joined_integrated_data,
    diagnosis == "control"
)

just_fibroblast_ssc_cells <- subset(
    just_fibroblast_cells,
    diagnosis == "SSC"
)

just_fibroblast_control_cells <- subset(
    just_fibroblast_cells,
    diagnosis == "control"
)

ssc_cells$singleR.labels_fine <- sanitize_label_for_AggregateExpression(ssc_cells$singleR.labels_fine)
control_cells$singleR.labels_fine <- sanitize_label_for_AggregateExpression(control_cells$singleR.labels_fine)
just_fibroblast_ssc_cells$singleR.labels_fine <- sanitize_label_for_AggregateExpression(just_fibroblast_ssc_cells$singleR.labels_fine)
just_fibroblast_control_cells$singleR.labels_fine <- sanitize_label_for_AggregateExpression(just_fibroblast_control_cells$singleR.labels_fine)


run_DE(ssc_cells, args$cell_type, "ssc_cells_single_type_vs_all_others", args$min_cells)
run_DE(control_cells, args$cell_type, "control_cells_single_type_vs_all_others", args$min_cells)
run_DE(just_fibroblast_ssc_cells, args$cell_type, "just_fibroblast_ssc_cells_single_type_vs_all_fibroblasts", args$min_cells)
run_DE(just_fibroblast_control_cells, args$cell_type, "just_fibroblast_control_cells_single_type_vs_all_fibroblasts", args$min_cells)
