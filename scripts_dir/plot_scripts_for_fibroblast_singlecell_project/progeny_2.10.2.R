library(progeny)
library(dplyr)
library(Seurat)
library(ggplot2)
library(tidyr)
library(readr)
library(pheatmap)
library(tibble)
library(argparse)

mapping <- read.csv("cell_type_mapping.csv")


sanitize_label_for_filename <- function(label) {
    gsub("[/\\:*?\"<>|]", "_", label)
}

sanitize_label_for_AggregateExpression <- function(label) {
    gsub("[/\\:*?\"<>|_]", "-", label)
}

# Create an ArgumentParser object
parser <- ArgumentParser(description = 'Integrate multiple scRNA-seq datasets into one Seurat object')
#parser$add_argument("--min_cells", type="integer", help="Min cells per celltype for diff express", default = 5)
parser$add_argument("--outdir", help="directory to write to", required = TRUE)
#parser$add_argument("--file_pref", help="file prefix", required = TRUE)
parser$add_argument("--diagnosis", help="diagnosis", required = TRUE)
# Parse the command-line arguments
args <- parser$parse_args()

print(args)

outdir <- file.path(args$outdir, args$diagnosis)
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
joined_integrated_data <- readRDS("rds_files/cytotrace_seurat.rds")
wanted_groups = c("fibroblasts", "pericyte")

just_subset <- subset(joined_integrated_data, cell_group %in% wanted_groups)
just_subset <- subset(just_subset, !is.na(CytoTRACE2_Potency))
just_subset$m_labels <- mapping$m_labels[match(just_subset$singleR.labels_fine, mapping$singleR.labels_fine)]
just_subset$new_cat_labels <- mapping$new_cat_labels[match(just_subset$singleR.labels_fine, mapping$singleR.labels_fine)]
just_subset <- subset(just_subset, diagnosis == args$diagnosis)

## We compute the Progeny activity scores and add them to our Seurat object
## as a new assay called Progeny. 
just_subset <- progeny(just_subset, scale=FALSE, organism="Human", top=500, perm=1, 
    return_assay = TRUE)

## We can now directly apply Seurat functions in our Progeny scores. 
## For instance, we scale the pathway activity scores. 
just_subset <- Seurat::ScaleData(just_subset, assay = "progeny") 

## We transform Progeny scores into a data frame to better handling the results
progeny_scores_df <- 
    as.data.frame(t(GetAssayData(just_subset, slot = "scale.data", 
        assay = "progeny"))) %>%
    rownames_to_column("Cell") %>%
    gather(Pathway, Activity, -Cell) 

## Extract cell metadata for the join
cell_metadata <- data.frame(
    Cell = colnames(just_subset),
    CellType = just_subset$CytoTRACE2_Potency,
    stringsAsFactors = FALSE
)

## We match Progeny scores with the cell clusters.
progeny_scores_df <- inner_join(progeny_scores_df, cell_metadata, by = "Cell")

filenm <- file.path(outdir, "progeny_scores_df.csv")
write.csv(progeny_scores_df, filenm)

## We summarize the Progeny scores by cellpopulation
summarized_progeny_scores <- progeny_scores_df %>% 
    group_by(Pathway, CellType) %>%
    summarise(avg = mean(Activity), std = sd(Activity))


## We prepare the data for the plot
summarized_progeny_scores_df <- summarized_progeny_scores %>%
    dplyr::select(-std) %>%   
    spread(Pathway, avg) %>%
    data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE) 

filenm <- file.path(outdir, "summarized_progeny_scores_df.csv")
write.csv(summarized_progeny_scores_df, filenm)

paletteLength = 100
myColor = colorRampPalette(c("Darkblue", "white","red"))(paletteLength)

progenyBreaks = c(seq(min(summarized_progeny_scores_df), 0, 
                      length.out=ceiling(paletteLength/2) + 1),
                  seq(max(summarized_progeny_scores_df)/paletteLength, 
                      max(summarized_progeny_scores_df), 
                      length.out=floor(paletteLength/2)))
progeny_hmap = pheatmap(t(summarized_progeny_scores_df[,-1]),fontsize=14, 
                        fontsize_row = 10, 
                        color=myColor, breaks = progenyBreaks, 
                        main = "PROGENy (500)", angle_col = 45,
                        treeheight_col = 0,  border_color = NA)

plot_path <- file.path(outdir, "progeny_hmap.pdf")
ggsave(plot_path, progeny_hmap, width = 10, height = 8)
