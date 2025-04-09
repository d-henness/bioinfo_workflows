library(clusterProfiler)
library(argparse)
library(org.Hs.eg.db)
library(readr)
library(ggplot2)
library(GO.db)
library(GSEABase)
library(DOSE)
library(dplyr)
library(tidyr)

parser <- ArgumentParser()
parser$add_argument("filenm")
parser$add_argument("-o", "--output_pref", default = "test")
args <- parser$parse_args()

diff_express_data <- read.csv(args$filenm, stringsAsFactors = FALSE)

nrow(diff_express_data)
head(diff_express_data)

up.genes <- diff_express_data[diff_express_data$log2FoldChange > 1 & diff_express_data$padj < 0.05, 1]
down.genes <- diff_express_data[diff_express_data$log2FoldChange < -1 & diff_express_data$padj < 0.05, 1]
bkgd.genes <- diff_express_data[,1]

up.genes.entrez <- clusterProfiler::bitr(up.genes,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)

print("up regged genes")
nrow(up.genes.entrez)
head(up.genes.entrez)

keytypes(org.Hs.eg.db)
down.genes.entrez <- bitr(down.genes,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
bkgd.genes.entrez <- bitr(bkgd.genes,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)

print("down regged genes")
nrow(down.genes.entrez)
head(down.genes.entrez)

print("background genes")
nrow(bkgd.genes.entrez)
head(bkgd.genes.entrez)

egobp_up <- clusterProfiler::enrichGO(
        gene     = up.genes.entrez[[2]],
        universe = bkgd.genes.entrez[[2]],
        OrgDb    = org.Hs.eg.db,
        ont      = "BP",
        pAdjustMethod = "fdr",
        pvalueCutoff = 0.05, #p.adjust cutoff (https://github.com/GuangchuangYu/clusterProfiler/issues/104)
        readable = TRUE)

head(egobp_up,10)


barplot(egobp_up, showCategory = 20)
ggsave(paste0(args$output_pref, "_up_genes_barplot.png"))
dotplot(egobp_up, showCategory = 20)
ggsave(paste0(args$output_pref, "_up_genes_dotplot.png"))
goplot(egobp_up)
ggsave(paste0(args$output_pref, "_up_genes_GOplot.png"))


egobp_down <- clusterProfiler::enrichGO(
        gene     = down.genes.entrez[[2]],
        universe = bkgd.genes.entrez[[2]],
        OrgDb    = org.Hs.eg.db,
        ont      = "BP",
        pAdjustMethod = "fdr",
        pvalueCutoff = 0.05, #p.adjust cutoff (https://github.com/GuangchuangYu/clusterProfiler/issues/104)
        readable = TRUE)

head(egobp_down,10)

barplot(egobp_down, showCategory = 20)
ggsave(paste0(args$output_pref, "_down_genes_barplot.png"))
dotplot(egobp_down, showCategory = 20)
ggsave(paste0(args$output_pref, "_down_genes_dotplot.png"))
goplot(egobp_down)
ggsave(paste0(args$output_pref, "_down_genes_GOplot.png"))
