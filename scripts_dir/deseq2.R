# Load libraries
library(argparse)
library(DESeq2)
library(tximportData)
library(tximport)
library(readr)
library(biomaRt)
library(EnhancedVolcano)

parser = ArgumentParser()
parser$add_argument('-m', '--sample_mapping', required = TRUE, help = 'csv file containing sample to condition to path mapping')
parser$add_argument('-g', '--gene', action = 'store_true', help = 'use gene names instead of gene ids')
parser$add_argument('-r', '--min_reads', type = 'integer', default = 5, help = 'min reads to pass filter')
parser$add_argument('-p', '--min_proportion', type = 'double', default = 0.47, help = 'min proportion of samples having min reads to pass filter')
parser$add_argument('-o', '--out_pre', default = "", help = 'prefix of output files')
parser$add_argument('-O', '--out_dir', default = ".", help = 'directory to write files to')
parser$add_argument('-R', '--rsem', action = 'store_true', help = 'run in rsem mode')
parser$add_argument('--ref', default = 'control', help = 'which condition is the control')

args = parser$parse_args()
print(args)


test_version = biomaRt::listEnsemblArchives()
print(test_version)
v98 = test_version[test_version$version == 98,]
print(v98$url)
if (v98$version != 98){
  stop("Biomart no longer using Ensembl Genes 98.  Make sure to update kallisto results to newest version of Ensembl Genes and this script.")
}
mart = biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = v98$url)
if (args$gene){
  t2g = biomaRt::getBM(attributes = c("ensembl_transcript_id", "external_gene_name"), mart = mart)
  tx2gene = dplyr::rename(t2g, target_id = ensembl_transcript_id, ens_gene = external_gene_name)
} else{
  t2g = biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id"), mart = mart)
  tx2gene = dplyr::rename(t2g, target_id = ensembl_transcript_id, ens_gene = ensembl_gene_id)
}

head(tx2gene)

# CUSTOM TX2GENE FROM GRCh38.p13 #

if (!dir.exists(args$out_dir)){
  dir.create(args$out_dir)
}

#tx2gene = read.table(args$tran_2_gene, stringsAsFactors = F)

coldata = read.csv(args$sample_mapping)
num_samples = length(coldata[, 1])
conditions = sort(unique(coldata$condition))
conditions

# Import data
if (args$rsem){
    coldata$path = paste(coldata$path, ".isoforms.results", sep = "")
    txi.tsv = tximport(coldata$path, type = "rsem", tx2gene = tx2gene, ignoreTxVersion = TRUE)
} else{
    coldata$path = paste(coldata$path, "abundance.h5", sep = "")
    txi.tsv = tximport(coldata$path, type = "kallisto", tx2gene = tx2gene, ignoreTxVersion = TRUE)
}

dds = DESeqDataSetFromTximport(txi.tsv, colData = coldata, design = ~condition)
# filter out lowly expressed genes
dds = estimateSizeFactors(dds)
#count_data = counts(dds)
#colnames(count_data) = coldata$sample
#write.csv(count_data, file = "counts_cohort.csv")
#stop()
kept_idx = rowSums(counts(dds, normalized = TRUE) >= args$min_reads) >= (num_samples * args$min_proportion)

dds = DESeq(dds[kept_idx, ])

# make the condition in args$ref the reference
dds$condition = relevel(dds$condition, ref = args$ref)

res = results(dds)
print(res)
res = as.data.frame(res)
res = res[order(res$padj),]
write.csv(res, paste(args$out_dir, "/", args$out_pre, "_deseq2.csv", sep = ""), row.names = TRUE)

plot = EnhancedVolcano(res,
  lab = rownames(res),
  x = 'log2FoldChange',
  y = 'padj',
  pCutoff = 0.05,
  FCcutoff = 1)
ggsave(paste(args$out_dir, "/", args$out_pre, "_deseq2_volcano.pdf", sep = ""), )

plot = EnhancedVolcano(res,
  lab = rownames(res),
  x = 'log2FoldChange',
  y = 'padj',
  ylim = c(0, -log10(10e-20)),
  pCutoff = 0.05,
  FCcutoff = 1)
ggsave(paste(args$out_dir, "/", args$out_pre, "_deseq2_volcano_max20.pdf", sep = ""), )