# Load libraries
library(argparse)

parser = ArgumentParser()
parser$add_argument('-f', '--filenm', required = TRUE, help = 'tsv file with transcript ids in first coloumn and no header')

args = parser$parse_args()
print(args)

library(readr)
library(biomaRt)

test_version = biomaRt::listEnsemblArchives()
print(test_version)
v98 = test_version[test_version$version == 98,]
print(v98$url)
if (v98$version != 98){
  stop("Biomart no longer using Ensembl Genes 98.  Make sure to update kallisto results to newest version of Ensembl Genes and this script.")
}
mart = biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = v98$url)
t2g = biomaRt::getBM(attributes = c("ensembl_transcript_id", "external_gene_name"), mart = mart)
tx2gene = dplyr::rename(t2g, target_id = ensembl_transcript_id, ens_gene = external_gene_name)

head(tx2gene)

transcripts = read.table(args$filenm)

head(transcripts$V1)

indices = sapply(transcripts$V1, function(x) match(x, tx2gene$target_id))
head(indices)

write.table(tx2gene[indices,c("target_id", "ens_gene")], file = "gene_names.tsv", quote = FALSE, sep = '\t', col.names = NA)
