#!/bin/bash

source activate bcftools

for file in *vcf; do
  sample=$(grep 'tumor_sample' $file | sed 's/.*=//')
  echo $sample; parse_1_file=$(basename -s .vcf $file)_parsed.tsv
  echo -e "CHROM\tPOS\tREF\tALT\tVAF\tSYMBOL\tGene\tConsequence\tSIFT\tPolyPhen\tCondel\tLoFtool\tBLOSUM62\tProtein_position\tAmino_acids\tCodons" > $parse_1_file
  bcftools view -f 'PASS' -s $sample $file | bcftools +split-vep -f "%CHROM\t%POS\t%REF\t%ALT\t[ %AF]\t%SYMBOL\t%Gene\t%Consequence\t%SIFT\t%PolyPhen\t%Condel\t%LoFtool\t%BLOSUM62\t%Protein_position\t%Amino_acids\t%Codons\n" >> $parse_1_file
done


