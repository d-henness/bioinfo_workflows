configfile: "./ref.yaml"
include: "./mutect2_alt_bed_v4.6.1.0.snakefile"
include: "./Strelka.snakefile"

tumors = config["pairs"]


rule overlap_all:
    input:
        expand("easy_transfer/overlap/{tumor}_{f_score_thresh}_overlap_parsed.tsv", tumor = tumors, f_score_thresh = [0.8]),

rule compute_overlap:
    input:
        mutect2_vcf = rules.VEP.output.filtered_vcf,
        strelka_snv_vcf = rules.Strelka_execute.output.vcfs_snvs,
        strelka_indels_vcf = rules.Strelka_execute.output.vcfs_indels,
    output:
        normed_strelka = "overlap/normed/{tumor}/strelka_normed_{f_score_thresh}.vcf.gz",
        concat_strelka = "overlap/normed/{tumor}/strelka_concat_{f_score_thresh}.vcf.gz",
        sorted_strelka = "overlap/normed/{tumor}/strelka_sorted_{f_score_thresh}.vcf.gz",
        final_out_vcf = "easy_transfer/overlap/{tumor}_{f_score_thresh}_overlap.vcf",
        final_out_parsed = "easy_transfer/overlap/{tumor}_{f_score_thresh}_overlap_parsed.tsv",
    conda: "envs_dir/mutect2_strelka_overlap.yaml"
    params:
        ref_fasta = config["ref_fasta"],
    log: "overlap/logs/compute_overlap/{tumor}_{f_score_thresh}.log"
    benchmark: "overlap/benchmarks/compute_overlap/{tumor}_{f_score_thresh}.benchmark"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
    shell: r"""
echo "-------------------------------------------------------" > {log}
bcftools norm -f {params.ref_fasta} -Oz -o {output.normed_strelka} {input.strelka_indels_vcf} >> {log} 2>&1
echo "-------------------------------------------------------" >> {log}
bcftools index {output.normed_strelka} >> {log} 2>&1
echo "-------------------------------------------------------" >> {log}
bcftools concat -a -Oz -o {output.concat_strelka} {input.strelka_snv_vcf} {output.normed_strelka} >> {log} 2>&1
echo "-------------------------------------------------------" >> {log}
bcftools index {output.concat_strelka} >> {log} 2>&1
echo "-------------------------------------------------------" >> {log}
bcftools sort -Oz -o {output.sorted_strelka}  {output.concat_strelka} >> {log} 2>&1
echo "-------------------------------------------------------" >> {log}
bcftools index {output.sorted_strelka} >> {log} 2>&1

echo "-------------------------------------------------------" >> {log}
bcftools isec -p overlap/{wildcards.tumor}/ -f 'PASS' {input.mutect2_vcf} {output.sorted_strelka} >> {log} 2>&1

echo "-------------------------------------------------------" >> {log}
mv overlap/{wildcards.tumor}/0002.vcf {output.final_out_vcf} >> {log} 2>&1
echo "-------------------------------------------------------" >> {log}
TUMOR=$(bcftools view -h {output.final_out_vcf} | grep '^##tumor_sample=' | sed 's/.*=//')
printf 'chr\tpos\tref\talt\tvaf\tgnomAD_AF_joint\tgnomAD_AF_grpmax_joint\tSYMBOL\tGene\tConsequence\tSIFT\tPolyPhen\tCondel\tLoFtool\tBLOSUM62\tProtein_position\tAmino_acids\tCodons\n' > {output.final_out_parsed}
bcftools +split-vep -s "$TUMOR" \
    -f '%CHROM\t%POS\t%REF\t%ALT\t[%AF]\t%INFO/gnomADj_AF_joint\t%INFO/gnomADj_AF_grpmax_joint\t%SYMBOL\t%Gene\t%Consequence\t%SIFT\t%PolyPhen\t%Condel\t%LoFtool\t%BLOSUM62\t%Protein_position\t%Amino_acids\t%Codons\n' \
    -A tab \
    {output.final_out_vcf} >> {output.final_out_parsed} 2>> {log}

"""

