configfile: "./ref.yaml"
include: "./mutect2_alt_bed_v4.6.1.0.snakefile"
include: "./Strelka.snakefile"

tumors = config["pairs"]


rule overlap_all:
    input:
        expand("easy_transfer/overlap/{tumor}_{f_score_thresh}_overlap_parsed.tsv", tumor = tumors, f_score_thresh = [0.8]),

rule compute_overlap:
    input:
        mutect2_vcf = rules.VEP.output.vcf_out,
        strelka_snv_vcf = rules.Strelka_execute.output.vcfs_snvs,
        strelka_indels_vcf = rules.Strelka_execute.output.vcfs_indels,
    output:
        normed_mutect2 = "overlap/normed/{tumor}/mutect2_normed_{f_score_thresh}.vcf.gz",
        normed_strelka = "overlap/normed/{tumor}/strelka_normed_{f_score_thresh}.vcf.gz",
        concat_strelka = "overlap/normed/{tumor}/strelka_concat_{f_score_thresh}.vcf.gz",
        sorted_strelka = "overlap/normed/{tumor}/strelka_sorted_{f_score_thresh}.vcf.gz",
        final_out_vcf = "easy_transfer/overlap/{tumor}_{f_score_thresh}_overlap.vcf",
        final_out_parsed = "easy_transfer/overlap/{tumor}_{f_score_thresh}_overlap_parsed.tsv",
    conda: "envs_dir/mutect2_strelka_overlap.yaml"
    params:
        ref_fasta = config["ref_fasta"],
        bioinfo_workflows_path = config["bioinfo_workflows_path"],
    log: "overlap/logs/compute_overlap/{tumor}_{f_score_thresh}.log"
    benchmark: "overlap/benchmarks/compute_overlap/{tumor}_{f_score_thresh}.benchmark"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
    shell:"""
echo "-------------------------------------------------------" > {log}
bcftools norm -f {params.ref_fasta} -Oz -o {output.normed_mutect2} {input.mutect2_vcf} >> {log} 2>&1
echo "-------------------------------------------------------" >> {log}
bcftools index {output.normed_mutect2} >> {log} 2>&1

echo "-------------------------------------------------------" >> {log}
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
bcftools isec -p overlap/{wildcards.tumor}/ -f 'PASS' {output.normed_mutect2} {output.sorted_strelka} >> {log} 2>&1

echo "-------------------------------------------------------" >> {log}
mv overlap/{wildcards.tumor}/0002.vcf {output.final_out_vcf} >> {log} 2>&1
echo "-------------------------------------------------------" >> {log}
python3 {params.bioinfo_workflows_path}/scripts_dir/vep_vcf_parser.py \
    -f SYMBOL Gene Consequence SIFT PolyPhen Condel LoFtool BLOSUM62 Protein_position Amino_acids Codons \
    -v {output.final_out_vcf} -p > {output.final_out_parsed}

"""

