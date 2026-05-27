configfile: "./ref.yaml"
include: "./mutect2_alt_bed_v4.6.1.0_no_normal.snakefile"
include: "./Strelka_no_normal.snakefile"

tumors = config["pairs"]


rule overlap_all:
    input:
        expand("easy_transfer_no_normal/overlap_no_normal/{tumor}_{f_score_thresh}_overlap_parsed.vcf", tumor = tumors, f_score_thresh = [0.8]),

rule compute_overlap:
    input:
        mutect2_vcf = rules.VEP.output.filtered_vcf,
        strelka_vcf = rules.Strelka_execute.output.vcf,
    output:
        normed_mutect2 = "overlap_no_normal/normed/{tumor}/mutect2_normed_{f_score_thresh}.vcf.gz",
        normed_strelka = "overlap_no_normal/normed/{tumor}/strelka_normed_{f_score_thresh}.vcf.gz",
        final_out_vcf = "easy_transfer_no_normal/overlap_no_normal/{tumor}_{f_score_thresh}_overlap.vcf",
        final_out_parsed = "easy_transfer_no_normal/overlap_no_normal/{tumor}_{f_score_thresh}_overlap_parsed.vcf",
    conda: "envs_dir/mutect2_strelka_overlap.yaml"
    params:
        ref_fasta = config["ref_fasta"],
        bioinfo_workflows_path = config["bioinfo_workflows_path"],
    log: "overlap_no_normal/logs/compute_overlap/{tumor}_{f_score_thresh}.log"
    benchmark: "overlap_no_normal/benchmarks/compute_overlap/{tumor}_{f_score_thresh}.benchmark"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
    shell:"""
echo "-------------------------------------------------------" > {log}
bcftools norm -f {params.ref_fasta} -Oz -o {output.normed_mutect2} {input.mutect2_vcf} >> {log} 2>&1
echo "-------------------------------------------------------" >> {log}
bcftools index {output.normed_mutect2} >> {log} 2>&1

echo "-------------------------------------------------------" >> {log}
bcftools norm -f {params.ref_fasta} -Oz -o {output.normed_strelka} {input.strelka_vcf} >> {log} 2>&1
echo "-------------------------------------------------------" >> {log}
bcftools index {output.normed_strelka} >> {log} 2>&1


echo "-------------------------------------------------------" >> {log}
bcftools isec -p overlap_no_normal/{wildcards.tumor}/ -f 'PASS' {output.normed_mutect2} {output.normed_strelka} >> {log} 2>&1

echo "-------------------------------------------------------" >> {log}
mv overlap_no_normal/{wildcards.tumor}/0002.vcf {output.final_out_vcf} >> {log} 2>&1
echo "-------------------------------------------------------" >> {log}
python3 {params.bioinfo_workflows_path}/scripts_dir/vep_vcf_parser.py \
    -f SYMBOL Gene Consequence SIFT PolyPhen Condel LoFtool BLOSUM62 Protein_position Amino_acids Codons \
    -v {output.final_out_vcf} -p > {output.final_out_parsed}

"""

