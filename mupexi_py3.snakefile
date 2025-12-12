configfile: "/ref.yaml"
#include: "./mutect2_alt_bed_v4.6.1.0.snakefile"
#include: "./Strelka.snakefile"
include: "./hla-hd.snakefile"

tumors = config["pairs"]


rule mupexi_all:
    input:
        expand("mupexi/mupexi/{tumor}/{tumor}_{f_score_thresh}.mupexi", tumor = tumors, f_score_thresh = [0.8]),

rule compute_overlap_mupexi:
    input:
        mutect2_vcf = 'GATK_runs/{tumor}_{f_score_thresh}/Filter/{tumor}_{f_score_thresh}.vcf',
        strelka_snv_vcf = "Strelka_runs/{tumor}/results/variants/somatic.snvs.vcf.gz",
        strelka_indels_vcf = "Strelka_runs/{tumor}/results/variants/somatic.indels.vcf.gz",
#        strelka_snv_vcf = rules.Strelka_execute.output.vcfs_snvs,
#        strelka_indels_vcf = rules.Strelka_execute.output.vcfs_indels,
    output:
        normed_mutect2 = "mupexi/overlap/{tumor}/mutect2_normed_{f_score_thresh}.vcf.gz",
        normed_strelka = "mupexi/overlap/{tumor}/strelka_normed_{f_score_thresh}.vcf.gz",
        concat_strelka = "mupexi/overlap/{tumor}/strelka_concat_{f_score_thresh}.vcf.gz",
        sorted_strelka = "mupexi/overlap/{tumor}/strelka_sorted_{f_score_thresh}.vcf.gz",
        final_out_vcf = "mupexi/overlap/{tumor}/{tumor}_{f_score_thresh}.vcf"
    conda: "envs_dir/mupexi_py3.yaml"
    params:
        ref_fasta = config["ref_fasta"],
        bioinfo_workflows_path = config["bioinfo_workflows_path"],
    log: "mupexi/logs/compute_overlap/{tumor}_{f_score_thresh}.log"
    benchmark: "mupexi/benchmarks/compute_overlap/{tumor}_{f_score_thresh}.benchmark"
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
"""

rule mupexi_hla_1:
    input:
        vcf = rules.compute_overlap_mupexi.output.final_out_vcf,
        hlas = lambda wildcards: "hla-hd/" + config['merge_libs'][config['pairs'][wildcards.tumor]][0] + "/result/" + config['merge_libs'][config['pairs'][wildcards.tumor]][0] + "_final.result.txt"
    output:
        mupexi_out = 'mupexi/mupexi/{tumor}/{tumor}_{f_score_thresh}.mupexi',
        peptides = 'mupexi/mupexi/{tumor}/{tumor}_{f_score_thresh}_peptide_netMHCinput.txt',
        hla2s = 'mupexi/mupexi/{tumor}/{tumor}_{f_score_thresh}netMHCII.txt'
    conda: "envs_dir/mupexi_py3.yaml"
    log: "mupexi/logs/mupexi/{tumor}_{f_score_thresh}.log"
    benchmark: "mupexi/benchmarks/mupexi/{tumor}_{f_score_thresh}.benchmark"
    threads: 1
    params:
        scripts_dir = config["bioinfo_workflows_path"] + '/scripts_dir/',
        hla1 = 'mupexi/mupexi/{tumor}/{tumor}_{f_score_thresh}_hla1.txt',
        hla2 = 'mupexi/mupexi/{tumor}/{tumor}_{f_score_thresh}_hla2.txt',
        netMHCIIpan_path = '/home/dylan/netmhcs/netMHCIIpan-4.3/netMHCIIpan'
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
    shell:"""
python ~/bioinfo_workflows/scripts_dir/get_hlas_from_hla_hd.py {input.hlas} {params.hla1} {params.hla2}

python ~/githubs/MuPeXI/MuPeXI_mini.py -v {input.vcf}\
    -a $(cat {params.hla1})\
    -c ~/githubs/MuPeXI/config.ini\
    -t -d mupexi/mupexi/{wildcards.tumor}/\
    -p {wildcards.tumor}_{wildcards.f_score_thresh}\
    -l 9\
    -f &>> {log}

 {params.netMHCIIpan_path} -inptype 1 -f {output.peptides} -a $(cat {params.hla2}) 1> {output.hla2s} 2>> {log}
"""
