configfile: "{}/ref.yaml".format(workflow.basedir)

#include: "mutect2_alt_bed.snakefile"

tumors = config["pairs"]
#normals = [config["pairs"][tumor] for tumor in tumors]

rule signatureanalyzer_all:
  input:
    expand("signatureanalyzer/cosmic3_exome/{tumor}/signature_contributions.pdf", tumor = tumors),
#    expand("signatureanalyzer/cosmic3_DBS/{tumor}/signature_contributions.pdf", tumor = tumors),
    expand("signatureanalyzer/pcawg_SBS/{tumor}/signature_contributions.pdf", tumor = tumors),
    #expand("signatureanalyzer/cosmic3_ID/{tumor}/signature_contributions.pdf", tumor = tumors),

rule vcf2maf:
  input:
    # vcf = rules.VEP.output.vcf_out,
    # change to make cedar runs easier
    lambda wildcards: f"overlap_mutect2_vep_filtered/{wildcards.tumor}/{wildcards.tumor}.vcf"
  output:
    maf = "signatureanalyzer/{tumor}/{tumor}.maf",
  conda: "envs_dir/SignatureAnalyzer.yaml"
  log: "signatureanalyzer/log/{tumor}_vcf2maf.log"
  benchmark: "signatureanalyzer/benchmark/{tumor}_vcf2maf.benchmark"
  params:
    ref_fasta = config["ref_fasta"],
    normal = lambda wildcards: config['pairs'][wildcards.tumor],
  threads: 1
  resources:
    mem_mb = 4000
  shell:
    """
      vcf2maf.pl --input-vcf {input.vcf} --output-maf {output.maf} --ref-fasta {params.ref_fasta} --inhibit-vep --tumor-id {wildcards.tumor} --normal-id {params.normal} &> {log}
      sed -i '/^#/d' {output.maf} &>> {log}
    """

rule signatureanalyzer_cosmic3_exome:
  input:
    maf = rules.vcf2maf.output.maf,
  output:
    pdf = "signatureanalyzer/cosmic3_exome/{tumor}/signature_contributions.pdf",
  conda: "envs_dir/SignatureAnalyzer.yaml"
  log: "signatureanalyzer/log/{tumor}_signatureanalyzer_cosmic3_exome.log"
  benchmark: "signatureanalyzer/benchmark/{tumor}_signatureanalyzer_cosmic3_exome.benchmark"
  params:
    ref_hg_build = "ref/hg38.2bit",
    outdir = "signatureanalyzer/cosmic3_exome/{tumor}"
  threads: 1
  resources:
    mem_mb = lambda wildcards, attempt: attempt * 4000,
    time_min = lambda wildcards, attempt: attempt * 24 * 60,  # time in minutes
  shell:
    """
      signatureanalyzer -n 10 \
        --reference cosmic3_exome \
        --hg_build {params.ref_hg_build} \
        --objective poisson \
        --max_iter 30000 \
        --prior_on_H L1 \
        --prior_on_W L1 \
        -o {params.outdir} \
        {input.maf} &>> {log}
    """

rule signatureanalyzer_cosmic3_DBS:
  input:
    maf = rules.vcf2maf.output.maf,
  output:
    pdf = "signatureanalyzer/cosmic3_DBS/{tumor}/signature_contributions.pdf",
  conda: "envs_dir/SignatureAnalyzer.yaml"
  log: "signatureanalyzer/log/{tumor}_signatureanalyzer_cosmic3_DBS.log"
  benchmark: "signatureanalyzer/benchmark/{tumor}_signatureanalyzer_cosmic3_DBS.benchmark"
  params:
    ref_hg_build = "ref/hg38.2bit",
    outdir = "signatureanalyzer/cosmic3_DBS/{tumor}"
  threads: 1
  resources:
    mem_mb = lambda wildcards, attempt: attempt * 4000,
    time_min = lambda wildcards, attempt: attempt * 24 * 60,  # time in minutes
  shell:
    """
      signatureanalyzer -n 10 \
        --reference cosmic3_DBS \
        --hg_build {params.ref_hg_build} \
        --objective poisson \
        --max_iter 30000 \
        --prior_on_H L1 \
        --prior_on_W L1 \
        -o {params.outdir} \
        {input.maf} &>> {log}
    """

rule signatureanalyzer_pcawg_SBS:
  input:
    maf = rules.vcf2maf.output.maf,
  output:
    pdf = "signatureanalyzer/pcawg_SBS/{tumor}/signature_contributions.pdf",
  conda: "envs_dir/SignatureAnalyzer.yaml"
  log: "signatureanalyzer/log/{tumor}_signatureanalyzer_pcawg_SBS.log"
  benchmark: "signatureanalyzer/benchmark/{tumor}_signatureanalyzer_pcawg_SBS.benchmark"
  params:
    ref_hg_build = "ref/hg38.2bit",
    outdir = "signatureanalyzer/pcawg_SBS/{tumor}"
  threads: 1
  resources:
    mem_mb = lambda wildcards, attempt: attempt * 4000,
    time_min = lambda wildcards, attempt: attempt * 24 * 60,  # time in minutes
  shell:
    """
      signatureanalyzer -n 10 \
        --reference pcawg_SBS \
        --hg_build {params.ref_hg_build} \
        --objective poisson \
        --max_iter 30000 \
        --prior_on_H L1 \
        --prior_on_W L1 \
        -o {params.outdir} \
        {input.maf} &>> {log}
    """

rule signatureanalyzer_cosmic3_ID:
  input:
    maf = rules.vcf2maf.output.maf,
  output:
    pdf = "signatureanalyzer/cosmic3_ID/{tumor}/signature_contributions.pdf",
  conda: "envs_dir/SignatureAnalyzer.yaml"
  log: "signatureanalyzer/log/{tumor}_signatureanalyzer_cosmic3_ID.log"
  benchmark: "signatureanalyzer/benchmark/{tumor}_signatureanalyzer_cosmic3_ID.benchmark"
  params:
    ref_hg_build = "ref/hg38.2bit",
    outdir = "signatureanalyzer/cosmic3_ID/{tumor}"
  threads: 1
  resources:
    mem_mb = lambda wildcards, attempt: attempt * 4000,
    time_min = lambda wildcards, attempt: attempt * 24 * 60,  # time in minutes
  shell:
    """
      signatureanalyzer -n 10 \
        --reference cosmic3_ID \
        --hg_build {params.ref_hg_build} \
        --objective poisson \
        --max_iter 30000 \
        --prior_on_H L1 \
        --prior_on_W L1 \
        -o {params.outdir} \
        {input.maf} &>> {log}
    """

rule signatureanalyzer_polymerase_msi:
  input:
    maf = rules.vcf2maf.output.maf,
  output:
    pdf = "signatureanalyzer/polymerase_msi/{tumor}/signature_contributions.pdf",
  conda: "envs_dir/SignatureAnalyzer.yaml"
  log: "signatureanalyzer/log/{tumor}_signatureanalyzer_polymerase_msi.log"
  benchmark: "signatureanalyzer/benchmark/{tumor}_signatureanalyzer_polymerase_msi.benchmark"
  params:
    ref_hg_build = "ref/hg38.2bit",
    outdir = "signatureanalyzer/polymerase_msi/{tumor}"
  threads: 1
  resources:
    mem_mb = lambda wildcards, attempt: attempt * 4000,
    time_min = lambda wildcards, attempt: attempt * 24 * 60,  # time in minutes
  shell:
    """
      signatureanalyzer -n 10 \
        --reference polymerase_msi \
        --hg_build {params.ref_hg_build} \
        --objective poisson \
        --max_iter 30000 \
        --prior_on_H L1 \
        --prior_on_W L1 \
        -o {params.outdir} \
        {input.maf} &>> {log}
    """
