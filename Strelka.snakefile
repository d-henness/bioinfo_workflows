configfile: "{}/ref.yaml".format(workflow.basedir)

rule run_Strelka:
  input:
    expand("Strelka_runs/{tumor}/results/variants/somatic.snvs.vcf.gz", tumor = config["pairs"])

rule Strelka_config:
  input:
    tumor_bam = "runs/{tumor}/ApplyBQSR/recal.bam",
    normal_bam = lambda wildcards: "runs/" + config["pairs"][wildcards.tumor] + "/ApplyBQSR/recal.bam"
  output:
    out_run = "Strelka_runs/{tumor}/runWorkflow.py"
  conda:
    "envs_dir/Strelka_env.yaml"
  params:
    ref_fasta = config["ref_fasta"],
    bed_file = "/data/synapse1/shared/gatk4/gatk_data/hg38/remapped_S07604624_Padded_primaryOnly.bed.gz"
  threads: 1
  resources:
    mem_mb = 1000
  benchmark:
    "benchmarks/Strelka_config.{tumor}.txt"
  log:
    config_log = "Strelka_runs/{tumor}/config.log",
  shell:
    """
      configureStrelkaSomaticWorkflow.py  \
        --normalBam {input.normal_bam}  \
        --tumorBam {input.tumor_bam}  \
        --ref {params.ref_fasta}  \
        --exome \
        --callRegions {params.bed_file} \
        --runDir Strelka_runs/{wildcards.tumor} &> {log.config_log}
    """
rule Strelka_execute:
  input:
    rules.Strelka_config.output.out_run
  output:
    vcfs_snvs = "Strelka_runs/{tumor}/results/variants/somatic.snvs.vcf.gz",
    vcfs_indels = "Strelka_runs/{tumor}/results/variants/somatic.indels.vcf.gz",
  conda:
    "envs_dir/Strelka_env.yaml"
  threads: 16
  resources:
    mem_mb = 4000,
    mem_gb = 4
  benchmark:
    "benchmarks/Strelka_execute.{tumor}.txt"
  shell:
    """
      {input} --mode local --jobs {threads} --memGb {resources.mem_gb}
    """
