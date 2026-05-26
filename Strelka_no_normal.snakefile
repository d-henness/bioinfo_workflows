include: "pre_pro_af_merge_alt_bed.snakefile"
configfile: "{}/ref.yaml".format(workflow.basedir)

rule run_Strelka:
  input:
    expand("Strelka_runs_no_normal/{tumor}/results/variants/variants.vcf.gz", tumor = config["pairs"])

rule Strelka_config:
  input:
    tumor_bam = "GATK_runs/{tumor}/ApplyBQSR/{tumor}_recal.bam",
#    normal_bam = lambda wildcards: "GATK_runs/" + config["pairs"][wildcards.tumor] + "/ApplyBQSR/" + config["pairs"][wildcards.tumor] + "_recal.bam"
  output:
    out_run = "Strelka_runs_no_normal/{tumor}/runWorkflow.py"
  conda:
    "envs_dir/Strelka_env.yaml"
  params:
    ref_fasta = config["ref_fasta"],
    bed_file = config['alt_bed'] + ".gz",
  threads: 1
  resources:
    mem_mb = 1000
  benchmark:
    "benchmarks/Strelka_config_no_normal.{tumor}.txt"
  log:
    config_log = "Strelka_runs_no_normal/{tumor}/config.log",
  shell:
    """
      configureStrelkaGermlineWorkflow.py  \
        --bam {input.tumor_bam}  \
        --ref {params.ref_fasta}  \
        --exome \
        --callRegions {params.bed_file} \
        --runDir Strelka_runs_no_normal/{wildcards.tumor} &> {log.config_log}
    """
rule Strelka_execute:
  input:
    rules.Strelka_config.output.out_run
  output:
    vcf = "Strelka_runs_no_normal/{tumor}/results/variants/variants.vcf.gz",
  conda:
    "envs_dir/Strelka_env.yaml"
  threads: 1
  resources:
    mem_mb = 4000,
    mem_gb = 4
  benchmark:
    "benchmarks/Strelka_execute_no_normal.{tumor}.txt"
  shell:
    """
      {input} --mode local --jobs {threads} --memGb {resources.mem_gb}
    """
