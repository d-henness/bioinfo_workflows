include: "fastp_rna.snakefile"

configfile: "{}/ref.yaml".format(workflow.basedir)

rule rna_variants_all:
  input:
#    expand("GATK_runs/{library}/MergeBamAlignment/merge.bam", library = config["libraries"])
    expand("GATK_runs/{library}/ubam/unprocessed.bam", library = config["libraries"])

def get_sample_name(wildcards):
  for sample_name in config["merge_libs"]:
    for lib in config["merge_libs"][sample_name]:
      if lib == config["libraries"][wildcards.library]:
        print(sample_name)

rule ubam:
  input:
    fq1 = rules.fastp_paired.output.fq1_out,
    fq2 = rules.fastp_paired.output.fq2_out,
  output:
    "rna_GATK_runs/{library}/ubam/{library}_unprocessed.bam",
  conda:
    "envs_dir/pre_proc.yaml"
  resources:
    mem_mb = lambda wildcards, attempt: attempt * 5000,
    time_min = lambda wildcards, attempt: attempt * 24 * 60,	# time in minutes
    io = 1, # used to indicate that this is an io heavy job and should not have many running at once
  params:
    exclude_list = '',
    sample_name = get_sample_name
  log:
    ubam = "rna_GATK_runs/log/{library}_ubam.log"
  benchmark:
    "rna_GATK_runs/benchmarks/{library}_ubam.benchmark.txt"
  shell:
    "python3 {workflow.basedir}/scripts_dir/ubam.py {input} {wildcards.library} {output} 2> {log.ubam} {params.sample_name}"

rule star:
  input:
    fq1 = rules.fastp_paired.output.fq1_out,
    fq2 = rules.fastp_paired.output.fq2_out,
  output:
    "rna_GATK_runs/{library}/star/{library}_unprocessed.bam",
  conda:
    "envs_dir/pre_proc.yaml"
  threads: 8
  resources:
    mem_mb = lambda wildcards, attempt: attempt * 10000,
    time_min = lambda wildcards, attempt: attempt * 24 * 60,	# time in minutes
    io = 1, # used to indicate that this is an io heavy job and should not have many running at once
  params:
  log:
    ubam = "rna_GATK_runs/log/{library}_star.log"
  benchmark:
    "rna_GATK_runs/benchmarks/{library}_star.benchmark.txt"
  shell:
    "python3 {workflow.basedir}/scripts_dir/ubam.py {input} {wildcards.library} {output} 2> {log.ubam} {params.sample_name}"
