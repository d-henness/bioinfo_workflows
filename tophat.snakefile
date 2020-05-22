configfile: "{}/ref.yaml".format(workflow.basedir)
include: "fastp_rna.snakefile"

rule tophat_all:
  input:
    expand("fastp/{rna_lib}/{rna_lib}_signal.txt", rna_lib = config["rna_merge_libs"]),

rule tophat:
  input:
    fq1 = lambda wildcards: f"fastp/{config['rna_merge_libs'][wildcards.rna_lib][0]}/fastp/{config['rna_merge_libs'][wildcards.rna_lib][0]}_1.fq.gz",
    fq2 = lambda wildcards: f"fastp/{config['rna_merge_libs'][wildcards.rna_lib][0]}/fastp/{config['rna_merge_libs'][wildcards.rna_lib][0]}_2.fq.gz",
  output:
    signal = "fastp/{rna_lib}/{rna_lib}_signal.txt",
  conda:
    "envs_dir/tophat.yaml"
  resources:
    mem_mb = lambda wildcards, attempt: attempt * (6 * 1024),
    time_min = lambda wildcards, attempt: attempt * (1 * 60),	# time in minutes
  threads: 4
  params:
    index = config["bowtie2_cdna"],
    out_dir = lambda wildcards: f"tophat/{wildcards.rna_lib}"
  benchmark:
    "tophat/benchmark/{rna_lib}_tophat.benchmark"
  log:
    "tophat/logs/{rna_lib}_tophat.log"
  shell:
    """
        tophat \
            -o {params.out_dir} \
            -p {threads} \
            {params.index} \
            {input.fq1} {input.fq2} &> {log}
        echo "test" > {output.signal}
    """
