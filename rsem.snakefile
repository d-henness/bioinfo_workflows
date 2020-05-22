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
    genes = "rsem/{wildcards.rna_lib}/{wildcards.rna_lib}.genes.results",
    isoforms = "rsem/{wildcards.rna_lib}/{wildcards.rna_lib}.isoforms.results",
  conda:
    "envs_dir/rsem.yaml"
  resources:
    mem_mb = lambda wildcards, attempt: attempt * (4 * 1024),
    time_min = lambda wildcards, attempt: attempt * (1 * 60),	# time in minutes
  threads: 4
  params:
    index = config["rsem_index"],
    out_pre = lambda wildcards: f"rsem/{wildcards.rna_lib}/{wildcards.rna_lib}"
  benchmark:
    "rsem/benchmark/{rna_lib}_rsem.benchmark"
  log:
    "rsem/logs/{rna_lib}_rsem.log"
  shell:
    """
    rsem-calculate-expression -p 8 \
        --paired-end \
        --bowtie2 --bowtie2-path $CONDA_PREFIX/bin \
        --estimate-rspd \
        --append-names \
        {input.fq1} {input.fq2}  \
        {params.index} \
        {params.out_pre} &> {log}
    """
