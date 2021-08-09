configfile: "{}/ref.yaml".format(workflow.basedir)
include: "fastp_rna.snakefile"

rule rsem_all:
  input:
    expand("rsem/{rna_lib}/{rna_lib}.genes.results", rna_lib = config["rna_merge_libs"]),
    expand("rsem/{rna_lib}/{rna_lib}.transcript-sorted.bam", rna_lib = config["rna_merge_libs"]),
    expand("rsem/{rna_lib}/{rna_lib}.transcript-sorted.bam.bai", rna_lib = config["rna_merge_libs"]),

rule rsem:
  input:
    fq1 = lambda wildcards: f"fastp/{config['rna_merge_libs'][wildcards.rna_lib][0]}/fastp/{config['rna_merge_libs'][wildcards.rna_lib][0]}_1.fq.gz",
    fq2 = lambda wildcards: f"fastp/{config['rna_merge_libs'][wildcards.rna_lib][0]}/fastp/{config['rna_merge_libs'][wildcards.rna_lib][0]}_2.fq.gz",
  output:
    genes = "rsem/{rna_lib}/{rna_lib}.genes.results",
    isoforms = "rsem/{rna_lib}/{rna_lib}.isoforms.results",
    bam = "rsem/{rna_lib}/{rna_lib}.transcript.bam",
  conda:
    "envs_dir/rsem.yaml"
  resources:
    mem_mb = lambda wildcards, attempt: attempt * (4 * 1024),
    time_min = lambda wildcards, attempt: attempt * (24 * 60),	# time in minutes
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
    rsem-calculate-expression -p {threads} \
        --paired-end \
        --bowtie2 --bowtie2-path $CONDA_PREFIX/bin \
        --estimate-rspd \
        {input.fq1} {input.fq2}  \
        {params.index} \
        {params.out_pre} &> {log}
    """

rule sort_bam:
  input:
    bam = rules.rsem.output.bam
  output:
    bam = "rsem/{rna_lib}/{rna_lib}.transcript-sorted.bam",
  conda:
    "envs_dir/rsem.yaml"
  resources:
    mem_mb = lambda wildcards, attempt: attempt * (4 * 1024),
    time_min = lambda wildcards, attempt: attempt * (24 * 60),	# time in minutes
  threads: 4
  params:
    temp_dir = "/tmp/$SLURM_JOB_ID",
  benchmark:
    "rsem/benchmark/{rna_lib}_sort_bam.benchmark"
  log:
    "rsem/logs/{rna_lib}_sort_bam.log"
  shell:
    """
      if [[ !(-d {params.temp_dir}) ]]; then
        mkdir {params.temp_dir} &> {log}
      fi
      samtools sort -@ {threads} -O BAM -o {output.bam} -T {params.temp_dir}/tmp {input.bam} &>> {log}
    """

rule index_bam:
  input:
    bam = rules.sort_bam.output.bam
  output:
    bai = "rsem/{rna_lib}/{rna_lib}.transcript-sorted.bam.bai",
  conda:
    "envs_dir/rsem.yaml"
  resources:
    mem_mb = lambda wildcards, attempt: attempt * (4 * 1024),
    time_min = lambda wildcards, attempt: attempt * (24 * 60),	# time in minutes
  threads: 4
  benchmark:
    "rsem/benchmark/{rna_lib}_index_bam.benchmark"
  log:
    "rsem/logs/{rna_lib}_index_bam.log"
  shell:
    """
      samtools index -@ {threads} {input.bam} &>> {log}
    """
