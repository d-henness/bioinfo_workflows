configfile: "{}/ref.yaml".format(workflow.basedir)
include: "fastp_rna.snakefile"

rule find_circ_all:
  input:
    expand("find_circ/find_circ/{rna_lib}/{rna_lib}_splice_sites.bed", rna_lib = config["rna_merge_libs"]),

rule bowtie2_1:
  input:
    fq1 = lambda wildcards: f"fastp/{config['rna_merge_libs'][wildcards.rna_lib][0]}/fastp/{config['rna_merge_libs'][wildcards.rna_lib][0]}_1.fq.gz",
    fq2 = lambda wildcards: f"fastp/{config['rna_merge_libs'][wildcards.rna_lib][0]}/fastp/{config['rna_merge_libs'][wildcards.rna_lib][0]}_2.fq.gz",
  output:
    bam = "find_circ/bowtie2_1/{rna_lib}/{rna_lib}.bam",
  conda:
    "envs_dir/find_circ.yaml"
  resources:
    mem_mb = lambda wildcards, attempt: attempt * (4 * 1024),
    time_min = lambda wildcards, attempt: attempt * (24 * 60),	# time in minutes
  threads: 4
  params:
    index = config["rsem_index"],
    out_pre = lambda wildcards: f"rsem/{wildcards.rna_lib}/{wildcards.rna_lib}"
  benchmark:
    "find_circ/benchmark/bowtie2_1_{rna_lib}.log"
  log:
    "find_circ/log/bowtie2_1_{rna_lib}.log"
  shell:
    """
      bowtie2 -p{threads} \
      --very-sensitive \
      --score-min=C,-15,0 \
      --mm \
      -x {params.index} \
      -q \
      -U {input.fq1} 2> {log}  \
      | samtools view -hbuS - | samtools sort -o {output.bam} -
    """

rule unmapped2anchors:
  input:
    bam = rules.bowtie2_1.output.bam
  output:
    unmapped_bam = "find_circ/unmapped2anchors/{rna_lib}/{rna_lib}_unmapped.bam",
    fq_anchors = "find_circ/unmapped2anchors/{rna_lib}/{rna_lib}_anchors.fq.gz"
  conda:
    "envs_dir/find_circ.yaml"
  resources:
    mem_mb = lambda wildcards, attempt: attempt * (4 * 1024),
    time_min = lambda wildcards, attempt: attempt * (24 * 60),  # time in minutes
  threads: 1
  params:
    out_pre = lambda wildcards: f"rsem/{wildcards.rna_lib}/{wildcards.rna_lib}",
    path_to_find_circ = "/home/dylan/find_circ_test/github/find_circ_quick_dirty"
  benchmark:
    "find_circ/benchmark/unmapped2anchors_{rna_lib}.log"
  log:
    "find_circ/log/unmapped2anchors_{rna_lib}.log"
  shell:
    """
      samtools view -hf 4 {input.bam} | samtools view -Sb - 1> {output.unmapped_bam} 2> {log}
      python {params.path_to_find_circ}/unmapped2anchors.py {output.unmapped_bam} 2>> {log} | gzip > {output.fq_anchors}
    """

rule bowtie2_2:
  input:
    fq = rules.unmapped2anchors.output.fq_anchors,
  output:
    bam = "find_circ/find_circ/{rna_lib}/{rna_lib}_anchors.bam",
  conda:
    "envs_dir/find_circ.yaml"
  resources:
    mem_mb = lambda wildcards, attempt: attempt * (4 * 1024),
    time_min = lambda wildcards, attempt: attempt * (24 * 60),	# time in minutes
  threads: 4
  params:
    index = "/data/shared/hg38/bowtie_index_find_circ/Homo_sapiens_assembly38",
    path_to_find_circ = "/home/dylan/find_circ_test/github/find_circ",
    out_dir = "find_circ/find_circ/{rna_lib}",
    fa = config["ref_fasta"],
  benchmark:
    "find_circ/benchmark/bowtie2_2_{rna_lib}.log"
  log:
    "find_circ/log/bowtie2_2_{rna_lib}.log"
  shell:
    """
      if [[ !(-d {params.out_dir}/) ]]; then
        mkdir {params.out_dir} &> {log}
      fi

      bowtie2 -p {threads} \
      --score-min=C,-15,0 \
      --reorder \
      --mm \
      -x {params.index} \
      -q \
      -U {input.fq} 2> {log} | samtools view -hbo {output.bam} -

    """

rule find_circ:
  input:
    bam = rules.bowtie2_2.output.bam,
  output:
    bed = "find_circ/find_circ/{rna_lib}/{rna_lib}_splice_sites.bed",
  conda:
    "envs_dir/find_circ.yaml"
  resources:
    mem_mb = lambda wildcards, attempt: attempt * (4 * 1024),
    time_min = lambda wildcards, attempt: attempt * (24 * 60),	# time in minutes
  threads: 4
  params:
    index = "/data/shared/hg38/bowtie_index_find_circ/Homo_sapiens_assembly38",
    path_to_find_circ = "/home/dylan/find_circ_test/github/find_circ",
    out_dir = "find_circ/find_circ/{rna_lib}",
    fa = config["ref_fasta"],
  benchmark:
    "find_circ/benchmark/find_circ_{rna_lib}.log"
  log:
    "find_circ/log/find_circ_{rna_lib}.log"
  shell:
    """
      samtools view -h {input.bam} | python {params.path_to_find_circ}/find_circ.py \
      -G {params.fa} \
      -p {wildcards.rna_lib}_ \
      -n {wildcards.rna_lib} \
      -s {params.out_dir}/{wildcards.rna_lib}_stats.txt \
      -R {params.out_dir}/{wildcards.rna_lib}_spliced_reads.fa 1> {output.bed} 2>> {log}
    """


    #rule sort_bam:
    #  input:
    #    bam = rules.rsem.output.bam
    #  output:
    #    bam = "rsem/{rna_lib}/{rna_lib}.transcript-sorted.bam",
    #  conda:
    #    "envs_dir/rsem.yaml"
    #  resources:
    #    mem_mb = lambda wildcards, attempt: attempt * (4 * 1024),
    #    time_min = lambda wildcards, attempt: attempt * (24 * 60),	# time in minutes
    #  threads: 4
    #  params:
    #    temp_dir = "/tmp/$SLURM_JOB_ID",
    #  benchmark:
    #    "rsem/benchmark/{rna_lib}_sort_bam.benchmark"
    #  log:
    #    "rsem/logs/{rna_lib}_sort_bam.log"
    #  shell:
    #    """
    #      if [[ !(-d {params.temp_dir}) ]]; then
    #        mkdir {params.temp_dir} &> {log}
    #      fi
    #      samtools sort -@ {threads} -O BAM -o {output.bam} -T {params.temp_dir}/tmp {input.bam} &>> {log}
    #    """
    #
    #rule index_bam:
    #  input:
    #    bam = rules.sort_bam.output.bam
    #  output:
    #    bai = "rsem/{rna_lib}/{rna_lib}.transcript-sorted.bam.bai",
    #  conda:
    #    "envs_dir/rsem.yaml"
    #  resources:
    #    mem_mb = lambda wildcards, attempt: attempt * (4 * 1024),
    #    time_min = lambda wildcards, attempt: attempt * (24 * 60),	# time in minutes
    #  threads: 4
    #  benchmark:
    #    "rsem/benchmark/{rna_lib}_index_bam.benchmark"
    #  log:
    #    "rsem/logs/{rna_lib}_index_bam.log"
    #  shell:
    #    """
    #      samtools index -@ {threads} {input.bam} &>> {log}
    #    """
