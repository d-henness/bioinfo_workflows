configfile: "{}/ref.yaml".format(workflow.basedir)

include: "pre_pro_bf_merge_alt_bed.snakefile"

merged = config["merge_libs"]

rule pre_proc:
  input:
    expand("GATK_runs/{merge}/ApplyBQSR/{merge}_recal.bam", merge = merged)

def MarkDuplicates_input(wildcards):
  bam_files = []
  for library in config["merge_libs"][wildcards.merge]:
    bam_files.append("".join(["GATK_runs/",library,"/MergeBamAlignment/merge.bam"]))
# f strings stopped working here for some reason
#    bam_files.append(f"GATK_runs/{test}/MergeBamAlignment/merge.bam")
  return bam_files

def MarkDuplicates_input_string(wildcards):
  bam_files = []
  for library in config["merge_libs"][wildcards.merge]:
    bam_files.append("".join(["--INPUT GATK_runs/",library,"/MergeBamAlignment/merge.bam"]))
# f strings stopped working here for some reason
#    bam_files.append(f"--INPUT GATK_runs/{library}/MergeBamAlignment/merge.bam")
  return ' '.join(bam_files)

rule MarkDuplicates:
  input:
    MarkDuplicates_input
  output:
    bam = temp("GATK_runs/{merge}/MarkDuplicates/dup.bam"),
    metrics = "GATK_runs/{merge}/MarkDuplicates/metrics.duplicate_metrics",
  conda:
    "envs_dir/pre_proc.yaml"
  log:
    "GATK_runs/{merge}/MarkDuplicates/dup.log"
  params:
    ref_fasta = config["ref_fasta"],
    java_opts = config["MarkDuplicates.java_opt"],
    libs_string = MarkDuplicates_input_string,
    exclude_list = '',
    temp_dir = "/tmp/$SLURM_JOB_ID/{merge}/MarkDuplicates/tmp",
  resources:
    mem_mb = lambda wildcards, attempt: attempt * (int(config["MarkDuplicates.java_opt"].strip("-Xmx")) + 1000),
    time_min = lambda wildcards, attempt: attempt * 24 * 60,	# time in minutes
    io = 1, # used to indicate that this is an io heavy job and should not have many running at once
  benchmark:
    "benchmarks/{merge}.MarkDuplicates.benchmark.txt"
  shell:
    """
      mkdir -p {params.temp_dir}

      gatk --java-options "{params.java_opts}"	MarkDuplicates \
      {params.libs_string} \
      --OUTPUT {output.bam} \
      --METRICS_FILE {output.metrics} \
      --VALIDATION_STRINGENCY SILENT \
      --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
      --ASSUME_SORT_ORDER "queryname" \
      --TAGGING_POLICY All \
      --CREATE_MD5_FILE true	\
      --TMP_DIR {params.temp_dir} \
      &> {log}
    """
rule SortAndFixTags:
  input:
    rules.MarkDuplicates.output.bam,
  output:
    out_bam = temp("GATK_runs/{merge}/SortAndFixTags/sorted.bam"),
    out_bai = temp("GATK_runs/{merge}/SortAndFixTags/sorted.bai"),
  conda:
    "envs_dir/pre_proc.yaml"
  log:
    log1 = "GATK_runs/{merge}/SortAndFixTags/sorted_sam.log",
    log2 = "GATK_runs/{merge}/SortAndFixTags/sorted_tags.log",
  params:
    ref_fasta = config["ref_fasta"],
    java_opts_sort = config["SortAndFixTags.java_opt_sort"],
    java_opts_fix = config["SortAndFixTags.java_opt_fix"],
    exclude_list = '',
    temp_dir = "/tmp/$SLURM_JOB_ID/{merge}/SortAndFixTags/tmp",
  resources:
    mem_mb = lambda wildcards, attempt: attempt * (int(config["SortAndFixTags.java_opt_sort"].strip("-Xmx")) + int(config["SortAndFixTags.java_opt_fix"].strip("-Xmx"))),
    time_min = lambda wildcards, attempt: attempt * 24 * 60,	# time in minutes
  benchmark:
    "benchmarks/{merge}.SortAndFixTags.benchmark.txt"
  shell:
    """
      mkdir -p {params.temp_dir}

      gatk --java-options "{params.java_opts_sort}" SortSam \
      --INPUT {input} \
      --OUTPUT /dev/stdout \
      --SORT_ORDER "coordinate" \
      --CREATE_INDEX false \
      --CREATE_MD5_FILE false \
      --TMP_DIR {params.temp_dir} \
      2> {log.log1} \
      | \
      gatk --java-options "{params.java_opts_fix}" SetNmAndUqTags \
      --INPUT /dev/stdin \
      --OUTPUT {output.out_bam} \
      --CREATE_INDEX true \
      --CREATE_MD5_FILE true \
			--REFERENCE_SEQUENCE {params.ref_fasta} \
      --TMP_DIR {params.temp_dir} \
      &> {log.log2}
    """

rule BaseRecalibrator:
  input:
    bam = rules.SortAndFixTags.output.out_bam,
    index = rules.SortAndFixTags.output.out_bai,
  output:
    recal = temp("GATK_runs/{merge}/BaseRecalibrator/recal.csv"),
#   temp_dir = temp (directory("GATK_runs/{merge}/BaseRecalibrator/tmp/")),
  conda:
    "envs_dir/pre_proc.yaml"
  log:
    "GATK_runs/{merge}/BaseRecalibrator/recal.log",
  params:
    dbSNP_vcf = config["dbSNP_vcf"],
    dbSNP_vcf_index = config["dbSNP_vcf_index"],
    ref_fasta = config["ref_fasta"],
    a1000G = config["a1000G"],
    known_indels = config["known_indels"],
    a1000G_index = config["a1000G_index"],
    known_indels_index = config["known_indels_index"],
    interval = f"-L {config['alt_bed']}" if config['alt_bed'] is not None else "",
    java_opts = config["BaseRecalibrator.java_opt"],
    exclude_list = '',
    temp_dir = "/tmp/$SLURM_JOB_ID/{merge}/BaseRecalibrator/tmp",
  resources:
    mem_mb = lambda wildcards, attempt: attempt * (int(config["BaseRecalibrator.java_opt"].strip("-Xmx"))),
    time_min = lambda wildcards, attempt: attempt * 24 * 60,	# time in minutes
  benchmark:
    "benchmarks/{merge}.BaseRecalibrator.benchmark.txt"
  shell:
    """
      mkdir -p {params.temp_dir}

      gatk --java-options "{params.java_opts}" BaseRecalibrator \
      -R {params.ref_fasta} \
      -I {input.bam} \
      {params.interval} \
      --use-original-qualities \
      -O {output.recal} \
      --known-sites {params.dbSNP_vcf} \
      --known-sites {params.known_indels} \
      --known-sites {params.a1000G} \
      &> {log}
    """

rule ApplyBQSR:
  input:
    bam = rules.SortAndFixTags.output.out_bam,
    bam_index = rules.SortAndFixTags.output.out_bai,
    recal = rules.BaseRecalibrator.output.recal
  output:
    bam_out = "GATK_runs/{merge}/ApplyBQSR/{merge}_recal.bam",
    bai_out = "GATK_runs/{merge}/ApplyBQSR/{merge}_recal.bam.bai",
    stats = "GATK_runs/{merge}/ApplyBQSR/{merge}_recal_stats.txt",
    stats_wo_bed = "GATK_runs/{merge}/ApplyBQSR/{merge}_recal_stats_wo_bed.txt",
  conda:
    "envs_dir/pre_proc.yaml"
  log:
    "GATK_runs/{merge}/ApplyBQSR/recal.log"
  params:
    dbSNP_vcf = config["dbSNP_vcf"],
    dbSNP_vcf_index = config["dbSNP_vcf_index"],
    ref_fasta = config["ref_fasta"],
    a1000G = config["a1000G"],
    known_indels = config["known_indels"],
    a1000G_index = config["a1000G_index"],
    known_indels_index = config["known_indels_index"],
    interval = f"-t {config['alt_bed']}" if config['alt_bed'] != None else "--reference {config[\"ref_fasta\"]}",
    java_opts = config["ApplyBQSR.java_opt"],
    temp_dir = "/tmp/$SLURM_JOB_ID/{merge}/ApplyBQSR/tmp",
    exclude_list = ''
  resources:
    mem_mb = lambda wildcards, attempt: attempt * (int(config["ApplyBQSR.java_opt"].strip("-Xmx"))),
    time_min = lambda wildcards, attempt: attempt * 24 * 60,	# time in minutes
  benchmark:
    "benchmarks/{merge}.ApplyBQSR.benchmark.txt"
  shell:
    """
      mkdir -p {params.temp_dir}

      gatk --java-options "{params.java_opts}" ApplyBQSR \
        -R {params.ref_fasta} \
        -I {input.bam} \
        -O {output.bam_out} \
        -bqsr {input.recal} \
        --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 \
        --add-output-sam-program-record \
        --create-output-bam-md5 \
        --use-original-qualities \
        &> {log}
      if [[ -e GATK_runs/{wildcards.merge}/ApplyBQSR/{wildcards.merge}_recal.bai ]]; then
        mv GATK_runs/{wildcards.merge}/ApplyBQSR/{wildcards.merge}_recal.bai GATK_runs/{wildcards.merge}/ApplyBQSR/{wildcards.merge}_recal.bam.bai
      fi
      samtools stats {params.interval} {output.bam_out} > {output.stats}
      samtools stats {output.bam_out} > {output.stats_wo_bed}
    """
