configfile: "{}/ref.yaml".format(workflow.basedir)

include: "pre_pro_bf_merge.snakefile"

merged = config["merge_libs"]

rule pre_proc:
  input:
    expand("GATK_runs/{merge}/ApplyBQSR/recal.bam", merge = merged)

def MarkDuplicates_input(wildcards):
  bam_files = []
  for library in config["merge_libs"][wildcards.merge]:
    bam_files.append(f"GATK_runs/{library}/MergeBamAlignment/merge.bam")
  return bam_files

def MarkDuplicates_input_string(wildcards):
  bam_files = []
  for library in config["merge_libs"][wildcards.merge]:
    bam_files.append(f"--INPUT GATK_runs/{library}/MergeBamAlignment/merge.bam")
  return ' '.join(bam_files)

rule MarkDuplicates:
  input:
    MarkDuplicates_input
  output:
    bam = temp("GATK_runs/{merge}/MarkDuplicates/dup.bam"),
    metrics = temp("GATK_runs/{merge}/MarkDuplicates/metrics.duplicate_metrics"),
    temp_dir = temp(directory("GATK_runs/{merge}/MarkDuplicates/tmp/"))
  conda:
    "envs_dir/pre_proc.yaml"
  log:
    "GATK_runs/{merge}/MarkDuplicates/dup.log"
  params:
    ref_fasta = config["ref_fasta"],
    java_opts = config["MarkDuplicates.java_opt"],
    libs_string = MarkDuplicates_input_string,
    exclude_list = ''
  resources:
    mem_mb = lambda wildcards, attempt: attempt * (int(config["MarkDuplicates.java_opt"].strip("-Xmx")) + 1000),
    time_min = lambda wildcards, attempt: attempt * 24 * 60,	# time in minutes
  benchmark:
    "benchmarks/{merge}.MarkDuplicates.benchmark.txt"
  shell:
    """
      gatk --java-options "{params.java_opts}"	MarkDuplicates \
      {params.libs_string} \
      --OUTPUT {output.bam} \
      --METRICS_FILE {output.metrics} \
      --VALIDATION_STRINGENCY SILENT \
      --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
      --ASSUME_SORT_ORDER "queryname" \
      --CREATE_MD5_FILE true	\
      --TMP_DIR {output.temp_dir} \
      &> {log}
    """
rule SortAndFixTags:
  input:
    rules.MarkDuplicates.output.bam,
  output:
    out_bam = temp("GATK_runs/{merge}/SortAndFixTags/sorted.bam"),
    out_bai = temp("GATK_runs/{merge}/SortAndFixTags/sorted.bai"),
    temp_dir = temp(directory("GATK_runs/{merge}/SortAndFixTags/tmp/"))
  conda:
    "envs_dir/pre_proc.yaml"
  log:
    log1 = "GATK_runs/{merge}/SortAndFixTags/sorted_sam.log",
    log2 = "GATK_runs/{merge}/SortAndFixTags/sorted_tags.log",
  params:
    ref_fasta = config["ref_fasta"],
    java_opts_sort = config["SortAndFixTags.java_opt_sort"],
    java_opts_fix = config["SortAndFixTags.java_opt_fix"],
    exclude_list = ''
  resources:
    mem_mb = lambda wildcards, attempt: attempt * (int(config["SortAndFixTags.java_opt_sort"].strip("-Xmx")) + int(config["SortAndFixTags.java_opt_fix"].strip("-Xmx")) + 1000),
    time_min = lambda wildcards, attempt: attempt * 24 * 60,	# time in minutes
  benchmark:
    "benchmarks/{merge}.SortAndFixTags.benchmark.txt"
  shell:
    """
      gatk --java-options "{params.java_opts_sort}" SortSam \
      --INPUT {input} \
      --OUTPUT /dev/stdout \
      --SORT_ORDER "coordinate" \
      --CREATE_INDEX false \
      --CREATE_MD5_FILE false \
      --TMP_DIR {output.temp_dir} \
      2> {log.log1} \
      | \
      gatk --java-options "{params.java_opts_fix}" SetNmAndUqTags \
      --INPUT /dev/stdin \
      --OUTPUT {output.out_bam} \
      --CREATE_INDEX true \
      --CREATE_MD5_FILE true \
			--REFERENCE_SEQUENCE {params.ref_fasta} \
      --TMP_DIR {output.temp_dir} \
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
    interval = config["exom_padded"],
    java_opts = config["BaseRecalibrator.java_opt"],
    exclude_list = ''
  resources:
    mem_mb = lambda wildcards, attempt: attempt * (int(config["BaseRecalibrator.java_opt"].strip("-Xmx")) + 1000),
    time_min = lambda wildcards, attempt: attempt * 24 * 60,	# time in minutes
  benchmark:
    "benchmarks/{merge}.BaseRecalibrator.benchmark.txt"
  shell:
    """
      gatk --java-options "{params.java_opts}" BaseRecalibrator \
      -R {params.ref_fasta} \
      -I {input.bam} \
      -L {params.interval} \
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
    bam_out = "GATK_runs/{merge}/ApplyBQSR/recal.bam",
    bai_out = "GATK_runs/{merge}/ApplyBQSR/recal.bam.bai",
  conda:
    "envs_dir/pre_proc.yaml"
  log:
    "GATK_runs/{merge}/ApplyBQSR/recal.log",
  params:
    dbSNP_vcf = config["dbSNP_vcf"],
    dbSNP_vcf_index = config["dbSNP_vcf_index"],
    ref_fasta = config["ref_fasta"],
    a1000G = config["a1000G"],
    known_indels = config["known_indels"],
    a1000G_index = config["a1000G_index"],
    known_indels_index = config["known_indels_index"],
    java_opts = config["ApplyBQSR.java_opt"],
    exclude_list = ''
  resources:
    mem_mb = lambda wildcards, attempt: attempt * (int(config["ApplyBQSR.java_opt"].strip("-Xmx")) + 1000),
    time_min = lambda wildcards, attempt: attempt * 24 * 60,	# time in minutes
  benchmark:
    "benchmarks/{merge}.ApplyBQSR.benchmark.txt"
  shell:
    """
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
      if [[ -e GATK_runs/{wildcards.merge}/ApplyBQSR/recal.bai ]]; then
        mv GATK_runs/{wildcards.merge}/ApplyBQSR/recal.bai GATK_runs/{wildcards.merge}/ApplyBQSR/recal.bam.bai
      fi
    """
