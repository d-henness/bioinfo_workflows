
configfile: "{}/ref.yaml".format(workflow.basedir)

include: "fastp_rna.snakefile"
rule rna_variants_all:
  input:
#    expand("GATK_runs/{library}/MergeBamAlignment/merge.bam", library = config["libraries"])
    expand("rna_GATK_runs/{rna_lib}/ApplyBQSR/{rna_lib}_recal.bam", rna_lib = config["rna_libraries"]),

def get_sample_name(wildcards):
  for sample_name in config["rna_merge_libs"]:
    for lib in config["rna_merge_libs"][sample_name]:
      if lib == config["rna_libraries"][wildcards.rna_lib]:
        print(sample_name)

#def MarkDuplicates_input(wildcards):
#  bam_files = []
#  for library in config["merge_libs"][wildcards.merge]:
#    bam_files.append(f"rna_GATK_runs/{library}/MergeBamAlignment/merge.bam")
#  return bam_files
#
#def MarkDuplicates_input_string(wildcards):
#  bam_files = []
#  for library in config["merge_libs"][wildcards.merge]:
#    bam_files.append(f"--INPUT /tmp/$SLURM_JOB_ID/{wildcards.merge}/MarkDuplicates/{library}.bam")
#  return ' '.join(bam_files)


rule ubam:
  input:
    fq1 = rules.fastp_paired.output.fq1_out,
    fq2 = rules.fastp_paired.output.fq2_out,
  output:
    "rna_GATK_runs/{rna_lib}/ubam/{rna_lib}_unprocessed.bam",
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
    ubam = "rna_GATK_runs/log/{rna_lib}_ubam.log"
  benchmark:
    "rna_GATK_runs/benchmarks/{rna_lib}_ubam.benchmark.txt"
  shell:
    "python3 {workflow.basedir}/scripts_dir/ubam.py {input} {wildcards.rna_lib} {output} 2> {log.ubam} {params.sample_name}"

rule star:
  input:
    fq1 = rules.fastp_paired.output.fq1_out,
    fq2 = rules.fastp_paired.output.fq2_out,
  output:
    "rna_GATK_runs/{rna_lib}/star/{rna_lib}.Aligned.sortedByCoord.out.bam",
  conda:
    "envs_dir/rna_variants.yaml"
  threads: 16
  resources:
    mem_mb = lambda wildcards, attempt: attempt * 40000,
    time_min = lambda wildcards, attempt: attempt * 24 * 60,	# time in minutes
    io = 1, # used to indicate that this is an io heavy job and should not have many running at once
  params:
    star_genomeDir = config["star_genomeDir"],
  log:
    "rna_GATK_runs/log/{rna_lib}_star.log"
  benchmark:
    "rna_GATK_runs/benchmarks/{rna_lib}_star.benchmark.txt"
  shell:
    """
      if [[ !(-d rna_GATK_runs/{wildcards.rna_lib}/star) ]]; then
        mkdir -p rna_GATK_runs/{wildcards.rna_lib}/star
      fi

      STAR \
        --runThreadN {threads} \
        --genomeDir {params.star_genomeDir} \
        --readFilesIn {input.fq1} {input.fq2} \
        --readFilesCommand "gunzip -c" \
        --outSAMtype BAM SortedByCoordinate \
        --outTmpDir /tmp/$SLURM_JOB_ID \
        --twopassMode Basic \
        --limitBAMsortRAM 8000000000 \
        --outFileNamePrefix rna_GATK_runs/{wildcards.rna_lib}/star/{wildcards.rna_lib}. &> {log}
    """

rule rna_MergeBamAlignment:
  input:
    ubam = rules.ubam.output,
    aligned_bam = rules.star.output,
  output:
    "rna_GATK_runs/{rna_lib}/MergeBamAlignment/{rna_lib}_merge.bam"
  conda:
    "envs_dir/pre_proc.yaml"
  log:
    "rna_GATK_runs/log/{rna_lib}_MergeBamAlignment.log"
  params:
    ref_fasta = config["ensembl_95_fasta"],
    java_opts = config["MergeBamAlignment.java_opt"],
  resources:
    mem_mb = lambda wildcards, attempt: attempt * (int(config["MergeBamAlignment.java_opt"].strip("-Xmx")) + 1000),
    time_min = lambda wildcards, attempt: attempt * 24 * 60,	# time in minutes
  benchmark:
    "rna_GATK_runs/benchmarks/{rna_lib}_MergeBamAlignment.benchmark.txt"
  shell:
    """
      gatk --java-options "{params.java_opts}" MergeBamAlignment \
      --REFERENCE_SEQUENCE {params.ref_fasta} \
      --UNMAPPED_BAM {input.ubam} \
      --ALIGNED_BAM {input.aligned_bam} \
      --OUTPUT {output} \
      --INCLUDE_SECONDARY_ALIGNMENTS false \
      --PAIRED_RUN false \
      --VALIDATION_STRINGENCY SILENT \
      --MAX_RECORDS_IN_RAM 2000000 \
      &> {log}
		"""

rule rna_MarkDuplicates:
  input:
    rules.rna_MergeBamAlignment.output
  output:
    bam = "rna_GATK_runs/{rna_lib}/MarkDuplicates/{rna_lib}_dup.bam",
    metrics = "rna_GATK_runs/{rna_lib}/MarkDuplicates/{rna_lib}_metrics.duplicate_metrics",
  conda:
    "envs_dir/pre_proc.yaml"
  log:
    "rna_GATK_runs/log/{rna_lib}_MarkDuplicates.log"
  params:
    ref_fasta = config["ensembl_95_fasta"],
    java_opts = config["MarkDuplicates.java_opt"],
    exclude_list = '',
    temp_dir = "/tmp/$SLURM_JOB_ID/{rna_lib}/MarkDuplicates/tmp",
  resources:
    mem_mb = lambda wildcards, attempt: attempt * (int(config["MarkDuplicates.java_opt"].strip("-Xmx")) + 1000),
    time_min = lambda wildcards, attempt: attempt * 24 * 60,	# time in minutes
    io = 1, # used to indicate that this is an io heavy job and should not have many running at once
  benchmark:
    "rna_GATK_runs/benchmarks/{rna_lib}_MarkDuplicates.benchmark.txt"
  shell:
    """
      for file in {input}; do
        lib=$(echo $file | sed -E 's,rna_GATK_runs/([^/]+)/.*,\\1,')
        mkdir -p {params.temp_dir}
        cp $file /tmp/$SLURM_JOB_ID/{wildcards.rna_lib}/MarkDuplicates/$lib.bam
      done

      gatk --java-options "{params.java_opts}"	MarkDuplicates \
      --INPUT=/tmp/$SLURM_JOB_ID/{wildcards.rna_lib}/MarkDuplicates/$lib.bam \
      --OUTPUT {output.bam} \
      --METRICS_FILE {output.metrics} \
      --VALIDATION_STRINGENCY SILENT \
      --TAGGING_POLICY All \
      --TMP_DIR {params.temp_dir} \
      &> {log}
    """

rule setup_GATK3:
  output:
    signal = "GATK3_setup/signal.txt"
  conda:
    "envs_dir/pre_proc_GATK3.yaml"
  log:
    "GATK3_setup/setup.log"
  params:
    GATK3_bz2 = config["GATK3_bz2"],
  resources:
    mem_mb = lambda wildcards, attempt: attempt * (int(config["MarkDuplicates.java_opt"].strip("-Xmx")) + 1000),
    time_min = lambda wildcards, attempt: attempt * 24 * 60,	# time in minutes
    io = 1, # used to indicate that this is an io heavy job and should not have many running at once
  shell:
    """
      gatk3-register {params.GATK3_bz2} &> {log}
      if [[ !(-d GATK3_setup) ]]; then
        mkdir GATK3_setup
      fi
      echo "setup GATK3 successfully" > {output.signal}
    """


rule rna_SplitNCigarReads:
  input:
    signal = rules.setup_GATK3.output.signal,
    bam = rules.rna_MarkDuplicates.output.bam
  output:
    bam = "rna_GATK_runs/{rna_lib}/SplitNCigarReads/{rna_lib}.bam",
    index = "rna_GATK_runs/{rna_lib}/SplitNCigarReads/{rna_lib}.bia",
  conda:
    "envs_dir/pre_proc_GATK3.yaml"
  log:
    "rna_GATK_runs/log/{rna_lib}_SplitNCigarReads.log"
  params:
    ref_fasta = config["ensembl_95_fasta"],
    exclude_list = '',
  resources:
    mem_mb = lambda wildcards, attempt: attempt * (int(config["MarkDuplicates.java_opt"].strip("-Xmx")) + 1000),
    time_min = lambda wildcards, attempt: attempt * 24 * 60,	# time in minutes
    io = 1, # used to indicate that this is an io heavy job and should not have many running at once
  benchmark:
    "rna_GATK_runs/benchmarks/{rna_lib}_SplitNCigarReads.benchmark.txt"
  shell:
    """
      samtools index {input.bam}
      java -Xmx{resources.mem_mb}m -jar $CONDA_PREFIX/opt/gatk-3.8/GenomeAnalysisTK.jar \
        -T SplitNCigarReads \
        -R {params.ref_fasta} \
        -I {input.bam} \
        -o {output.bam} \
        -rf ReassignOneMappingQuality \
        -RMQF 255 \
    		-RMQT 60 \
        -U ALLOW_N_CIGAR_READS &> {log}
    """

rule rna_BaseRecalibrator:
  input:
    bam = rules.rna_SplitNCigarReads.output.bam,
    index = rules.rna_SplitNCigarReads.output.index,
  output:
    recal = "rna_GATK_runs/{rna_lib}/BaseRecalibrator/{rna_lib}_recal.csv",
#   temp_dir = temp (directory("GATK_runs/{merge}/BaseRecalibrator/tmp/")),
  conda:
    "envs_dir/pre_proc.yaml"
  log:
    "rna_GATK_runs/log/{rna_lib}_BaseRecalibrator.log"
  params:
    dbSNP_vcf = config["dbSNP_vcf"],
    dbSNP_vcf_index = config["dbSNP_vcf_index"],
    ref_fasta = config["ref_fasta"],
    a1000G = config["a1000G"],
    known_indels = config["known_indels"],
    a1000G_index = config["a1000G_index"],
    known_indels_index = config["known_indels_index"],
    interval = config["alt_bed"],
    java_opts = config["BaseRecalibrator.java_opt"],
    exclude_list = ''
  resources:
    mem_mb = lambda wildcards, attempt: attempt * (int(config["BaseRecalibrator.java_opt"].strip("-Xmx")) + 1000),
    time_min = lambda wildcards, attempt: attempt * 24 * 60,	# time in minutes
  benchmark:
    "rna_GATK_runs/benchmarks/{rna_lib}_BaseRecalibrator.benchmark.txt"
  shell:
    """
      gatk --java-options "{params.java_opts} -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:+PrintFlagsFinal \
            -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps -XX:+PrintGCDetails \
            -Xloggc:gc_log.log" BaseRecalibrator \
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
    bam = rules.rna_SplitNCigarReads.output.bam,
    index = rules.rna_SplitNCigarReads.output.index,
    recal = rules.rna_BaseRecalibrator.output.recal
  output:
    bam_out = "rna_GATK_runs/{rna_lib}/ApplyBQSR/{rna_lib}_recal.bam",
    bai_out = "rna_GATK_runs/{rna_lib}/ApplyBQSR/{rna_lib}_recal.bam.bai",
    stats = "rna_GATK_runs/{rna_lib}/ApplyBQSR/{rna_lib}_recal.stats",
  conda:
    "envs_dir/pre_proc.yaml"
  log:
    "rna_GATK_runs/log/{rna_lib}_ApplyBQSR.log"
  params:
    dbSNP_vcf = config["dbSNP_vcf"],
    dbSNP_vcf_index = config["dbSNP_vcf_index"],
    ref_fasta = config["ref_fasta"],
    a1000G = config["a1000G"],
    known_indels = config["known_indels"],
    a1000G_index = config["a1000G_index"],
    known_indels_index = config["known_indels_index"],
    java_opts = config["ApplyBQSR.java_opt"],
    interval = config["alt_bed"],
    exclude_list = ''
  resources:
    mem_mb = lambda wildcards, attempt: attempt * (int(config["ApplyBQSR.java_opt"].strip("-Xmx")) + 1000),
    time_min = lambda wildcards, attempt: attempt * 24 * 60,	# time in minutes
  benchmark:
    "rna_GATK_runs/benchmarks/{rna_lib}_ApplyBQSR.benchmark.txt"
  shell:
    """
      gatk --java-options "-XX:+PrintFlagsFinal -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps \
            -XX:+PrintGCDetails -Xloggc:gc_log.log \
            -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 {params.java_opts}" ApplyBQSR \
        -R {params.ref_fasta} \
        -I {input.bam} \
        -O {output.bam_out} \
        -bqsr {input.recal} \
        --add-output-sam-program-record \
        --use-original-qualities \
        &> {log}
      if [[ -e GATK_runs/{wildcards.rna_lib}/ApplyBQSR/{wildcards.rna_lib}_recal.bai ]]; then
        mv GATK_runs/{wildcards.rna_lib}/ApplyBQSR/{wildcards.rna_lib}_recal.bai GATK_runs/{wildcards.rna_lib}/ApplyBQSR/{wildcards.rna_lib}_recal.bam.bai
      fi
      samtools stats -t {params.interval} {output.bam_out} > {output.stats}
    """
