configfile: "{}/ref.yaml".format(workflow.basedir)

rule pre_proc:
  input:
    expand("runs/{library}/ApplyBQRS/recal.bam", library=config["libraries"])

def ubam_input(wildcards):
  return config["libraries"][wildcards.library]

rule ubam:
  input:
    ubam_input
  output:
    "runs/{library}/ubam/unprocessed.bam"
  conda:
    "envs_dir/pre_proc.yaml"
  shell:
    "python3 {workflow.basedir}/scripts_dir/ubam.py {input} {wildcards.library} {output}"

rule SamToFastqAndBwaMem:
  input:
    rules.ubam.output
  output:
    "runs/{library}/SamToFastqAndBwaMem/bwa.bam"
  conda:
    "envs_dir/pre_proc.yaml"
  log:
    bwa = "runs/{library}/SamToFastqAndBwaMem/bwa.log",
    picard = "runs/{library}/SamToFastqAndBwaMem/picard.log",
    sam = "runs/{library}/SamToFastqAndBwaMem/sam.log",
  params:
    ref_fasta = config["ref_fasta"],
    java_opts = config["SamToFastqAndBwaMem.java_opt"]
  threads: 16
  shell:
    "picard SamToFastq INPUT={input} FASTQ=/dev/stdout INTERLEAVE=true NON_PF=true 2> {log.picard} | bwa mem -K 100000000 -p -v 3 -t {threads} -Y {params.ref_fasta} /dev/stdin - 2> {log.bwa} | samtools view -1 - > {output} 2> {log.sam}"


rule MergeBamAlignment:
  input:
    ubam = rules.ubam.output,
    bwa_bam = rules.SamToFastqAndBwaMem.output,
  output:
    "runs/{library}/MergeBamAlignment/merge.bam"
  conda:
    "envs_dir/pre_proc.yaml"
  log:
    "runs/{library}/MergeBamAlignment/merge.log"
  params:
    ref_fasta = config["ref_fasta"],
    java_opts = config["MergeBamAlignment.java_opt"]
  shell:
    """
      gatk --java-options "{params.java_opts}" MergeBamAlignment \
      --VALIDATION_STRINGENCY SILENT \
      --EXPECTED_ORIENTATIONS FR \
      --ATTRIBUTES_TO_RETAIN X0 \
      --ALIGNED_BAM {input.bwa_bam} \
      --UNMAPPED_BAM {input.ubam} \
      --OUTPUT {output} \
      --REFERENCE_SEQUENCE {params.ref_fasta} \
      --PAIRED_RUN true \
      --SORT_ORDER "unsorted" \
      --IS_BISULFITE_SEQUENCE false \
      --ALIGNED_READS_ONLY false \
      --CLIP_ADAPTERS false \
      --MAX_RECORDS_IN_RAM 2000000 \
      --ADD_MATE_CIGAR true \
      --MAX_INSERTIONS_OR_DELETIONS -1 \
      --PRIMARY_ALIGNMENT_STRATEGY MostDistant \
      --UNMAPPED_READ_STRATEGY COPY_TO_TAG \
      --ALIGNER_PROPER_PAIR_FLAGS true \
			--UNMAP_CONTAMINANT_READS true \
      &> {log}
		"""
# program group deleted

rule MarkDuplicates:
  input:
    rules.MergeBamAlignment.output,
  output:
    bam = "runs/{library}/MarkDuplicates/dup.bam",
    metrics = "runs/{library}/MarkDuplicates/metrics.duplicate_metrics",
  conda:
    "envs_dir/pre_proc.yaml"
  log:
    "runs/{library}/MarkDuplicates/dup.log"
  params:
    ref_fasta = config["ref_fasta"],
    java_opts = config["MarkDuplicates.java_opt"]
  shell:
    """
      gatk --java-options "{params.java_opts}"	MarkDuplicates \
      --INPUT {input} \
      --OUTPUT {output.bam} \
      --METRICS_FILE {output.metrics} \
      --VALIDATION_STRINGENCY SILENT \
      --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
      --ASSUME_SORT_ORDER "queryname" \
      --CREATE_MD5_FILE true	\
      &> {log}
    """
rule SortAndFixTags:
  input:
    rules.MarkDuplicates.output.bam,
  output:
    out_bam = "runs/{library}/SortAndFixTags/sorted.bam",
    out_bai = "runs/{library}/SortAndFixTags/sorted.bai",
  conda:
    "envs_dir/pre_proc.yaml"
  log:
    log1 = "runs/{library}/SortAndFixTags/sorted_sam.log",
    log2 = "runs/{library}/SortAndFixTags/sorted_tags.log",
  params:
    ref_fasta = config["ref_fasta"],
    java_opts_sort = config["SortAndFixTags.java_opt_sort"],
    java_opts_fix = config["SortAndFixTags.java_opt_fix"]
  shell:
    """
      gatk --java-options "{params.java_opts_sort}" SortSam \
      --INPUT {input} \
      --OUTPUT /dev/stdout \
      --SORT_ORDER "coordinate" \
      --CREATE_INDEX false \
      --CREATE_MD5_FILE false \
      2> {log.log1} \
      | \
      gatk --java-options "{params.java_opts_fix}" SetNmAndUqTags \
      --INPUT /dev/stdin \
      --OUTPUT {output.out_bam} \
      --CREATE_INDEX true \
      --CREATE_MD5_FILE true \
			--REFERENCE_SEQUENCE {params.ref_fasta} \
      &> {log.log2}
    """

rule BaseRecalibrator:
  input:
    bam = rules.SortAndFixTags.output.out_bam,
    index = rules.SortAndFixTags.output.out_bai,
  output:
    recal = "runs/{library}/BaseRecalibrator/recal.csv",
  conda:
    "envs_dir/pre_proc.yaml"
  log:
    "runs/{library}/BaseRecalibrator/recal.log",
  params:
    dbSNP_vcf = config["dbSNP_vcf"],
    dbSNP_vcf_index = config["dbSNP_vcf_index"],
    ref_fasta = config["ref_fasta"],
    a1000G = config["a1000G"],
    known_indels = config["known_indels"],
    a1000G_index = config["a1000G_index"],
    known_indels_index = config["known_indels_index"],
    interval = config["exom_padded"],
    java_opts = config["BaseRecalibrator.java_opt"]
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

rule ApplyBQRS:
  input:
    bam = rules.SortAndFixTags.output.out_bam,
    bam_index = rules.SortAndFixTags.output.out_bai,
    recal = rules.BaseRecalibrator.output.recal
  output:
    bam_out = "runs/{library}/ApplyBQRS/recal.bam",
    bai_out = "runs/{library}/ApplyBQRS/recal.bai",
  conda:
    "envs_dir/pre_proc.yaml"
  log:
    "runs/{library}/ApplyBQRS/recal.log",
  params:
    dbSNP_vcf = config["dbSNP_vcf"],
    dbSNP_vcf_index = config["dbSNP_vcf_index"],
    ref_fasta = config["ref_fasta"],
    a1000G = config["a1000G"],
    known_indels = config["known_indels"],
    a1000G_index = config["a1000G_index"],
    known_indels_index = config["known_indels_index"],
    java_opts = config["ApplyBQSR.java_opt"]
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
    """
