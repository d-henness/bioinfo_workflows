configfile: "{}/ref.yaml".format(workflow.basedir)

rule make_alig:
  input:
#    expand("GATK_runs/{library}/MergeBamAlignment/merge.bam", library = config["libraries"])
    expand("GATK_runs/{library}/ubam/unprocessed.bam", library = config["libraries"])

def get_sample_name(wildcards):
  for sample_name in config["merge_libs"]:
    for lib in config["merge_libs"][sample_name]:
      if lib == config["libraries"][wildcards.library]:
        print(sample_name)

def ubam_input(wildcards):
  return config["libraries"][wildcards.library]

rule ubam:
  input:
    ubam_input
  output:
    temp("GATK_runs/{library}/ubam/unprocessed.bam")
  conda:
    "envs_dir/pre_proc.yaml"
  resources:
    mem_mb = lambda wildcards, attempt: attempt * 5000,
    time_min = lambda wildcards, attempt: attempt * 24 * 60,	# time in minutes
  params:
    exclude_list = '',
    sample_name = get_sample_name
  log:
    ubam = "GATK_runs/{library}/ubam/ubam.log"
  benchmark:
    "benchmarks/{library}.ubam.benchmark.txt"
  shell:
    "python3 {workflow.basedir}/scripts_dir/ubam.py {input} {wildcards.library} {output} 2> {log.ubam} {params.sample_name}"

rule SamToFastqAndBwaMem:
  input:
    rules.ubam.output
  output:
    temp("GATK_runs/{library}/SamToFastqAndBwaMem/bwa.bam")
  conda:
    "envs_dir/pre_proc.yaml"
  log:
    bwa = "GATK_runs/{library}/SamToFastqAndBwaMem/bwa.log",
    picard = "GATK_runs/{library}/SamToFastqAndBwaMem/picard.log",
    sam = "GATK_runs/{library}/SamToFastqAndBwaMem/sam.log",
  params:
    ref_fasta = config["ref_fasta"],
    java_opts = config["SamToFastqAndBwaMem.java_opt"],
    exclude_list = ''
  threads: 8
  resources:
    mem_mb = lambda wildcards, attempt: attempt * int(config["SamToFastqAndBwaMem.java_opt"].strip("-Xmx")) + 1000,
    time_min = lambda wildcards, attempt: attempt * 24 * 60,	# time in minutes
  benchmark:
    "benchmarks/{library}.SamToFastqAndBwaMem.benchmark.txt"
  shell:
    "picard {params.java_opts} SamToFastq INPUT={input} FASTQ=/dev/stdout INTERLEAVE=true NON_PF=true 2> {log.picard} | bwa mem -K 100000000 -p -v 3 -t {threads} -Y {params.ref_fasta} /dev/stdin - 2> {log.bwa} | samtools view -1 - > {output} 2> {log.sam}"


rule MergeBamAlignment:
  input:
    ubam = rules.ubam.output,
    bwa_bam = rules.SamToFastqAndBwaMem.output,
  output:
    temp("GATK_runs/{library}/MergeBamAlignment/merge.bam")
  conda:
    "envs_dir/pre_proc.yaml"
  log:
    "GATK_runs/{library}/MergeBamAlignment/merge.log"
  params:
    ref_fasta = config["ref_fasta"],
    java_opts = config["MergeBamAlignment.java_opt"],
    exclude_list = ''
  resources:
    mem_mb = lambda wildcards, attempt: attempt * (int(config["MergeBamAlignment.java_opt"].strip("-Xmx")) + 1000),
    time_min = lambda wildcards, attempt: attempt * 24 * 60,	# time in minutes
  benchmark:
    "benchmarks/{library}.MergeBamAlignment.benchmark.txt"
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
