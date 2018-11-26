configfile: "{}/ref.yaml".format(workflow.basedir)

def ubam_input(wildcards):
  return config["RNA_libs"][wildcards.library]

rule all:
  input:
    expand("runs/{library}/StarAlign/{library}.Aligned.sortedByCoord.out.bam", library = config["RNA_libs"])

rule ubam:
  input:
    ubam_input
  output:
    "runs/{library}/ubam/unprocessed.bam"
  conda:
    "envs_dir/pre_proc.yaml"
  resources:
    mem_mb = 5000
  benchmark:
    "benchmarks/{library}.ubam.benchmark.txt"
  shell:
    "python3 {workflow.basedir}/scripts_dir/ubam.py {input} {wildcards.library} {output}"

rule SamToFastq:
  input:
    rules.ubam.output
  output:
    fastq1 = "runs/{library}/SamToFastq/fastq_1.gz",
    fastq2 = "runs/{library}/SamToFastq/fastq_2.gz",
  conda:
    "envs_dir/pre_proc.yaml"
  resources:
    mem_mb = 5000
  log:
    picard = "runs/{library}/SamToFastq/picard.log",
  params:
    ref_fasta = config["ref_fasta_rna"],
    java_opts = config["SamToFastqAndBwaMem.java_opt"]
  resources:
    mem_mb = int(config["SamToFastqAndBwaMem.java_opt"].strip("-Xmx")) + 1000
  benchmark:
    "benchmarks/{library}.SamToFastq.benchmark.txt"
  shell:
    """
    picard {params.java_opts} SamToFastq \
      INPUT={input} \
      VALIDATION_STRINGENCY=SILENT \
      FASTQ={output.fastq1} \
      SECOND_END_FASTQ={output.fastq2} 2> {log.picard}     
    """

rule StarGenerateReferences:
  output:
    "/home/dylan/workflow_script/star_genome_index/star-HUMAN-refs.tar.gz"
  conda:
    "envs_dir/pre_proc.yaml"
  params:
    ref_fasta = config["ref_fasta_rna"],
    ref_gtf = config["ref_gtf_rna"],
    star_default_overhang = 100
  threads: 
    8
  resources:
    mem_mb = 100 * 1024 # 100G
  benchmark:
    "benchmarks/stargeneratereferences.benchmark.txt"
  shell:
    """
      if [[ ! (-d STAR2_5) ]]; then
        mkdir STAR2_5
      fi
      
      STAR \
      --runMode genomeGenerate \
      --genomeDir STAR2_5 \
      --genomeFastaFiles {params.ref_fasta} \
      --sjdbGTFfile {params.ref_gtf} \
      --sjdbOverhang {params.star_default_overhang} \
      --runThreadN {threads}
      
      ls STAR2_5
      
      tar -zcvf star-HUMAN-refs.tar.gz STAR2_5
      if [[ ! (-d /home/dylan/workflow_script/star_genome_index) ]]; then
        mkdir /home/dylan/workflow_script/star_genome_index -p
      fi
      mv star-HUMAN-refs.tar.gz /home/dylan/workflow_script/star_genome_index/
		"""

rule StarAlign:
  input:
    rna_index = rules.StarGenerateReferences.output,
    fastq1 = rules.SamToFastq.output.fastq1,
    fastq2 = rules.SamToFastq.output.fastq1,
  output:
    out_bam = "runs/{library}/StarAlign/{library}.Aligned.sortedByCoord.out.bam"
  conda:
    "envs_dir/pre_proc.yaml"
  params:
    ref_fasta = config["ref_fasta_rna"],
    ref_gtf = config["ref_gtf_rna"],
    star_default_overhang = 100,
    rna_ref_loc = "/home/dylan/workflow_script/STAR2_5"
  threads: 8
  resources:
    mem_mb = 45 * 1024 # 45G
  benchmark:
    "benchmarks/{library}.StarAlign.benchmark.txt"
  shell:
    """
			STAR \
			--genomeDir {params.rna_ref_loc} \
			--runThreadN {threads} \
			--readFilesIn {input.fastq1} {input.fastq2} \
			--readFilesCommand "gunzip -c" \
			--sjdbOverhang {params.star_default_overhang} \
			--outSAMtype BAM SortedByCoordinate \
			--twopassMode Basic \
			--limitBAMsortRAM {resources.mem_mb}000000 \
			--limitOutSJcollapsed 1000000 \
      --outFileNamePrefix runs/{wildcards.library}/StarAlign/{wildcards.library}.
    """
# --limitOutSJcollapsed 1000000 is default value

#rule MergeBamAlignment:
#  input:
#    ubam = rules.ubam.output,
#    bwa_bam = rules.SamToFastqAndBwaMem.output,
#  output:
#    "runs/{library}/MergeBamAlignment/merge.bam"
#  conda:
#    "envs_dir/pre_proc.yaml"
#  resources:
#    mem_mb = 5000
#  log:
#    "runs/{library}/MergeBamAlignment/merge.log"
#  params:
#    ref_fasta = config["ref_fasta"],
#    java_opts = config["MergeBamAlignment.java_opt"]
#  resources:
#    mem_mb = int(config["MergeBamAlignment.java_opt"].strip("-Xmx")) + 1000
#  benchmark:
#    "benchmarks/{library}.MergeBamAlignment.benchmark.txt"
#  shell:
#    """
#      gatk --java-options "{params.java_opts}" MergeBamAlignment \
#      --VALIDATION_STRINGENCY SILENT \
#      --EXPECTED_ORIENTATIONS FR \
#      --ATTRIBUTES_TO_RETAIN X0 \
#      --ALIGNED_BAM {input.bwa_bam} \
#      --UNMAPPED_BAM {input.ubam} \
#      --OUTPUT {output} \
#      --REFERENCE_SEQUENCE {params.ref_fasta} \
#      --PAIRED_RUN true \
#      --SORT_ORDER "unsorted" \
#      --IS_BISULFITE_SEQUENCE false \
#      --ALIGNED_READS_ONLY false \
#      --CLIP_ADAPTERS false \
#      --MAX_RECORDS_IN_RAM 2000000 \
#      --ADD_MATE_CIGAR true \
#      --MAX_INSERTIONS_OR_DELETIONS -1 \
#      --PRIMARY_ALIGNMENT_STRATEGY MostDistant \
#      --UNMAPPED_READ_STRATEGY COPY_TO_TAG \
#      --ALIGNER_PROPER_PAIR_FLAGS true \
#			--UNMAP_CONTAMINANT_READS true \
#      &> {log}
#		"""
