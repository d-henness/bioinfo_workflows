
configfile: "{}/ref.yaml".format(workflow.basedir)

include: "fastp_rna.snakefile"
rule rna_variants_all:
  input:
#    expand("GATK_runs/{library}/MergeBamAlignment/merge.bam", library = config["libraries"])
    expand("rna_GATK_runs/{rna_lib}/MergeBamAlignment/{rna_lib}_merge.bam", rna_lib = config["rna_libraries"]),

def get_sample_name(wildcards):
  for sample_name in config["merge_libs"]:
    for lib in config["merge_libs"][sample_name]:
      if lib == config["rna_libraries"][wildcards.rna_lib]:
        print(sample_name)

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
  threads: 8
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
    ref_fasta = config["ref_fasta"],
    java_opts = config["MergeBamAlignment.java_opt"],
  resources:
    mem_mb = lambda wildcards, attempt: attempt * (int(config["MergeBamAlignment.java_opt"].strip("-Xmx")) + 1000),
    time_min = lambda wildcards, attempt: attempt * 24 * 60,	# time in minutes
  benchmark:
    "rna_GATK_runs/benchmarks/{rna_lib}_MergeBamAlignment.benchmark.txt"
  shell:
    """
      gatk --java-options "{params.java_opts}" MergeBamAlignment \
      --REFERENCE_SEQUENCE /data/shared/hg38/ftp.ensembl.org/pub/release-95/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
      --UNMAPPED_BAM {input.ubam} \
      --ALIGNED_BAM {input.aligned_bam} \
      --OUTPUT {output} \
      --INCLUDE_SECONDARY_ALIGNMENTS false \
      --PAIRED_RUN false \
      --VALIDATION_STRINGENCY SILENT \
      --MAX_RECORDS_IN_RAM 2000000 \
      &> {log}
		"""
