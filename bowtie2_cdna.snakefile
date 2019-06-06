configfile: "{}/ref.yaml".format(workflow.basedir)

def bowtie2_input(wildcards):
  if len(config['libraries'][wildcards.rna_lib]) == 1:
    return f"-U bowtie2_cdna/{wildcards.rna_lib}/fastp/{wildcards.rna_lib}_1.fq.gz"
  elif len(config['libraries'][wildcards.rna_lib]) == 2:
    return f"-1 bowtie2_cdna/{wildcards.rna_lib}/fastp/{wildcards.rna_lib}_1.fq.gz -2 bowtie2_cdna/{wildcards.rna_lib}/fastp/{wildcards.rna_lib}_2.fq.gz"

def fastp_paired_input(wildcards):
  if len(config['libraries'][wildcards.rna_lib]) == 2:
    return [config['libraries'][wildcards.rna_lib][0], config['libraries'][wildcards.rna_lib][1]]
  else:
    return ""

def fastp_unpaired_input(wildcards):
  if len(config['libraries'][wildcards.rna_lib]) == 1:
    return f"{config['libraries'][wildcards.rna_lib][0]}"
  else:
    return ""

rule all:
  input:
    expand("bowtie2_cdna/{rna_lib}/{rna_lib}.sam", rna_lib = config["libraries"])

rule fastp_paired:
  input:
    fastp_paired_input
  output:
    fq1_out = "bowtie2_cdna/{rna_lib}/fastp/{rna_lib}_1.fq.gz",
    fq2_out = "bowtie2_cdna/{rna_lib}/fastp/{rna_lib}_2.fq.gz",
    signal = "bowtie2_cdna/{rna_lib}/fastp/{rna_lib}_signal.txt",
  conda:
    "envs_dir/bowtie2.yaml"
  resources:
    mem_mb = lambda wildcards, attempt: attempt * 5000,
    time_min = lambda wildcards, attempt: attempt * 24 * 60,  # time in minutes
  threads: 1
  benchmark:
    "bowtie2_cdna/benchmark/{rna_lib}_fastp.log"
  log:
    "bowtie2_cdna/log/{rna_lib}_fastp.log"
  shell:
    """
      fastp -i {input[0]} -I {input[1]} -o {output.fq1_out} -O {output.fq2_out} &> {log}
      echo "finished" > {output.signal}
    """

rule fastp_unpaired:
  input:
    fastp_unpaired_input
  output:
    fq1_out = "bowtie2_cdna/{rna_lib}/fastp/{rna_lib}_1.fq.gz",
    signal = "bowtie2_cdna/{rna_lib}/fastp/{rna_lib}_signal.txt",
  conda:
    "envs_dir/bowtie2.yaml"
  resources:
    mem_mb = lambda wildcards, attempt: attempt * 5000,
    time_min = lambda wildcards, attempt: attempt * 24 * 60,  # time in minutes
  threads: 1
  benchmark:
    "bowtie2_cdna/benchmark/{rna_lib}_fastp_up.log"
  log:
    "bowtie2_cdna/log/{rna_lib}_fastp_up.log"
  shell:
    """
      fastp -i {input} -o {output.fq1_out} --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA &> {log}
      echo "finished" > {output.signal}
    """

rule bowtie2:
  input:
    signal = "bowtie2_cdna/{rna_lib}/fastp/{rna_lib}_signal.txt",
  output:
    sam_file = "bowtie2_cdna/{rna_lib}/{rna_lib}.sam"
  conda:
    "envs_dir/bowtie2.yaml"
  params:
    regions = config["bowtie2_cdna"],
    input_cmd = bowtie2_input
  resources:
    mem_mb = lambda wildcards, attempt: attempt * 5000,
    time_min = lambda wildcards, attempt: attempt * 24 * 60,  # time in minutes
  threads: 1
  benchmark:
    "bowtie2_cdna/benchmark/{rna_lib}_bowtie2.log"
  log:
    "bowtie2_cdna/log/{rna_lib}_bowtie2.log"
  shell:
    """
      bowtie2 \
        -x {params.regions} \
        {params.input_cmd} \
        -S {output.sam_file} &> {log}
    """
