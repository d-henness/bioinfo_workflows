configfile: "{}/ref.yaml".format(workflow.basedir)

def hisat2_input(wildcards):
  if len(config['libraries'][wildcards.rna_lib]) == 1:
    return f"-U hisat2_rna_WG/{wildcards.rna_lib}/fastp/{wildcards.rna_lib}_1.fq.gz"
  elif len(config['libraries'][wildcards.rna_lib]) == 2:
    return f"-1 hisat2_rna_WG/{wildcards.rna_lib}/fastp/{wildcards.rna_lib}_1.fq.gz -2 hisat2_rna_WG/{wildcards.rna_lib}/fastp/{wildcards.rna_lib}_2.fq.gz"

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
    expand("hisat2_rna_WG/{rna_lib}/{rna_lib}.sam", rna_lib = config["libraries"])

rule fastp_paired:
  input:
    fastp_paired_input
  output:
    fq1_out = "hisat2_rna_WG/{rna_lib}/fastp/{rna_lib}_1.fq.gz",
    fq2_out = "hisat2_rna_WG/{rna_lib}/fastp/{rna_lib}_2.fq.gz",
    signal = "hisat2_rna_WG/{rna_lib}/fastp/{rna_lib}_signal.txt",
  conda:
    "envs_dir/hisat2.yaml"
  resources:
    mem_mb = lambda wildcards, attempt: attempt * 5000,
    time_min = lambda wildcards, attempt: attempt * 24 * 60,  # time in minutes
  threads: 1
  benchmark:
    "hisat2_rna_WG/benchmark/{rna_lib}_fastp.log"
  log:
    "hisat2_rna_WG/log/{rna_lib}_fastp.log"
  shell:
    """
      fastp -i {input[0]} -I {input[1]} -o {output.fq1_out} -O {output.fq2_out} &> {log}
      echo "finished" > {output.signal}
    """

rule fastp_unpaired:
  input:
    fastp_unpaired_input
  output:
    fq1_out = "hisat2_rna_WG/{rna_lib}/fastp/{rna_lib}_1.fq.gz",
    signal = "hisat2_rna_WG/{rna_lib}/fastp/{rna_lib}_signal.txt",
  conda:
    "envs_dir/hisat2.yaml"
  resources:
    mem_mb = lambda wildcards, attempt: attempt * 5000,
    time_min = lambda wildcards, attempt: attempt * 24 * 60,  # time in minutes
  threads: 1
  benchmark:
    "hisat2_rna_WG/benchmark/{rna_lib}_fastp_up.log"
  log:
    "hisat2_rna_WG/log/{rna_lib}_fastp_up.log"
  shell:
    """
      fastp -i {input} -o {output.fq1_out} --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA &> {log}
      echo "finished" > {output.signal}
    """

rule hisat2:
  input:
    signal = "hisat2_rna_WG/{rna_lib}/fastp/{rna_lib}_signal.txt",
  output:
    sam_file = "hisat2_rna_WG/{rna_lib}/{rna_lib}.sam"
  conda:
    "envs_dir/hisat2.yaml"
  params:
    regions = config["hisat2_rna_WG"],
    input_cmd = hisat2_input
  resources:
    mem_mb = lambda wildcards, attempt: attempt * 5000,
    time_min = lambda wildcards, attempt: attempt * 24 * 60,  # time in minutes
  threads: 1
  benchmark:
    "hisat2_rna_WG/benchmark/{rna_lib}_hisat2.log"
  log:
    "hisat2_rna_WG/log/{rna_lib}_hisat2.log"
  shell:
    """
      hisat2 \
        -x {params.regions} \
        {params.input_cmd} \
        -S {output.sam_file} &> {log}
    """
