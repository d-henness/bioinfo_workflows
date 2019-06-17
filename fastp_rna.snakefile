configfile: "{}/ref.yaml".format(workflow.basedir)

def fastp_paired_input(wildcards):
  if len(config['rna_libraries'][wildcards.rna_lib]) == 2:
    return [config['rna_libraries'][wildcards.rna_lib][0], config['rna_libraries'][wildcards.rna_lib][1]]
  else:
    return ""

def fastp_unpaired_input(wildcards):
  if len(config['rna_libraries'][wildcards.rna_lib]) == 1:
    return f"{config['rna_libraries'][wildcards.rna_lib][0]}"
  else:
    return ""

rule fastp_all:
  input:
    expand("kallisto/{rna_lib}/fastp/{rna_lib}_signal.txt", rna_lib = config["rna_libraries"]),

rule fastp_paired:
  input:
    fastp_paired_input
  output:
    fq1_out = "kallisto/{rna_lib}/fastp/{rna_lib}_1.fq.gz",
    fq2_out = "kallisto/{rna_lib}/fastp/{rna_lib}_2.fq.gz",
    signal = "kallisto/{rna_lib}/fastp/{rna_lib}_signal.txt",
  conda:
    "envs_dir/kallisto.yaml"
  resources:
    mem_mb = lambda wildcards, attempt: attempt * (2 * 1024),
    time_min = lambda wildcards, attempt: attempt * (1 * 60),	# time in minutes
  threads: 1
  params:
    adapter_sequence = config["illumina_adapter"],
  benchmark:
    "kallisto/benchmark/{rna_lib}_fastp.benchmark"
  log:
    "kallisto/log/{rna_lib}_fastp.log"
  shell:
    """
      fastp -i {input[0]} -I {input[1]} -o {output.fq1_out} -O {output.fq2_out} --adapter_sequence={params.adapter_sequence} &> {log}
      echo "finished" > {output.signal}
    """

rule fastp_unpaired:
  input:
    fastp_unpaired_input
  output:
    fq1_out = "kallisto/{rna_lib}/fastp/{rna_lib}_1.fq.gz",
    signal = "kallisto/{rna_lib}/fastp/{rna_lib}_signal.txt",
  conda:
    "envs_dir/kallisto.yaml"
  resources:
    mem_mb = lambda wildcards, attempt: attempt * (2 * 1024),
    time_min = lambda wildcards, attempt: attempt * (1 * 60),	# time in minutes
  threads: 1
  params:
    adapter_sequence = config["illumina_adapter"],
  benchmark:
    "kallisto/benchmark/{rna_lib}_fastp_up.benchmark"
  log:
    "kallisto/log/{rna_lib}_fastp_up.log"
  shell:
    """
      fastp -i {input} -o {output.fq1_out} --adapter_sequence={params.adapter_sequence} &> {log}
      echo "finished" > {output.signal}
    """
