configfile: "{}/ref.yaml".format(workflow.basedir)

def kallisto_input(wildcards):
  if len(config['libraries'][wildcards.rna_lib]) == 1:
    return f"--single kallisto/{wildcards.rna_lib}/fastp/{wildcards.rna_lib}_1.fq.gz -l 150 -s 50"
  elif len(config['libraries'][wildcards.rna_lib]) == 2:
    return f"kallisto/{wildcards.rna_lib}/fastp/{wildcards.rna_lib}_1.fq.gz kallisto/{wildcards.rna_lib}/fastp/{wildcards.rna_lib}_2.fq.gz"

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
    expand("kallisto/{rna_lib}/kallisto/abundance.h5", rna_lib = config["libraries"])

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
    mem_mb = 4000
  threads: 1
  benchmark:
    "kallisto/benchmark/{rna_lib}_fastp.benchmark"
  log:
    "kallisto/log/{rna_lib}_fastp.log"
  shell:
    """
      fastp -i {input[0]} -I {input[1]} -o {output.fq1_out} -O {output.fq2_out} &> {log}
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
    mem_mb = 4000
  threads: 1
  benchmark:
    "kallisto/benchmark/{rna_lib}_fastp_up.benchmark"
  log:
    "kallisto/log/{rna_lib}_fastp_up.log"
  shell:
    """
      fastp -i {input} -o {output.fq1_out} --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA &> {log}
      echo "finished" > {output.signal}
    """

rule kallisto:
  input:
    signal = "kallisto/{rna_lib}/fastp/{rna_lib}_signal.txt",
  output:
    abundance = "kallisto/{rna_lib}/kallisto/abundance.h5"
  conda:
    "envs_dir/kallisto.yaml"
  params:
    index = "/data/shared/hg38/ftp.ensembl.org/pub/release-95/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.kallisto.index",
    input_cmd = kallisto_input
  resources:
    mem_mb = 16000
  threads: 1
  benchmark:
    "kallisto/benchmark/{rna_lib}_kallisto.benchmark"
  log:
    quant = "kallisto/log/{rna_lib}_kallisto_quant.log",
    h5dump = "kallisto/log/{rna_lib}_kallisto_h5dump.log",
  shell:
    """
      kallisto quant -b 500 \
        --index={params.index} \
        {params.input_cmd} \
        --output-dir=kallisto/{wildcards.rna_lib}/kallisto \
        --threads={threads} &> {log.quant}

      kallisto h5dump \
        --output-dir=kallisto/{wildcards.rna_lib}/kallisto \
        kallisto/{wildcards.rna_lib}/kallisto/abundance.h5 &> {log.h5dump}
    """