import sys
configfile: "{}/ref.yaml".format(workflow.basedir)

include: "./fastp.snakefile"

def make_normal_set(tumors):
  normals = []
  for tumor in tumors:
    if tumors[tumor] not in normals:
      normals.append(tumors[tumor])
  return normals

def make_razers3_input(wildcards, read):
    if len(config['merge_libs'][wildcards.tumor]) > 1:
        sys.exit("This workflow does not handle merge_libs with more than one lib")
    lib = config['merge_libs'][wildcards.tumor][0]
    return "fastp_dna/" + lib + "/fastp/" + lib + "_" + str(read) + ".fq.gz"
    
        


rule run_optitype:
  input:
#   expand("mupexi_runs/{lib}/optitype/{lib}_result.tsv", lib = noramls),
    expand("optitype/{tumor}/optitype/{tumor}_result.tsv", tumor = config['pairs']),

rule razers3_1:
  input:
    fq1_in = lambda wildcards: make_razers3_input(wildcards, 1),
  output:
    fq1_out = temp("optitype/{tumor}/razers3/{tumor}_1_filtered.fq"),
    bam1_txt = temp("optitype/{tumor}/razers3/all_bam_1.txt"),
    bam1_out = temp("optitype/{tumor}/razers3/{tumor}_1_filtered.bam"),
  conda:
    "./envs_dir/optitype.yaml"
  log:
    razers3_log_1 = "optitype/log/{tumor}_razers3_1.log",
  benchmark: "optitype/benchmark/{tumor}_razers3_1.benchmark"
  params:
    out_dir = "optitype/{tumor}/razers3",
  threads: 2
  resources:
    mem_mb = lambda wildcards, attempt: attempt * (16 * 1024),
    time_min = lambda wildcards, attempt: attempt * (24 * 60),	# time in minutes
  shell:
    """
      gzip -dc {input.fq1_in} | parallel -j {threads} --pipe -L4 -N1000000 "cat > {params.out_dir}/{wildcards.tumor}_1_{{#}}.fastq; razers3 -i 95 -m 1 -dr 0 -o {params.out_dir}/{wildcards.tumor}_1_{{#}}.bam $CONDA_PREFIX/bin/data/hla_reference_dna.fasta {params.out_dir}/{wildcards.tumor}_1_{{#}}.fastq; rm {params.out_dir}/{wildcards.tumor}_1_{{#}}.fastq"
      ls {params.out_dir}/{wildcards.tumor}_1_*.bam > {output.bam1_txt}
      samtools cat -b {output.bam1_txt} > {output.bam1_out}
      samtools bam2fq {output.bam1_out} > {output.fq1_out}
    """

rule razers3_2:
  input:
    fq2_in = lambda wildcards: make_razers3_input(wildcards, 2),
  output:
    fq2_out = temp("optitype/{tumor}/razers3/{tumor}_2_filtered.fq"),
    bam2_txt = temp("optitype/{tumor}/razers3/all_bam_2.txt"),
    bam2_out = temp("optitype/{tumor}/razers3/{tumor}_2_filtered.bam"),
  conda:
    "./envs_dir/optitype.yaml"
  log:
    razers3_log_2 = "optitype/log/{tumor}_razers3_2.log",
  benchmark: "optitype/benchmark/{tumor}_razers3_2.benchmark"
  params:
    out_dir = "optitype/{tumor}/razers3",
  threads: 2
  resources:
    mem_mb = lambda wildcards, attempt: attempt * (16 * 1024),
    time_min = lambda wildcards, attempt: attempt * (24 * 60),	# time in minutes
  shell:
    """
      gzip -dc {input.fq2_in} | parallel -j {threads} --pipe -L4 -N1000000 "cat > {params.out_dir}/{wildcards.tumor}_2_{{#}}.fastq; razers3 -i 95 -m 1 -dr 0 -o {params.out_dir}/{wildcards.tumor}_2_{{#}}.bam $CONDA_PREFIX/bin/data/hla_reference_dna.fasta {params.out_dir}/{wildcards.tumor}_2_{{#}}.fastq; rm {params.out_dir}/{wildcards.tumor}_2_{{#}}.fastq"
      ls {params.out_dir}/{wildcards.tumor}_2_*.bam > {output.bam2_txt}
      samtools cat -b {output.bam2_txt} > {output.bam2_out}
      samtools bam2fq {output.bam2_out} > {output.fq2_out}
    """

rule optitype:
  input:
    fq1_in = rules.razers3_1.output.fq1_out,
    fq2_in = rules.razers3_2.output.fq2_out,
  output:
    hlas = "optitype/{tumor}/optitype/{tumor}_result.tsv",
  conda:
    "./envs_dir/optitype.yaml"
  log:
    optitype_log = "optitype/log/{tumor}_optitype.log",
  benchmark: "optitype/benchmark/{tumor}_optitype.benchmark"
  params:
    out_dir = "optitype/{tumor}/optitype",
  threads: 1
  resources:
    mem_mb = lambda wildcards, attempt: attempt * 16000,
    time_min = lambda wildcards, attempt: attempt * 24 * 60,	# time in minutes
  shell:
    """
      if [[ !(-d {params.out_dir}) ]]; then
        mkdir -p {params.out_dir}
      fi
      cp $CONDA_PREFIX/bin/config.ini {params.out_dir}/config.ini
      sed -i 's/razers3.*/razers3=razers3/' {params.out_dir}/config.ini
      sed -i 's/threads.*/threads={threads}/' {params.out_dir}/config.ini

      OptiTypePipeline.py -i {input.fq1_in} {input.fq2_in} --dna -v -o {params.out_dir} -p {wildcards.tumor} -c {params.out_dir}/config.ini 2> {log.optitype_log}
    """
