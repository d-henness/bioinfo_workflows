configfile: "{}/ref.yaml".format(workflow.basedir)


rule mixcr_all:
    input:
      expand("mixcr/{sample}/mixcr/clones.clna", sample = config["libraries"])

def mixcr_input(wildcards):
  if len(config['libraries'][wildcards.sample]) == 1:
    return "mixcr/"+ wildcards.sample + "/fastp/" + wildcards.sample + "_1.fq.gz"
  elif len(config['libraries'][wildcards.sample]) == 2:
    return "mixcr/" + wildcards.sample + "/fastp/" + wildcards.sample + "_1.fq.gz mixcr/" + wildcards.sample + "/fastp/" + wildcards.sample + "_2.fq.gz"


def fastp_paired_input(wildcards):
  if len(config['libraries'][wildcards.sample]) == 2:
    return [config['libraries'][wildcards.sample][0], config['libraries'][wildcards.sample][1]]
  else:
    return ""

def fastp_unpaired_input(wildcards):
  if len(config['libraries'][wildcards.sample]) == 1:
    return f"{config['libraries'][wildcards.sample][0]}"
  else:
    return ""

rule fastp_paired:
  input:
    fastp_paired_input
  output:
    fq1_out = "mixcr/{sample}/fastp/{sample}_1.fq.gz",
    fq2_out = "mixcr/{sample}/fastp/{sample}_2.fq.gz",
    signal = "mixcr/{sample}/fastp/{sample}_signal.txt",
  conda:
    "envs_dir/mixcr.yaml"
  resources:
    mem_mb = lambda wildcards, attempt: attempt * (6 * 1024),
    time_min = lambda wildcards, attempt: attempt * (1 * 60),	# time in minutes
  threads: 4 # master, writer, worker, reader
  params:
    adapter_sequence = config["illumina_adapter"],
    temp_dir = "/tmp/$SLURM_JOB_ID/{sample}/fastp_paired/tmp",
    fq1 = "{sample}_1.fq.gz",
    fq2 = "{sample}_2.fq.gz",
  benchmark: "mixcr/benchmark/{sample}_fastp.benchmark"
  log:
    overall = "mixcr/log/{sample}_fastp.log",
    json = "mixcr/{sample}/fastp/{sample}_log.json",
    html = "mixcr/{sample}/fastp/{sample}_log.html",
  shell:
    """
      fastp -i {input[0]} \
        -I {input[1]} \
        -o {output.fq1_out} \
        -O {output.fq2_out} \
        -w 1 \
        -j {log.json} \
        -h {log.html} \
        --trim_poly_g \
        --length_required 30 \
        --qualified_quality_phred 20 \
        --cut_right --cut_right_mean_quality 20 \
        --detect_adapter_for_pe &> {log.overall}
      echo "finished" > {output.signal}
    """


# will need to be changed to deal with io issues
#rule fastp_unpaired:
#  input:
#    fastp_unpaired_input
#  output:
#    fq1_out = "mixcr/{sample}/fastp/{sample}_1.fq.gz",
#    signal = "mixcr/{sample}/fastp/{sample}_signal.txt",
#  conda:
#    "envs_dir/mixcr.yaml"
#  resources:
#    mem_mb = lambda wildcards, attempt: attempt * 5000,
#    time_min = lambda wildcards, attempt: attempt * 24 * 60,  # time in minutes
#  threads: 1
#  benchmark:
#    "mixcr/benchmark/{sample}_fastp_up.benchmark"
#  log:
#    "mixcr/log/{sample}_fastp_up.log"
#  shell:
#    """
#      fastp -i {input} -o {output.fq1_out} --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA &> {log}
#      echo "finished" > {output.signal}
#    """

alignment_threads = 1
rule mixcr_analyze:
  input:
    signal = "mixcr/{sample}/fastp/{sample}_signal.txt",
  output:
    clones = "mixcr/{sample}/mixcr/clones.clna",
  conda: "envs_dir/mixcr.yaml",
  resources:
    mem_mb = lambda wildcards, attempt: attempt * 7 * 1024,
    time_min = lambda wildcards, attempt: attempt * 24 * 60,  # time in minutes
  params:
    input_cmd = mixcr_input,
  threads: 1 + alignment_threads
  benchmark: "mixcr/benchmark/{sample}_mixcr_align.benchmark"
  log: "mixcr/log/{sample}_mixcr_align.log",
  shell:
    """
    mixcr analyze exome-seq\
        --threads {threads}
        --species hsa\
        --assemble-longest-contigs\
        {params.input_cmd}\
        mixcr/{wildcards.sample}/

    """

#alignment_threads = 4
#rule mixcr_align:
#  input:
#    signal = "mixcr/{sample}/fastp/{sample}_signal.txt",
#  output:
#    alignments = "mixcr/{sample}/mixcr/alignments.vdjca",
#  conda: "envs_dir/mixcr.yaml",
#  resources:
#    mem_mb = lambda wildcards, attempt: attempt * 7 * 1024,
#    time_min = lambda wildcards, attempt: attempt * 24 * 60,  # time in minutes
#  params:
#    input_cmd = mixcr_input,
#  threads: 1 + alignment_threads
#  benchmark: "mixcr/benchmark/{sample}_mixcr_align.benchmark"
#  log: "mixcr/log/{sample}_mixcr_align.log",
#  shell:
#    """
#      mixcr align --species hsa -t {alignment_threads} -r {log} {params.input_cmd} {output.alignments}
#    """
#
#rule mixcr_assemble:
#  input:
#    alignments = rules.mixcr_align.output.alignments
#  output:
#    clones = "mixcr/{sample}/mixcr/clones.clna",
#  threads: 4
#  log: "mixcr/log/{sample}_mixcr_assemble.log",
#  benchmark: "mixcr/benchmark/{sample}_mixcr_assemble.benchmark"
#  resources:
#    mem_mb = lambda wildcards, attempt: attempt * 1 * 1024,
#    time_min = lambda wildcards, attempt: attempt * 24 * 60,  # time in minutes
#  conda: "envs_dir/mixcr.yaml"
#  shell:
#    "mixcr assemble -t {threads} -r {log} -a {input.alignments} {output.clones}"
#
#rule mixcr_exportclones:
#  input:
#    clones = rules.mixcr_assemble.output.clones,
#  output:
#    tra_summary = "mixcr/{sample}/mixcr/{sample}_tra_clones.txt",
#    trb_summary = "mixcr/{sample}/mixcr/{sample}_trb_clones.txt",
#    trg_summary = "mixcr/{sample}/mixcr/{sample}_trg_clones.txt",
#    all_summary = "mixcr/{sample}/mixcr/{sample}_clones.txt",
#  conda: "envs_dir/mixcr.yaml"
#  log: "mixcr/log/{sample}_mixcr_exportclones.log",
#  benchmark: "mixcr/benchmark/{sample}_mixcr_exportclones.benchmark"
#  threads: 1
#  resources:
#    mem_mb = lambda wildcards, attempt: attempt * 1 * 1024,
#    time_min = lambda wildcards, attempt: attempt * 24 * 60,  # time in minutes
#  shell:
#    """
#      mixcr exportClones -c TRA {input.clones} {output.tra_summary}
#      mixcr exportClones -c TRB {input.clones} {output.trb_summary}
#      mixcr exportClones -c TRG {input.clones} {output.trg_summary}
#      mixcr exportClones {input.clones} {output.all_summary}
#    """
