configfile: "{}/ref.yaml".format(workflow.basedir)

def fastp_paired_input(wildcards):
  if len(config['libraries'][wildcards.library]) == 2:
    return [config['libraries'][wildcards.library][0], config['libraries'][wildcards.library][1]]
  else:
    return ""

def fastp_unpaired_input(wildcards):
  if len(config['libraries'][wildcards.library]) == 1:
    return f"{config['libraries'][wildcards.library][0]}"
  else:
    return ""

rule fastp_all:
  input:
    expand("fastp_dna/{library}/fastp/{library}_signal.txt", library = config["libraries"]),

rule fastp_paired:
  input:
    fastp_paired_input
  output:
    fq1_out = "fastp_dna/{library}/fastp/{library}_1.fq.gz",
    fq2_out = "fastp_dna/{library}/fastp/{library}_2.fq.gz",
    signal = "fastp_dna/{library}/fastp/{library}_signal.txt",
  conda:
    "envs_dir/kallisto.yaml"
  resources:
    mem_mb = lambda wildcards, attempt: attempt * (20 * 1024),
    time_min = lambda wildcards, attempt: attempt * (1 * 60),	# time in minutes
  threads: 4 # master, writer, worker, reader
  params:
    adapter_sequence = config["illumina_adapter"],
    temp_dir = "/tmp/$SLURM_JOB_ID/{library}/fastp_paired/tmp",
    fq1 = "{library}_1.fq.gz",
    fq2 = "{library}_2.fq.gz",
  benchmark:
    "fastp_dna/benchmark/{library}_fastp.benchmark"
  log:
    overall = "fastp_dna/log/{library}_fastp.log",
    json = "fastp_dna/{library}/fastp/{library}_log.json",
    html = "fastp_dna/{library}/fastp/{library}_log.html",
  shell:
    """

      fastp -i {input[0]} \
        -I {input[1]} \
        -o {output.fq1_out} \
        -O {output.fq2_out} \
        -w 1 \
        -j {log.json} \
        -h {log.html} \
        -D \
        --dup_calc_accuracy 5 \
        --trim_poly_g \
        --length_required 30 \
        --qualified_quality_phred 20 \
        --cut_right --cut_right_mean_quality 20 \
        --detect_adapter_for_pe &> {log.overall}
      echo "finished" > {output.signal}

      fastqc {output.fq1_out} {output.fq2_out} -o fastp_dna/{wildcards.library}/fastp/
    """

# will need to be changed to deal with io issues
#rule fastp_unpaired:
#  input:
#    fastp_unpaired_input
#  output:
#    fq1_out = "fastp/{rna_lib}/fastp/{rna_lib}_1.fq.gz",
#    signal = "fastp/{rna_lib}/fastp/{rna_lib}_signal.txt",
#  conda:
#    "envs_dir/kallisto.yaml"
#  resources:
#    mem_mb = lambda wildcards, attempt: attempt * (6 * 1024),
#    time_min = lambda wildcards, attempt: attempt * (1 * 60),	# time in minutes
#  threads: 1
#  params:
#    adapter_sequence = config["illumina_adapter"],
#  benchmark:
#    "fastp/benchmark/{rna_lib}_fastp_up.benchmark"
#  log:
#    "fastp/log/{rna_lib}_fastp_up.log"
#  shell:
#    """
#      fastp -i {input} -o {output.fq1_out} --adapter_sequence={params.adapter_sequence} &> {log}
#      echo "finished" > {output.signal}
#    """
