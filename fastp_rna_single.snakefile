configfile: "{}/ref.yaml".format(workflow.basedir)

rule fastp_all:
  input:
    expand("fastp/{rna_lib}/fastp/{rna_lib}_signal.txt", rna_lib = config["rna_libraries"]),

rule fastp_single:
  input:
    fq1 = lambda wildcards: config["rna_libraries"][wildcards.rna_lib][0],
  output:
    fq1_out = "fastp/{rna_lib}/fastp/{rna_lib}_1.fq.gz",
    signal = "fastp/{rna_lib}/fastp/{rna_lib}_signal.txt",
  conda:
    "envs_dir/kallisto.yaml"
  resources:
    mem_mb = lambda wildcards, attempt: attempt * (6 * 1024),
    time_min = lambda wildcards, attempt: attempt * (1 * 60),	# time in minutes
  threads: 4 # master, writer, worker, reader
  params:
    adapter_sequence = config["illumina_adapter"],
    temp_dir = "/tmp/$SLURM_JOB_ID/{rna_lib}/fastp_single/tmp",
    fq1 = "{rna_lib}_1.fq.gz",
  benchmark:
    "fastp/benchmark/{rna_lib}_fastp.benchmark"
  log:
    overall = "fastp/log/{rna_lib}_fastp.log",
    json = "fastp/{rna_lib}/fastp/{rna_lib}_log.json",
    html = "fastp/{rna_lib}/fastp/{rna_lib}_log.html",
  shell:
    """
      mkdir -p {params.temp_dir}
      cp {input.fq1} {params.temp_dir}/{params.fq1}

      echo copy complete

      fastp -i {params.temp_dir}/{params.fq1} \
        -o {output.fq1_out} \
        -w 1 \
        -j {log.json} \
        -h {log.html} &> {log.overall}
      echo "finished" > {output.signal}
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
