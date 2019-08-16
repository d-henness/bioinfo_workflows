configfile: "{}/ref.yaml".format(workflow.basedir)

rule download_all:
  input:
    [f"{sample}/{config['samples'][sample]}_1.fastq.gz" for sample in config['samples']]


rule download:
  output: "{sample}/{srr_num}_1.fastq.gz"
  log: "fastq_downloads/log/{sample}_{srr_num}_download.log"
  benchmark: "fastq_downloads/benchmark/{sample}_{srr_num}_download.txt"
  conda: "envs_dir/download_srr.yaml"
  threads: 6
  resources:
    mem_mb = lambda wildcards, attempt: attempt * 32 * 1024,
    time_min = lambda wildcards, attempt: attempt * 60,
  shell:
    """
      fasterq-dump -e {threads} -t /dev/shm -O {wildcards.sample} {wildcards.srr_num}
      pigz -p {threads} {wildcards.sample}/{wildcards.srr_num}*
    """
