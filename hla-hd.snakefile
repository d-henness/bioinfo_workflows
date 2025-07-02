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
    
        


rule run_hla_hd:
  input:
    expand("hla-hd{lib}/result/{lib}_final.result.txt", lib = config["libraries"]),

rule hla_hd:
  input:
    fq1 = rules.fastp_paired.output.fq1,
    fq2 = rules.fastp_paired.output.fq2,
  output:
    results = "hla-hd{lib}/result/{lib}_final.result.txt"
  conda:
    "./envs_dir/hla-hd.yaml"
  log:
    "hla-hd/log/{lib}.log",
  benchmark: "hla-hd/benchmark/{lib}.benchmark"
  threads: 4
  params:
    min_read_len = 75,
    freq_data_path = "/home/dylan/HLA-HD/hlahd.1.7.1/freq_data",
    sample_name = "{lib}",
    output_dir = "hla-hd",
  resources:
    mem_mb = lambda wildcards, attempt: attempt * (2 * 1024),
    time_min = lambda wildcards, attempt: attempt * (24 * 60),	# time in minutes
  shell:
    """
        mkdir -p {output_dir}

        hlahd.sh -t {threads} \
            -m {params.min_read_len} \
            -f {params.freq_data_path} \
            {input.fq1} {input.fq2} \
            /home/dylan/HLA-HD/hlahd.1.7.1/HLA_gene.split.txt /home/dylan/HLA-HD/hlahd.1.7.1/dictionary/ {params.sample} {params.output_dir} &> {log}
    """
