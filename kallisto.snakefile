configfile: "{}/ref.yaml".format(workflow.basedir)

include: "fastp_rna.snakefile"

# this will only work on paired fastq inputs
def kallisto_input(wildcards):
  return " ".join([f"fastp/{lib}/fastp/{lib}_1.fq.gz fastp/{lib}/fastp/{lib}_2.fq.gz" for lib in config["rna_merge_libs"][wildcards.rna_lib]])

rule kallisto_all:
  input: "kallisto/mapping.csv"

rule kallisto:
  input:
    signal = lambda wildcards: [f"fastp/{lib}/fastp/{lib}_signal.txt" for lib in config['rna_merge_libs'][wildcards.rna_lib]]
  output:
    abundance_h5 = "kallisto/{rna_lib}/kallisto/abundance.h5",
    abundance_tsv = "kallisto/{rna_lib}/kallisto/abundance.tsv",
    mean_exp = "kallisto/{rna_lib}/kallisto/{rna_lib}_mean_exp.tsv",
  conda:
    "envs_dir/kallisto.yaml"
  params:
    index = config["kallisto_index_98"],
    input_cmd = kallisto_input,
    bioinfo_workflows_path = config["bioinfo_workflows_path"],
  resources:
    mem_mb = lambda wildcards, attempt: attempt * (4 * 1024),
    time_min = lambda wildcards, attempt: attempt * (24 * 60),	# time in minutes
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
        --output-dir=kallisto/{wildcards.rna_lib}/kallisto/h5dump \
        kallisto/{wildcards.rna_lib}/kallisto/abundance.h5 &> {log.h5dump}

      python3 {params.bioinfo_workflows_path}/scripts_dir/get_mean_exp.py kallisto/{wildcards.rna_lib}/kallisto/h5dump > {output.mean_exp}
    """

rule organize:
  input:
    tsvs = expand("kallisto/{rna_lib}/kallisto/abundance.tsv", rna_lib = config["rna_merge_libs"]),
  output:
    mapping_csv = "kallisto/mapping.csv"
  resources:
    mem_mb = 1024,
    time_min = 10	# time in minutes
  threads: 1
  params:
    rna_libs = lambda wildcards: [f"\"{lib},{config['stage'][lib]},kallisto/{lib}/kallisto/\"," for lib in config['stage']],
  benchmark:
    "kallisto/benchmark/organize.benchmark"
  shell:
    """
      /usr/bin/env python3 -c '
array = [{params.rna_libs}]
print("sample,condition,path")
for lib in array:
  print(lib)
      ' > {output.mapping_csv}
    """
