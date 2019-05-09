configfile: "{}/ref.yaml".format(workflow.basedir)

include: "mutect2.snakefile"
include: "Strelka.snakefile"
#include: "TitanCNA.snakefile"

def parse_cnv_and_vcf_cnv_input_list(wildcards):
  input_list = []
  for i, sample in enumerate(config["phylowgs_samples"][wildcards.tumor]):
    input_list.append(f"phylowgs/{sample}/pre_pro_cnv/cnvs_filter.txt")
  return input_list

def parse_cnv_and_vcf_vcf_input_list(wildcards):
  input_list = []
  for i, sample in enumerate(config["phylowgs_samples"][wildcards.tumor]):
    input_list.append(f"phylowgs/{sample}/prepro_vfcs/overlap_mut_and_strel/0003.vcf")
  return input_list

def parse_cnv_and_vcf_cnv_command_string(wildcards):
  commmand_string = []
  for i, sample in enumerate(config["phylowgs_samples"][wildcards.tumor]):
    commmand_string.append(f"--cnvs S{i}=phylowgs/{sample}/pre_pro_cnv/cnvs_filter.txt")
  return " ".join(commmand_string)

def parse_cnv_and_vcf_vcf_command_string(wildcards):
  commmand_string = []
  for i, sample in enumerate(config["phylowgs_samples"][wildcards.tumor]):
    commmand_string.append(f"--vcf-type S{i}=strelka")
  for i, sample in enumerate(config["phylowgs_samples"][wildcards.tumor]):
    commmand_string.append(f"S{i}=phylowgs/{sample}/prepro_vfcs/overlap_mut_and_strel/0003.vcf")
  return " ".join(commmand_string)

rule phylowgs_all:
  input:
#    expand("prepro_vfcs/{tumor}/pre_pro_mutect/mutect.split.filter.vcf.gz", tumor = config["pairs"]),
#    expand("prepro_vfcs/{tumor}/pre_pro_strelka/strelka.split.filter.vcf.gz", tumor = config["pairs"]),
#    expand("prepro_vfcs/{tumor}/overlap_mut_and_strel/0003.vcf", tumor = config["pairs"]),
#   expand("prepro_vfcs/{tumor}/overlap_mut_and_var/0003.vcf", tumor = config["pairs"]),
#   expand("prepro_vfcs/{tumor}/overlap_strel_and_var/0003.vcf", tumor = config["pairs"]),
#    expand("run_phylo/{tumor}", tumor = config["pairs"])
#   expand("phylowgs/{tumor}/prepro_vfcs/overlap_mut_and_strel/0002.vcf", tumor = config["pairs"]),
#   expand("phylowgs/{tumor}/pre_pro_cnv/cnvs_filter.txt", tumor = config["pairs"]),
#    expand("phylowgs/{sample_pairs}/parse_cnv_and_vcf/cnv_data.txt", tumor = config["pairs"], sample_pairs = config["phylowgs_samples"])
#   expand("phylowgs/{sample_pairs}/run_phylo/", tumor = config["pairs"], sample_pairs = config["phylowgs_samples"]),
    expand("phylowgs/{sample_pairs}/visualize_data/{sample_pairs}_summ.json.gz", sample_pairs = config["phylowgs_samples"])

# rule pre_pro_varscan:
#   input:
#     vcf = lambda wildcards : "varscan_runs/" + wildcards.tumor + "/varscan/" + wildcards.tumor + "_snp.vcf",
#   output:
#     split_filter_vfc_zip = "phylowgs/{tumor}/prepro_vfcs/pre_pro_varscan/varscan.split.filter.vcf.gz",
#     split_filter_vfc = temp("phylowgs/{tumor}/prepro_vfcs/pre_pro_varscan/out.split.filter.vcf"),
#   params:
#     ref_fasta = config["ref_fasta"],
#   conda:
#     "envs_dir/phylowgs.yaml"
#   resources:
#     mem_mb = 4000
#   shell:
#     """
#       bcftools filter -i "TYPE='snp' && FILTER='PASS'" {input.vcf} -o {output.split_filter_vfc}
#       bgzip -c {output.split_filter_vfc} > {output.split_filter_vfc_zip}
#       bcftools index -f {output.split_filter_vfc_zip}
#     """
    
rule pre_pro_mutect:
  input:
    vcf = lambda wildcards : f"runs/{wildcards.tumor}/FilterByOrientationBias_{config['pairs'][wildcards.tumor]}/{wildcards.tumor}.vcf",
  output:
    split_filter_vfc_zip = "phylowgs/{tumor}/prepro_vfcs/pre_pro_mutect/mutect.split.filter.vcf.gz",
    split_filter_vfc = temp("phylowgs/{tumor}/prepro_vfcs/pre_pro_mutect/out.split.filter.vcf"),
    split_vfc = temp("phylowgs/{tumor}/prepro_vfcs/pre_pro_mutect/out.split.vcf"),
  params:
    ref_fasta = config["ref_fasta"],
  conda:
    "envs_dir/phylowgs.yaml"
  resources:
    mem_mb = 4000
  shell:
    """
      bcftools norm -f {params.ref_fasta} -m - -o {output.split_vfc} {input.vcf}
      bcftools filter -i "TYPE='snp' && FILTER='PASS'" {output.split_vfc} -o {output.split_filter_vfc}
      bgzip -c {output.split_filter_vfc} > {output.split_filter_vfc_zip}
      bcftools index -f {output.split_filter_vfc_zip}
    """

rule pre_pro_strelka:
  input:
    vcf = rules.Strelka_execute.output.vcfs_snvs
  output:
    split_filter_vfc_zip = "phylowgs/{tumor}/prepro_vfcs/pre_pro_strelka/strelka.split.filter.vcf.gz",
    split_filter_vfc = temp("phylowgs/{tumor}/prepro_vfcs/pre_pro_strelka/out.split.filter.vcf"),
  conda:
    "envs_dir/phylowgs.yaml"
  resources:
    mem_mb = 4000
  shell:
    """
      bcftools filter -i "TYPE='snp' && FILTER='PASS'" {input.vcf} -o {output.split_filter_vfc}
      bgzip -c {output.split_filter_vfc} > {output.split_filter_vfc_zip}
      bcftools index -f {output.split_filter_vfc_zip}
    """

rule make_overlap_mut_and_strel:
  input:
    rules.pre_pro_mutect.output.split_filter_vfc_zip,
    rules.pre_pro_strelka.output.split_filter_vfc_zip,
  output:
    overlap_vcf = "phylowgs/{tumor}/prepro_vfcs/overlap_mut_and_strel/0003.vcf"
  conda:
    "envs_dir/phylowgs.yaml"
  params:
    out_dir = "phylowgs/{tumor}/prepro_vfcs/overlap_mut_and_strel"
  resources:
    mem_mb = 4000
  shell:
    """
      bcftools isec -p {params.out_dir} {input}
    """

#rule make_overlap_mut_and_var:
#  input:
#    rules.pre_pro_mutect.output.split_filter_vfc_zip,
#    rules.pre_pro_varscan.output.split_filter_vfc_zip
#  output:
#    overlap_vcf = "prepro_vfcs/{tumor}/overlap_mut_and_var/0003.vcf"
#  conda:
#    "envs_dir/phylowgs.yaml"
#  params:
#    out_dir = "prepro_vfcs/{tumor}/overlap_mut_and_var"
#  resources:
#    mem_mb = 4000
#  shell:
#    """
#      bcftools isec -p {params.out_dir} {input}
#    """
#
#rule make_overlap_strel_and_var:
#  input:
#    rules.pre_pro_strelka.output.split_filter_vfc_zip,
#    rules.pre_pro_varscan.output.split_filter_vfc_zip
#  output:
#    overlap_vcf = "prepro_vfcs/{tumor}/overlap_strel_and_var/0003.vcf"
#  conda:
#    "envs_dir/phylowgs.yaml"
#  params:
#    out_dir = "prepro_vfcs/{tumor}/overlap_strel_and_var"
#  resources:
#    mem_mb = 4000
#  shell:
#    """
#      bcftools isec -p {params.out_dir} {input}
#    """
#

# rule VEP:
#   input:
#     vcf_in = rules.make_overlap_mut_and_strel.output.overlap_vcf,
#   output:
#     vcf_out = "prepro_vfcs/{tumor}/VEP/{tumor}_VEP.vcf",
#     summary = "prepro_vfcs/{tumor}/VEP/{tumor}_VEP.vcf_summary.html"
#   conda: "envs_dir/phylowgs.yaml"
#   log: "phylowgs/log/{tumor}_VEP.log"
#   benchmark: "phylowgs/benchmark/{tumor}_VEP.benchmark"
#   params:
#     ref_fasta = config["ref_fasta"],
#     vep_cache = config["vep_cache"],
# #    vep_plugins = config["vep_plug_dir"],
#   threads: 1
#   resources:
#     mem_mb = 4000
#   shell:
#     """
#       vep \
#         --input_file {input.vcf_in} \
#         --output_file {output.vcf_out} \
#         --format vcf \
#         --vcf \
#         --symbol \
#         --terms SO \
#         --tsl \
#         --hgvs \
#         --hgvsg \
#         --fasta {params.ref_fasta} \
#         --offline --cache --dir_cache {params.vep_cache} \
#         --plugin Downstream \
#         --pick \
#         --sift b \
#         --polyphen b \
#         --transcript_version &> {log}
#     """
rule parse_cnvs:
  input:
    cnv = "results/titan/hmm/optimalClusterSolution.txt"
  output:
    pre_filter = "phylowgs/{tumor}/pre_pro_cnv/cnvs.txt",
    post_filter = "phylowgs/{tumor}/pre_pro_cnv/cnvs_filter.txt"
  conda:
    "envs_dir/phylowgs.yaml"
  resources:
    mem_mb = 4000
  shell:
    '''
      PURITY=$(cat {input.cnv} | sed 's/\,[^\t]*//' | grep '{wildcards.tumor}_cluster' | awk '{{print $6}}')
      TITANFILE1=$(grep '{wildcards.tumor}_cluster' {input.cnv} | sed 's/.*results\//results\//; s/$/\.segs.txt/')
      TITANFILE2=$(echo $TITANFILE1 | sed 's/.*\///; s/^/phylowgs\/{wildcards.tumor}\/pre_pro_cnv\//')
      if [[ !(-d phylowgs/{wildcards.tumor}/pre_pro_cnv) ]]; then
        mkdir -p phylowgs/{wildcards.tumor}/pre_pro_cnv
      fi
      cp $TITANFILE1 phylowgs/{wildcards.tumor}/pre_pro_cnv
      sed -i 's/\.bp\./\(bp\)/g' $TITANFILE2
      sed -i 's/Clonal_Cluster/Clonal_Frequency/' $TITANFILE2
      python2 $CONDA_PREFIX/share/phylowgs/parser/parse_cnvs.py -f titan --cnv-output {output.pre_filter} -c $PURITY $TITANFILE2
      python3 {workflow.basedir}/scripts_dir/filter_cnvs.py {output.pre_filter} {output.post_filter}
    '''

rule parse_cnv_and_vcf:
  input:
    cnv_input = parse_cnv_and_vcf_cnv_input_list,
    vcf_input = parse_cnv_and_vcf_vcf_input_list,
  output:
    out_cnv = "phylowgs/{tumor}/parse_cnv_and_vcf/cnv_data.txt",
    out_var = "phylowgs/{tumor}/parse_cnv_and_vcf/ssm_data.txt",
    out_param = "phylowgs/{tumor}/parse_cnv_and_vcf/params.json"
  conda:
    "envs_dir/phylowgs.yaml"
  resources:
    mem_mb = 4000
  params:
    cnv_input = parse_cnv_and_vcf_cnv_command_string,
    vcf_input = parse_cnv_and_vcf_vcf_command_string,
  shell:
    """
      $CONDA_PREFIX/share/phylowgs/parser/create_phylowgs_inputs.py \
        --output-cnvs {output.out_cnv}  \
        --output-variants {output.out_var}  \
        --output-params {output.out_param}  \
        {params.cnv_input}  \
        {params.vcf_input}
    """

rule run_phylo:
  input:
    cnv = rules.parse_cnv_and_vcf.output.out_cnv,
    var = rules.parse_cnv_and_vcf.output.out_var,
    param = rules.parse_cnv_and_vcf.output.out_param,
  output:
    out_zip = "phylowgs/{tumor}/run_phylo/trees.zip"
  conda:
    "envs_dir/phylowgs.yaml"
  params:
    out_dir = "phylowgs/{tumor}/run_phylo/"
  threads: 4
  resources:
    mem_mb = 4000
  benchmark:
    "phylowgs/benchmarks/{tumor}.run_phylo.benchmark"
  log:
    "phylowgs/logs/{tumor}.run_phylo.log"
  shell:
    """
      python2 {workflow.basedir}/phylowgs/multievolve.py  \
        --num-chains 4  \
        --ssms {input.var}  \
        --cnvs {input.cnv}  \
        -O {params.out_dir} &> {log}
    """

rule visualize_data:
  input:
    zip_file = rules.run_phylo.output.out_zip
  output:
    tree_summary = "phylowgs/{tumor}/visualize_data/{tumor}_summ.json.gz",
    mutlist = "phylowgs/{tumor}/visualize_data/{tumor}_muts.json.gz",
    mutass = "phylowgs/{tumor}/visualize_data/{tumor}_mutass.json.gz",
  conda:
    "envs_dir/phylowgs.yaml"
  params:
    dataset_names = lambda wildcards: [sample for sample in config["phylowgs_samples"][wildcards.tumor]]
  threads: 1
  resources:
    mem_mb = 4000
  benchmark:
    "phylowgs/benchmarks/{tumor}.visualize_data.benchmark"
  log:
    "phylowgs/logs/{tumor}.visualize_data.log"
  shell:
    """
    if [[ !(-d phylowgs/{wildcards.tumor}/visualize_data) ]]; then
      mkdir -p phylowgs/{wildcards.tumor}/visualize_data
    fi

    python2 {workflow.basedir}/phylowgs/write_results.py  \
      {wildcards.tumor}  \
      {input.zip_file}  \
      {output.tree_summary}  \
      {output.mutlist}  \
      {output.mutass} &> {log}
    """
