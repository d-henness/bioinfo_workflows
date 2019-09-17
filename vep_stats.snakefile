configfile: "{}/ref.yaml".format(workflow.basedir)

include: "mutect2_alt_bed.snakefile"
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
    input_list.append(f"phylowgs/{sample}/prepro_vcfs/overlap_mut_and_strel/0003.vcf")
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
    commmand_string.append(f"S{i}=phylowgs/{sample}/prepro_vcfs/overlap_mut_and_strel/0003.vcf")
  return " ".join(commmand_string)

rule vep_stats_all:
  input:
    expand("vep_stats/{tumor}/VEP/{tumor}_snp_VEP_parsed.csv", tumor = config["pairs"]),
    expand("vep_stats/{tumor}/VEP/{tumor}_indel_VEP_parsed.csv", tumor = config["pairs"]),


rule pre_pro_mutect:
  input:
    vcf = lambda wildcards : f"GATK_runs/{wildcards.tumor}/FilterByOrientationBias/{wildcards.tumor}.vcf",
  output:
    split_filter_indel_vcf = temp("vep_stats/{tumor}/prepro_vcfs/pre_pro_mutect/{tumor}.indel.filter.vcf"),
    split_filter_snp_vcf = temp("vep_stats/{tumor}/prepro_vcfs/pre_pro_mutect/{tumor}.snp.filter.vcf"),
    split_vcf = temp("vep_stats/{tumor}/prepro_vcfs/pre_pro_mutect/out.split.vcf"),
  params:
    ref_fasta = config["ref_fasta"],
  conda:
    "envs_dir/phylowgs.yaml"
  resources:
    mem_mb = 4000
  shell:
    """
      bcftools norm -f {params.ref_fasta} -m - -o {output.split_vcf} {input.vcf}
      bcftools filter -i "FILTER='PASS' && TYPE='indel'" {output.split_vcf} -o {output.split_filter_indel_vcf}
      bcftools filter -i "FILTER='PASS' && TYPE='snp'" {output.split_vcf} -o {output.split_filter_snp_vcf}
    """

rule VEP_snp:
  input:
    vcf_in = rules.pre_pro_mutect.output.split_filter_snp_vcf,
  output:
    vcf_out = "vep_stats/{tumor}/VEP/{tumor}_snp_VEP.vcf",
    vcf_out_zip = "vep_stats/{tumor}/VEP/{tumor}_snp_VEP.vcf.gz",
    summary = "vep_stats/{tumor}/VEP/{tumor}_snp_VEP.vcf_summary.html",
    parsed_output = "vep_stats/{tumor}/VEP/{tumor}_snp_VEP_parsed.csv",
  conda: "envs_dir/phylowgs.yaml"
  log: "vep_stats/log/{tumor}_snp_VEP.log"
  benchmark: "vep_stats/benchmark/{tumor}_snp_VEP.benchmark"
  params:
    ref_fasta = config["ref_fasta"],
    vep_cache = config["vep_cache"],
    vep_plugins = config["vep_plug_dir"],
    dbNSFP_config = config["dbNSFP_config"],
    condel_config = config["condel_config"],
    loftool_config = config["loftool_config"],
    bioinfo_workflows_path = config["bioinfo_workflows_path"],
  threads: 1
  resources:
    mem_mb = 4000
  shell:
    """
      vep \
        --input_file {input.vcf_in} \
        --output_file {output.vcf_out} \
        --format vcf \
        --vcf \
        --symbol \
        --terms SO \
        --tsl \
        --hgvs \
        --hgvsg \
        --fasta {params.ref_fasta} \
        --offline --cache --dir_cache {params.vep_cache} \
        --dir_plugins {params.vep_plugins} \
        --plugin Downstream \
        --plugin dbNSFP,{params.dbNSFP_config},ALL \
        --plugin Condel,{params.condel_config},b,2 \
        --plugin LoFtool,{params.loftool_config} \
        --plugin Blosum62 \
        --pick \
        --sift b \
        --polyphen b \
        --transcript_version &> {log}
      bgzip -c {output.vcf_out} > {output.vcf_out_zip}
      bcftools index -f {output.vcf_out_zip}

      python3 {params.bioinfo_workflows_path}/scripts_dir/vep_vcf_parser.py \
        -f SYMBOL Gene Consequence SIFT PolyPhen Condel LoFtool BLOSUM62 \
        -v {output.vcf_out} > {output.parsed_output}

      if [[ !(-d easy_transfer) ]]; then
        mkdir easy_transfer
      fi

      cp {output.vcf_out} {output.parsed_output} easy_transfer
    """


rule VEP_indel:
  input:
    vcf_in = rules.pre_pro_mutect.output.split_filter_indel_vcf,
  output:
    vcf_out = "vep_stats/{tumor}/VEP/{tumor}_indel_VEP.vcf",
    vcf_out_zip = "vep_stats/{tumor}/VEP/{tumor}_indel_VEP.vcf.gz",
    parsed_output = "vep_stats/{tumor}/VEP/{tumor}_indel_VEP_parsed.csv",
    summary = "vep_stats/{tumor}/VEP/{tumor}_indel_VEP.vcf_summary.html"
  conda: "envs_dir/phylowgs.yaml"
  log: "vep_stats/log/{tumor}_indel_VEP.log"
  benchmark: "vep_stats/benchmark/{tumor}_indel_VEP.benchmark"
  params:
    ref_fasta = config["ref_fasta"],
    vep_cache = config["vep_cache"],
    vep_plugins = config["vep_plug_dir"],
    dbNSFP_config = config["dbNSFP_config"],
    condel_config = config["condel_config"],
    loftool_config = config["loftool_config"],
    bioinfo_workflows_path = config["bioinfo_workflows_path"],
  threads: 1
  resources:
    mem_mb = 4000
  shell:
    """
      vep \
        --input_file {input.vcf_in} \
        --output_file {output.vcf_out} \
        --format vcf \
        --vcf \
        --symbol \
        --terms SO \
        --tsl \
        --hgvs \
        --hgvsg \
        --fasta {params.ref_fasta} \
        --offline --cache --dir_cache {params.vep_cache} \
        --dir_plugins {params.vep_plugins} \
        --plugin Downstream \
        --plugin Condel,{params.condel_config},b,2 \
        --plugin LoFtool,{params.loftool_config} \
        --plugin Blosum62 \
        --pick \
        --sift b \
        --polyphen b \
        --transcript_version &> {log}
      bgzip -c {output.vcf_out} > {output.vcf_out_zip}
      bcftools index -f {output.vcf_out_zip}

      python3 {params.bioinfo_workflows_path}/scripts_dir/vep_vcf_parser.py \
        -f SYMBOL Gene Consequence SIFT PolyPhen Condel LoFtool BLOSUM62 \
        -v {output.vcf_out} > {output.parsed_output}

      if [[ !(-d easy_transfer) ]]; then
        mkdir easy_transfer
      fi

      cp {output.vcf_out} {output.parsed_output} easy_transfer
    """

#rule VEP_overlap_snp:
#  input:
#    lambda wildcards: [f"vep_stats/{tumor}/VEP/{tumor}_snp_VEP.vcf.gz" for tumor in config["groups"][wildcards.group]]
#  output:
#    out_vcf = "vep_stats/{group}/VEP_overlap_snp/0003.vcf"
#  conda:
#    "envs_dir/phylowgs.yaml"
#  params:
#    out_dir = "vep_stats/{group}/VEP_overlap_snp"
#  resources:
#    mem_mb = 4000
#  shell:
#    """
#      if [[ !(-d "{params.out_dir}") ]]; then
#        mkdir -p "{params.out_dir}"
#      fi
#      bcftools isec -n+2 -p {params.out_dir} {input}
#    """
#
#rule VEP_overlap_indel:
#  input:
#    lambda wildcards: [f"vep_stats/{tumor}/VEP/{tumor}_indel_VEP.vcf.gz" for tumor in config["groups"][wildcards.group]]
#  output:
#    out_vcf = "vep_stats/{group}/VEP_overlap_indel/0003.vcf"
#  conda:
#    "envs_dir/phylowgs.yaml"
#  params:
#    out_dir = "vep_stats/{group}/VEP_overlap_indel"
#  resources:
#    mem_mb = 4000
#  shell:
#    """
#      if [[ !(-d "{params.out_dir}") ]]; then
#        mkdir -p "{params.out_dir}"
#      fi
#      bcftools isec -n+2 -p {params.out_dir} {input}
#    """
