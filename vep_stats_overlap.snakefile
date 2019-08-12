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
    expand("vep_stats/{tumor}/VEP/{tumor}_snp_VEP.vcf", tumor = config["pairs"]),
    expand("vep_stats/{tumor}/VEP/{tumor}_indel_VEP.vcf", tumor = config["pairs"]),
    expand("vep_stats/{group}/VEP_overlap_snp/0003.vcf", group = config["groups"]),
    expand("vep_stats/{group}/VEP_overlap_indel/0003.vcf", group = config["groups"]),


rule pre_pro_mutect:
  input:
    vcf = lambda wildcards : f"runs/{wildcards.tumor}/FilterByOrientationBias_{config['pairs'][wildcards.tumor]}/{wildcards.tumor}.vcf",
  output:
    split_filter_vcf_zip = "vep_stats/{tumor}/prepro_vcfs/pre_pro_mutect/mutect.split.filter.vcf.gz",
    split_filter_vcf = temp("vep_stats/{tumor}/prepro_vcfs/pre_pro_mutect/out.split.filter.vcf"),
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
      bcftools filter -i "FILTER='PASS'" {output.split_vcf} -o {output.split_filter_vcf}
      bgzip -c {output.split_filter_vcf} > {output.split_filter_vcf_zip}
      bcftools index -f {output.split_filter_vcf_zip}
    """

rule pre_pro_strelka:
  input:
    vcf_snp = rules.Strelka_execute.output.vcfs_snvs,
    vcf_indel = rules.Strelka_execute.output.vcfs_indels
  output:
    split_filter_vcf_zip = "vep_stats/{tumor}/prepro_vcfs/pre_pro_strelka/strelka.split.filter.vcf.gz",
    split_filter_vcf = temp("vep_stats/{tumor}/prepro_vcfs/pre_pro_strelka/out.split.filter.vcf"),
    norm_vcf_indel = temp("vep_stats/{tumor}/prepro_vcfs/pre_pro_strelka/norm_vcf_indel.vcf"),
    filter_vcf_indel = temp("vep_stats/{tumor}/prepro_vcfs/pre_pro_strelka/filter_vcf_indel.vcf"),
    filter_vcf_indel_zip = "vep_stats/{tumor}/prepro_vcfs/pre_pro_strelka/filter_vcf_indel.vcf.gz",
  params:
    ref_fasta = config["ref_fasta"],
  conda:
    "envs_dir/phylowgs.yaml"
  resources:
    mem_mb = 4000
  shell:
    """
      bcftools norm -f {params.ref_fasta} -m - -o {output.norm_vcf_indel} {input.vcf_indel}
      bcftools filter -i "FILTER='PASS'" {output.norm_vcf_indel} -o {output.filter_vcf_indel}
      bgzip -c {output.filter_vcf_indel} > {output.filter_vcf_indel_zip}
      bcftools index -f {output.filter_vcf_indel_zip}

      bcftools filter -i "FILTER='PASS'" {input.vcf_snp} -o {output.split_filter_vcf}
      bgzip -c {output.split_filter_vcf} > {output.split_filter_vcf_zip}
      bcftools index -f {output.split_filter_vcf_zip}
    """

rule make_overlap_mut_and_strel:
  input:
    mut = rules.pre_pro_mutect.output.split_filter_vcf_zip,
    strel_snp = rules.pre_pro_strelka.output.split_filter_vcf_zip,
    strel_indel = rules.pre_pro_strelka.output.filter_vcf_indel_zip,
  output:
    overlap_vcf_snp = "vep_stats/{tumor}/prepro_vcfs/overlap_mut_and_strel_snp/0003.vcf",
    overlap_vcf_indel = "vep_stats/{tumor}/prepro_vcfs/overlap_mut_and_strel_indel/0003.vcf",
    overlap_vcf_snp_mutect2 = "vep_stats/{tumor}/prepro_vcfs/overlap_mut_and_strel_snp/0002.vcf",
    overlap_vcf_indel_mutect2 = "vep_stats/{tumor}/prepro_vcfs/overlap_mut_and_strel_indel/0002.vcf",
  conda:
    "envs_dir/phylowgs.yaml"
  params:
    out_dir_snp = "vep_stats/{tumor}/prepro_vcfs/overlap_mut_and_strel_snp",
    out_dir_indel = "vep_stats/{tumor}/prepro_vcfs/overlap_mut_and_strel_indel",
  resources:
    mem_mb = 4000
  shell:
    """
      bcftools isec -p {params.out_dir_snp} {input.mut} {input.strel_snp}
      bcftools isec -p {params.out_dir_indel} {input.mut} {input.strel_indel}
    """

rule VEP_snp:
  input:
    vcf_in = rules.make_overlap_mut_and_strel.output.overlap_vcf_snp,
  output:
    vcf_out = "vep_stats/{tumor}/VEP/{tumor}_snp_VEP.vcf",
    vcf_out_zip = "vep_stats/{tumor}/VEP/{tumor}_snp_VEP.vcf.gz",
    summary = "vep_stats/{tumor}/VEP/{tumor}_snp_VEP.vcf_summary.html"
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
    """


rule VEP_indel:
  input:
    vcf_in = rules.make_overlap_mut_and_strel.output.overlap_vcf_indel,
  output:
    vcf_out = "vep_stats/{tumor}/VEP/{tumor}_indel_VEP.vcf",
    vcf_out_zip = "vep_stats/{tumor}/VEP/{tumor}_indel_VEP.vcf.gz",
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
    """

rule VEP_overlap_snp:
  input:
    lambda wildcards: [f"vep_stats/{tumor}/VEP/{tumor}_snp_VEP.vcf.gz" for tumor in config["groups"][wildcards.group]]
  output:
    out_vcf = "vep_stats/{group}/VEP_overlap_snp/0003.vcf"
  conda:
    "envs_dir/phylowgs.yaml"
  params:
    out_dir = "vep_stats/{group}/VEP_overlap_snp"
  resources:
    mem_mb = 4000
  shell:
    """
      if [[ !(-d "{params.out_dir}") ]]; then
        mkdir -p "{params.out_dir}"
      fi
      bcftools isec -n+2 -p {params.out_dir} {input}
    """

rule VEP_overlap_indel:
  input:
    lambda wildcards: [f"vep_stats/{tumor}/VEP/{tumor}_indel_VEP.vcf.gz" for tumor in config["groups"][wildcards.group]]
  output:
    out_vcf = "vep_stats/{group}/VEP_overlap_indel/0003.vcf"
  conda:
    "envs_dir/phylowgs.yaml"
  params:
    out_dir = "vep_stats/{group}/VEP_overlap_indel"
  resources:
    mem_mb = 4000
  shell:
    """
      if [[ !(-d "{params.out_dir}") ]]; then
        mkdir -p "{params.out_dir}"
      fi
      bcftools isec -n+2 -p {params.out_dir} {input}
    """
