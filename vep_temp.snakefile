configfile: "{}/ref.yaml".format(workflow.basedir)

include: "mutect2_alt_bed.snakefile"

rule vep_all:
  input:
    expand("vep_temp/{tumor}/VEP_temp/{tumor}_VEP.vcf", tumor = config["pairs"]),

rule VEP:
  input:
    vcf_in = rules.FilterByOrientationBias.output.vcf,
  output:
    vcf_out = "vep_temp/{tumor}/VEP_temp/{tumor}_VEP.vcf",
    vcf_pass_only = "vep_temp/{tumor}/VEP_temp/{tumor}_pass_only.vcf",
  conda: "envs_dir/vep_temp.yaml"
  log: "vep_temp/log/{tumor}_VEP.log"
  benchmark: "vep_temp/benchmark/{tumor}_snp_VEP.benchmark"
  params:
    ref_fasta = config["ref_fasta"],
    vep_cache = config["vep_cache_temp"],
    vep_plugins = config["vep_plug_dir"],
    dbNSFP_config = config["dbNSFP_config"],
    condel_config = config["condel_config"],
    loftool_config = config["loftool_config"],
    bioinfo_workflows_path = config["bioinfo_workflows_path"],
  threads: 1
  resources:
    mem_mb = 6000
  shell:
    """
      bcftools filter -i "FILTER='PASS'" {input.vcf_in} -o {output.vcf_pass_only}
      vep \
        --input_file {output.vcf_pass_only} \
        --output_file {output.vcf_out} \
        --format vcf \
        --symbol \
        --offline --cache --dir_cache {params.vep_cache} \
        --quiet \
        --species "homo_sapiens" \
        --force_overwrite \
        --transcript_version &> {log}

      if [[ !(-d easy_transfer_temp) ]]; then
        mkdir easy_transfer_temp
      fi

      cp {output.vcf_out} easy_transfer_temp
    """
