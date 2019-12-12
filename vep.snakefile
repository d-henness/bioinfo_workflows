configfile: "{}/ref.yaml".format(workflow.basedir)

rule vep_all:
  input:
    expand("vep/{tumor}/VEP/{tumor}_VEP.vcf", tumor = config["tumors"]),

rule VEP:
  input:
    vcf_in = lambda wildcards: config["tumors"][wildcards.tumor],
  output:
    vcf_pass_only = "vep/{tumor}/VEP/{tumor}_pass_only.vcf",
    vcf_out = "vep/{tumor}/VEP/{tumor}_VEP.vcf",
    vcf_out_zip = "vep/{tumor}/VEP/{tumor}_VEP.vcf.gz",
    summary = "vep/{tumor}/VEP/{tumor}_VEP.vcf_summary.html",
    parsed_output = "vep/{tumor}/VEP/{tumor}_VEP_parsed.csv",
  conda: "envs_dir/phylowgs.yaml"
  log: "vep/log/{tumor}_VEP.log"
  benchmark: "vep/benchmark/{tumor}_snp_VEP.benchmark"
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
      bcftools filter -i "FILTER='PASS'" {input.vcf_in} -o {output.vcf_pass_only}
      vep \
        --input_file {output.vcf_pass_only} \
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
