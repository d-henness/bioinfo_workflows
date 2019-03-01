configfile: "{}/ref.yaml".format(workflow.basedir)

rule pVACtools_all:
  input:
    expand("pVACtools/{tumor}/add_coverage_data/{tumor}_bam_readcount_snv.tsv", tumor = config["pairs"])

rule make_env:
  input:
  output:
    signal = "pVACtools/make_env/signal.txt",
    Wildtypes = "pVACtools/make_env/Wildtype.pm"
  conda: "./envs_dir/pVACtools.yaml"
  log: "pVACtools/log/make_env.log"
  benchmark: "pVACtools/benchmark/make_env.benchmark"
  params:
    vep_version = config["vep_version"],
#    vep_plugins = config["vep_plug_dir"],
  threads: 1
  resources:
    mem_mb = 4000
# todo put in installation of VEP cashe
# todo put in installation of MHC classes
  shell:
    """
      echo "--------pvactools--------"
      pip install pvactools
      echo "--------vcf-annotation-tools--------"
      pip install vcf-annotation-tools
      echo "--------install_vep_plugin--------"
      pvacseq install_vep_plugin pVACtools/make_env
      echo "--------all done--------"
      echo "installed all dependencies" > {output.signal}
    """

rule VEP:
  input:
    vcf_in = lambda wildcards: "runs/{tumor}/FilterByOrientationBias_" + config['pairs'][wildcards.tumor] + "/out.vcf",
    signal = rules.make_env.output.signal
  output:
    vcf_out = "pVACtools/{tumor}/VEP/{tumor}_VEP.vcf",
    summary = "pVACtools/{tumor}/VEP/{tumor}_VEP.vcf_summary.html"
  conda: "./envs_dir/pVACtools.yaml"
  log: "pVACtools/log/{tumor}_VEP.log"
  benchmark: "pVACtools/benchmark/{tumor}_VEP.benchmark"
  params:
    ref_fasta = config["ref_fasta"],
    vep_cache = config["vep_cache"],
#    vep_plugins = config["vep_plug_dir"],
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
        --fasta {params.ref_fasta} \
        --offline --cache --dir_cache {params.vep_cache} \
        --plugin Downstream --plugin Wildtype --dir_plugins pVACtools/make_env \
        --pick \
        --transcript_version &> {log}
    """

rule add_coverage_data_p1:
  input:
    signal = rules.make_env.output.signal,
    vcf_in = rules.VEP.output.vcf_out,
  output:
    vt_out = "pVACtools/{tumor}/add_coverage_data/{tumor}_vt.vcf",
  conda: "./envs_dir/pVACtools.yaml"
  log:
    vt_log = "pVACtools/log/{tumor}_add_coverage_data_vt.log",
  benchmark: "pVACtools/benchmark/{tumor}_add_coverage_data.benchmark"
  params:
    ref_fasta = config["ref_fasta"],
  threads: 1
  resources:
    mem_mb = 4000
  shell:
    """
      vt decompose -s {input.vcf_in} -o {output.vt_out} &> {log.vt_log}
    """

rule add_coverage_data_p2:
  input:
    bam_in = "runs/{tumor}/ApplyBQSR/recal.bam",
    signal = rules.make_env.output.signal,
    vt_in = rules.add_coverage_data_p1.output.vt_out,
  output:
    bam_helper_snv_out = "pVACtools/{tumor}/add_coverage_data/{tumor}_bam_readcount_snv.tsv",
  conda: "./envs_dir/bam_helper.yaml"
  log:
    bam_helper_log = "pVACtools/log/{tumor}_add_coverage_data_bam_help.log",
  benchmark: "pVACtools/benchmark/{tumor}_add_coverage_data.benchmark"
  params:
    ref_fasta = config["ref_fasta"],
  threads: 1
  resources:
    mem_mb = 4000
  shell:
    """
      python /usr/local/bioinfo_workflows/scripts_dir/bam_readcount_helper.py \
        {input.vt_in} {wildcards.tumor}_0 {params.ref_fasta} {input.bam_in} pVACtools/{wildcards.tumor}/add_coverage_data &> {log.bam_helper_log}
    """
