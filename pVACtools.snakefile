include: "mutect2.snakefile"
configfile: "{}/ref.yaml".format(workflow.basedir)

rule pVACtools_all:
  input:
    expand(directory("pVACtools/{tumor}/split_vcf"), tumor = config["pairs"])
#    expand(directory("pVACtools/{tumor}/pVACseq"), tumor = config["pairs"])

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
      mhcflurry-downloads fetch
      echo "installed all dependencies" > {output.signal}
    """

rule VEP:
  input:
    vcf_in = rules.FilterByOrientationBias.output.vcf,
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
  benchmark: "pVACtools/benchmark/{tumor}_add_coverage_data_p1.benchmark"
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
    vcf_in = rules.VEP.output.vcf_out,
  output:
    bam_helper_snv_out = "pVACtools/{tumor}/add_coverage_data/{tumor}_bam_readcount_snv.tsv",
    bam_helper_indel_out = "pVACtools/{tumor}/add_coverage_data/{tumor}_bam_readcount_indel.tsv",
  conda: "./envs_dir/bam_helper.yaml"
  log:
    bam_helper_log = "pVACtools/log/{tumor}_add_coverage_data_bam_help.log",
  benchmark: "pVACtools/benchmark/{tumor}_add_coverage_data_p2.benchmark"
  params:
    ref_fasta = config["ref_fasta"],
  threads: 1
  resources:
    mem_mb = 4000
  shell:
    """
      tumor_sample=$(grep 'tumor_sample' {input.vcf_in} | sed 's/^.*=//')
      python {workflow.basedir}/scripts_dir/bam_readcount_helper.py \
        {input.vt_in} $tumor_sample {params.ref_fasta} {input.bam_in} pVACtools/{wildcards.tumor}/add_coverage_data &> {log.bam_helper_log}
      mv pVACtools/{wildcards.tumor}/add_coverage_data/$tumor_sample\_bam_readcount_snv.tsv pVACtools/{wildcards.tumor}/add_coverage_data/{wildcards.tumor}_bam_readcount_snv.tsv
      mv pVACtools/{wildcards.tumor}/add_coverage_data/$tumor_sample\_bam_readcount_indel.tsv pVACtools/{wildcards.tumor}/add_coverage_data/{wildcards.tumor}_bam_readcount_indel.tsv
    """

rule add_coverage_data_p3:
  input:
    signal = rules.make_env.output.signal,
    vcf_in = rules.VEP.output.vcf_out,
    vt_in = rules.add_coverage_data_p1.output.vt_out,
    snv_in = rules.add_coverage_data_p2.output.bam_helper_snv_out,
    indel_in = rules.add_coverage_data_p2.output.bam_helper_indel_out,
  output:
    snv_out = "pVACtools/{tumor}/add_coverage_data/{tumor}_snv_annotated.vcf",
    snv_and_indel_out = "pVACtools/{tumor}/add_coverage_data/{tumor}_snv_and_indel_annotated.vcf",
  conda: "./envs_dir/pVACtools.yaml"
  log:
    snv_log = "pVACtools/log/{tumor}_snv_annotated.log",
    indel_log = "pVACtools/log/{tumor}_snv_and_indel_annotated.log",
  benchmark: "pVACtools/benchmark/{tumor}_add_coverage_data_p3.benchmark"
  threads: 1
  resources:
    mem_mb = 4000
  shell:
    """
      tumor_sample=$(grep 'tumor_sample' {input.vcf_in} | sed 's/^.*=//')
      vcf-readcount-annotator {input.vt_in} {input.snv_in} DNA -s $tumor_sample -t snv -o {output.snv_out} &> {log.snv_log}
      vcf-readcount-annotator {output.snv_out} {input.indel_in} DNA -s $tumor_sample -t indel -o {output.snv_and_indel_out} &> {log.indel_log}
    """

rule split_vcf:
  input:
    signal = rules.make_env.output.signal,
    vcf_in = rules.add_coverage_data_p3.output.snv_and_indel_out
  output:
    out_dir = directory("pVACtools/{tumor}/split_vcf")
  conda: "./envs_dir/pVACtools.yaml"
  params:
  threads: 1
  resources:
    mem_mb = 16000
  shell:
    """
      python3 {workflow.basedir}/scripts_dir/split_up_vcf_by_line.py {input.vcf_in} {output.out_dir}/{wildcards.tumor} {workflow.basedir}/scripts_dir/valid_split_sites 300
    """

rule pVACseq:
  input:
    signal = rules.make_env.output.signal,
    vcf_in = rules.add_coverage_data_p3.output.snv_and_indel_out,
    vcf_dir = rules.split_vcf.output.out_dir,
  output:
    out_dir = directory("pVACtools/{tumor}/pVACseq")
  conda: "./envs_dir/pVACtools.yaml"
  log: "pVACtools/log/{tumor}_pVACseq.log"
  benchmark: "pVACtools/benchmark/{tumor}_pVACseq.benchmark"
  params:
    iedb = config["iedb-install-directory"],
  threads: 1
  resources:
    mem_mb = 4000
  shell:
    """
      i=0
      tumor_sample=$(grep 'tumor_sample' {input.vcf_in} | sed 's/^.*=//')
      for file in {input.vcf_dir}/*vcf; do
        pvacseq run \
          $file \
          $tumor_sample \
          HLA-A*01:01,HLA-A*32:01,HLA-B*08:01,HLA-B*14:01,HLA-C*07:01,HLA-C*08:02 \
          MHCflurry MHCnuggetsI MHCnuggetsII NNalign NetMHC PickPocket SMM SMMPMBEC SMMalign \
          {output.out_dir}/$i \
          -e 8,9,10 \
          --iedb-install-directory {params.iedb} \
          --pass-only &> {log}_$i
        let i=$i+1
      done
    """
