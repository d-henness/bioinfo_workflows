#include: "mutect2_alt_bed.snakefile"
#include: "kallisto.snakefile"

configfile: "{}/ref.yaml".format(workflow.basedir)

def make_normal_set(tumors):
  normals = []
  for tumor in tumors:
    if tumors[tumor] not in normals:
      normals.append(tumors[tumor])
  return normals

def kallisto_runs(tumor):
  try:
    return f"kallisto/{config['rna_pairs'][tumor]}/kallisto/{config['rna_pairs'][tumor]}_mean_exp.tsv"
  except:
    return f"mupexi_runs/{config['pairs'][tumor]}/optitype/{config['pairs'][tumor]}_result.tsv"

rule pVACtools_all:
  input:
    "pVACtools/make_env/signal.txt"
#    expand(directory("pVACtools/{tumor}/pVACseq"), tumor = config["pairs"])


rule make_env:
  input:
  output:
    signal = "pVACtools/make_env/signal.txt",
    Wildtypes = "pVACtools/make_env/Wildtype.pm",
    Downstream = "pVACtools/make_env/Downstream.pm",
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
#for computecanada
    """
      echo "--------tensorflow---------"
      echo "--------other dependencies---------"
      pip install \
        requests \
        PyYAML \
        connexion==1.4.2 \
        biopython \
        networkx \
        simanneal \
        wget \
        mhcflurry \
        mhcnuggets \
        pysam \
        pymp-pypi \
        connexion==1.4.2 \
        py-postgresql \
        watchdog \
        vaxrank \
        flask-cors \
        bokeh==0.13.0 \
        tornado==5.0.2 \
        swagger-spec-validator==2.1.0 \
        jsonschema==2.6.0

      echo "--------pvactools--------"
      pip install pvactools[API] --no-deps
      echo "--------vcf-annotation-tools--------"
      pip install vcf-annotation-tools
      echo "--------install_vep_plugin--------"
      pvacseq install_vep_plugin pVACtools/make_env
      cp ~/projects/def-rgni0001/hg38/vep_plugins/Downstream.pm pVACtools/make_env
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
#
##this should not be needed with mutect2 vcfs
##rule add_coverage_data_p1:
##  input:
##    signal = rules.make_env.output.signal,
##    vcf_in = rules.VEP.output.vcf_out,
##  output:
##    vt_out = "pVACtools/{tumor}/add_coverage_data/{tumor}_vt.vcf",
##  conda: "./envs_dir/pVACtools.yaml"
##  log:
##    vt_log = "pVACtools/log/{tumor}_add_coverage_data_vt.log",
##  benchmark: "pVACtools/benchmark/{tumor}_add_coverage_data_p1.benchmark"
##  params:
##    ref_fasta = config["ref_fasta"],
##  threads: 1
##  resources:
##    mem_mb = 4000
##  shell:
##    """
##      vt decompose -s {input.vcf_in} -o {output.vt_out} &> {log.vt_log}
##    """
##
##rule add_coverage_data_p2:
##  input:
##    bam_in = "runs/{tumor}/ApplyBQSR/recal.bam",
##    signal = rules.make_env.output.signal,
##    vt_in = rules.add_coverage_data_p1.output.vt_out,
##    vcf_in = rules.VEP.output.vcf_out,
##  output:
##    bam_helper_snv_out = "pVACtools/{tumor}/add_coverage_data/{tumor}_bam_readcount_snv.tsv",
##    bam_helper_indel_out = "pVACtools/{tumor}/add_coverage_data/{tumor}_bam_readcount_indel.tsv",
##  conda: "./envs_dir/bam_helper.yaml"
##  log:
##    bam_helper_log = "pVACtools/log/{tumor}_add_coverage_data_bam_help.log",
##  benchmark: "pVACtools/benchmark/{tumor}_add_coverage_data_p2.benchmark"
##  params:
##    ref_fasta = config["ref_fasta"],
##  threads: 1
##  resources:
##    mem_mb = 4000
##  shell:
##    """
##      tumor_sample=$(grep 'tumor_sample' {input.vcf_in} | sed 's/^.*=//')
##      python {workflow.basedir}/scripts_dir/bam_readcount_helper.py \
##        {input.vt_in} $tumor_sample {params.ref_fasta} {input.bam_in} pVACtools/{wildcards.tumor}/add_coverage_data &> {log.bam_helper_log}
##      mv pVACtools/{wildcards.tumor}/add_coverage_data/$tumor_sample\_bam_readcount_snv.tsv pVACtools/{wildcards.tumor}/add_coverage_data/{wildcards.tumor}_bam_readcount_snv.tsv
##      mv pVACtools/{wildcards.tumor}/add_coverage_data/$tumor_sample\_bam_readcount_indel.tsv pVACtools/{wildcards.tumor}/add_coverage_data/{wildcards.tumor}_bam_readcount_indel.tsv
##    """
##
##rule add_coverage_data_p3:
##  input:
##    signal = rules.make_env.output.signal,
##    vcf_in = rules.VEP.output.vcf_out,
##    vt_in = rules.add_coverage_data_p1.output.vt_out,
##    snv_in = rules.add_coverage_data_p2.output.bam_helper_snv_out,
##    indel_in = rules.add_coverage_data_p2.output.bam_helper_indel_out,
##  output:
##    snv_out = "pVACtools/{tumor}/add_coverage_data/{tumor}_snv_annotated.vcf",
##    snv_and_indel_out = "pVACtools/{tumor}/add_coverage_data/{tumor}_snv_and_indel_annotated.vcf",
##  conda: "./envs_dir/pVACtools.yaml"
##  log:
##    snv_log = "pVACtools/log/{tumor}_snv_annotated.log",
##    indel_log = "pVACtools/log/{tumor}_snv_and_indel_annotated.log",
##  benchmark: "pVACtools/benchmark/{tumor}_add_coverage_data_p3.benchmark"
##  threads: 1
##  resources:
##    mem_mb = 4000
##  shell:
##    """
##      tumor_sample=$(grep 'tumor_sample' {input.vcf_in} | sed 's/^.*=//')
##      vcf-readcount-annotator {input.vt_in} {input.snv_in} DNA -s $tumor_sample -t snv -o {output.snv_out} &> {log.snv_log}
##      vcf-readcount-annotator {output.snv_out} {input.indel_in} DNA -s $tumor_sample -t indel -o {output.snv_and_indel_out} &> {log.indel_log}
##    """
#
#rule add_expression_data:
#  input:
#    vcf_in = rules.VEP.output.vcf_out,
#    rna = lambda wildcards: kallisto_runs(wildcards.tumor),
#  output:
#    vcf_out = "pVACtools/{tumor}/add_expression_data/{tumor}.vcf"
#  conda: "./envs_dir/pVACtools.yaml"
#  log: "pVACtools/benchmark/{tumor}_add_coverage_data_p3.benchmark"
#  benchmark: "pVACtools/benchmark/{tumor}_add_coverage_data_p3.benchmark"
#  threads: 1
#  resources:
#    mem_mb = 4000
#  shell:
#    """
#      vcf-expression-annotator {input.vcf_in} {input.rna} kallisto
#    """
#
#
##rule split_vcf_files:
##  input:
##    vcf = rules.FilterByOrientationBias.output.vcf,
##  output:
##    chr1 = "mupexi_runs/{tumor}/mupexi/{tumor}_chr1.vcf",
##    chr2 = "mupexi_runs/{tumor}/mupexi/{tumor}_chr2.vcf",
##    chr3 = "mupexi_runs/{tumor}/mupexi/{tumor}_chr3.vcf",
##    chr4 = "mupexi_runs/{tumor}/mupexi/{tumor}_chr4.vcf",
##    chr5 = "mupexi_runs/{tumor}/mupexi/{tumor}_chr5.vcf",
##    chr6 = "mupexi_runs/{tumor}/mupexi/{tumor}_chr6.vcf",
##    chr7 = "mupexi_runs/{tumor}/mupexi/{tumor}_chr7.vcf",
##    chr8 = "mupexi_runs/{tumor}/mupexi/{tumor}_chr8.vcf",
##    chr9 = "mupexi_runs/{tumor}/mupexi/{tumor}_chr9.vcf",
##    chr10 = "mupexi_runs/{tumor}/mupexi/{tumor}_chr10.vcf",
##    chr11 = "mupexi_runs/{tumor}/mupexi/{tumor}_chr11.vcf",
##    chr12 = "mupexi_runs/{tumor}/mupexi/{tumor}_chr12.vcf",
##    chr13 = "mupexi_runs/{tumor}/mupexi/{tumor}_chr13.vcf",
##    chr14 = "mupexi_runs/{tumor}/mupexi/{tumor}_chr14.vcf",
##    chr15 = "mupexi_runs/{tumor}/mupexi/{tumor}_chr15.vcf",
##    chr16 = "mupexi_runs/{tumor}/mupexi/{tumor}_chr16.vcf",
##    chr17 = "mupexi_runs/{tumor}/mupexi/{tumor}_chr17.vcf",
##    chr18 = "mupexi_runs/{tumor}/mupexi/{tumor}_chr18.vcf",
##    chr19 = "mupexi_runs/{tumor}/mupexi/{tumor}_chr19.vcf",
##    chr20 = "mupexi_runs/{tumor}/mupexi/{tumor}_chr20.vcf",
##    chr21 = "mupexi_runs/{tumor}/mupexi/{tumor}_chr21.vcf",
##    chr22 = "mupexi_runs/{tumor}/mupexi/{tumor}_chr22.vcf",
##    chrX = "mupexi_runs/{tumor}/mupexi/{tumor}_chrX.vcf",
##    chrY = "mupexi_runs/{tumor}/mupexi/{tumor}_chrY.vcf",
##  log:
##    razers3_log_2 = "mupexi_runs/log/{tumor}_split_vcf_files.log",
##  benchmark: "mupexi_runs/benchmark/{tumor}_split_vcf_files.benchmark"
##  params:
##    bioinfo_workflows_path = config["bioinfo_workflows_path"],
##    out_dir = "mupexi_runs/{tumor}/mupexi",
##  threads: 1
##  resources:
##    mem_mb = lambda wildcards, attempt: attempt * 1000,
##    time_min = lambda wildcards, attempt: attempt * 10,	# time in minutes
##  shell:
##    """
##      python3 {params.bioinfo_workflows_path}/scripts_dir/split_up_vcf_by_chr.py {input.vcf} {params.out_dir} {wildcards.tumor}
##    """
##
##
##rule pVACseq:
##  input:
##    signal = rules.make_env.output.signal,
##    vcf_in = rules.add_coverage_data_p3.output.snv_and_indel_out,
##    vcf_dir = rules.split_vcf.output.out_dir,
##  output:
##    out_dir = directory("pVACtools/{tumor}/pVACseq")
##  conda: "./envs_dir/pVACtools.yaml"
##  log: "pVACtools/log/{tumor}_pVACseq.log"
##  benchmark: "pVACtools/benchmark/{tumor}_pVACseq.benchmark"
##  params:
##    iedb = config["iedb-install-directory"],
##  threads: 1
##  resources:
##    mem_mb = 4000
##  shell:
##    """
##      i=0
##      tumor_sample=$(grep 'tumor_sample' {input.vcf_in} | sed 's/^.*=//')
##      for file in {input.vcf_dir}/*vcf; do
##        pvacseq run \
##          $file \
##          $tumor_sample \
##          HLA-A*01:01,HLA-A*32:01,HLA-B*08:01,HLA-B*14:01,HLA-C*07:01,HLA-C*08:02 \
##          MHCflurry MHCnuggetsI MHCnuggetsII NNalign NetMHC PickPocket SMM SMMPMBEC SMMalign \
##          {output.out_dir}/$i \
##          -e 8,9,10 \
##          --iedb-install-directory {params.iedb} \
##          --pass-only &> {log}_$i
##        let i=$i+1
##      done
##    """
