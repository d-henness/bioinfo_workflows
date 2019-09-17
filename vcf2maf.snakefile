configfile: "{}/ref.yaml".format(workflow.basedir)

include: "mutect2_alt_bed.snakefile"

tumors = config["pairs"]
normals = [config["pairs"][tumor] for tumor in tumors]

rule vcf2maf_all:
  input:
    expand("maf/{tumor}/{tumor}.vep.maf", tumor = tumors)

rule vcf2maf:
  input:
    vcf = rules.FilterByOrientationBias.output.vcf,
  output:
    maf = "maf/{tumor}/{tumor}.vep.maf",
    tmp_dir = directory("maf/{tumor}/vep"),
  conda:
    "envs_dir/vcf2maf.yaml"
  threads: 4
  resources:
    mem_mb = lambda wildcards, attempt: attempt * 10000,
    time_min = lambda wildcards, attempt: attempt * 30,	# time in minutes
  log: "maf/logs/{tumor}_vcf2maf.log"
  params:
    vep_cache = config["vep_cache"],
    ref_fasta = config["ref_fasta"],
#    ref_fai = config["ref_fasta_index"],
#    ref_dict = config["ref_dict"],
#    gnomad = config["gnomad"],
#    interval = config["alt_bed"],
#    variants_con = config["variants_for_contamination"],
#    variants_con_index = config["variants_for_contamination_index"],
#    exclude_list = ''
  benchmark: "maf/benchmarks/{tumor}_vcf2maf_benchmark.txt"
  shell:
    """
      if [[ !(-d {output.tmp_dir}) ]]; then
        mkdir -p {output.tmp_dir}
      fi

      vcf2maf.pl \
        --input-vcf {input.vcf} \
        --output-maf {output.maf} \
        --vep-path $CONDA_PREFIX/bin \
        --vep-data {params.vep_cache} \
        --ref-fasta {params.ref_fasta} \
        --tumor-id {wildcards.tumor} \
        --vcf-tumor-id TUMOR \
        --vep-forks {threads} \
        --tmp-dir {output.tmp_dir} \
        --ncbi-build GRCh38 &> {log}
    """
