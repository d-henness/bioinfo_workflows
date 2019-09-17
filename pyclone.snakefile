configfile: "{}/ref.yaml".format(workflow.basedir)

#include: "TitanCNA_alt.snakefile"
include: "Strelka.snakefile"
include: "mutect2_alt_bed.snakefile"


plot_clusters_types = [
  "density",
  "parallel_coordinates",
  "scatter"
]

build_table_types = [
  'cluster',
  'loci',
]

plot_loci_types = [
  'density',
  'parallel_coordinates',
  'scatter',
  'similarity_matrix',
  'vaf_parallel_coordinates',
  'vaf_scatter'
]

rule pyclone_all:
#  input: expand("pyclone/{tumor}/plot_loci_{plot_loci_type}.png", tumor = config['pairs'], plot_loci_type = plot_loci_types)
  input:
    expand("pyclone/{tumor}/{tumor}_plot_clusters_{plot_type}.png", tumor = config['pairs'], plot_type = plot_clusters_types),
    expand("pyclone/{tumor}/{tumor}_build_tables_{table_type}.tsv", tumor = config['pairs'], table_type = build_table_types)

rule pre_pro_mutect:
  input:
    vcf = rules.FilterByOrientationBias.output.vcf #lambda wildcards : f"GATK_runs/{wildcards.tumor}/FilterByOrientationBias_{config['pairs'][wildcards.tumor]}/{wildcards.tumor}.vcf",
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

rule parse_inputs:
  input:
    vcf_snp = rules.make_overlap_mut_and_strel.output.overlap_vcf_snp_mutect2,
    vcf_indel = rules.make_overlap_mut_and_strel.output.overlap_vcf_indel_mutect2,
    cna_results = "results/titan/hmm/optimalClusterSolution.txt",
  output:
    tsv_out = "pyclone/{tumor}/{tumor}_parse_input.tsv",
    yaml_out = "pyclone/{tumor}/{tumor}_parse_input.yaml",
    yaml_run = "pyclone/{tumor}/{tumor}_pyclone_input.yaml",
  params:
    bioinfo_workflows_path = config["bioinfo_workflows_path"],
  conda: "envs_dir/pyclone.yaml"
  log: "pyclone/log/{tumor}_parse_inputs.log"
  benchmark: "pyclone/benchmarks/{tumor}_parse_inputs.benchmark"
  resources:
    mem_mb = 4000
  shell:
    """
      purity=$(grep '\s{wildcards.tumor}\s' {input.cna_results} | sed -E 's/\s+,/,/' | cut -f 6)

      titan_seg_file=$(grep '\s{wildcards.tumor}\s' {input.cna_results} | sed -E 's/.*\s([^ ]+$)/\\1/').segs.txt

      if [[ !(-d pyclone/{wildcards.tumor}) ]]; then
        mkdir -p pyclone/{wildcards.tumor}
      fi

      python3 {params.bioinfo_workflows_path}/scripts_dir/make_pyclone_input.py {input.vcf_snp} {input.vcf_indel} $titan_seg_file > {output.tsv_out}
      cp {output.tsv_out} {output.tsv_out}-$(date +%Y%m%d%H%M%S)
      PyClone build_mutations_file --in_file {output.tsv_out} --out_file {output.yaml_out} &> {log}
      python3 {params.bioinfo_workflows_path}/scripts_dir/make_pyclone_yaml.py $PWD/pyclone/{wildcards.tumor}/ {wildcards.tumor} $purity > {output.yaml_run}
    """

rule pyclone_run_analysis:
  input:
    yaml_run = rules.parse_inputs.output.yaml_run,
  output:
    out_dir = directory("pyclone/{tumor}/trace"),
  params:
    bioinfo_workflows_path = config["bioinfo_workflows_path"],
  conda: "envs_dir/pyclone.yaml"
  log: "pyclone/log/{tumor}_run_analysis.log"
  benchmark: "pyclone/benchmarks/{tumor}_run_analysis.benchmark"
  resources:
    mem_mb = lambda wildcards, attempt: attempt * 5000,
    time_min = lambda wildcards, attempt: attempt * 24 * 60,	# time in minutes
  shell:
    """
    PyClone run_analysis --config_file {input.yaml_run} &> {log}
    """

rule pyclone_plot_loci:
  input:
    yaml_run = rules.parse_inputs.output.yaml_run,
    trace_dir = rules.pyclone_run_analysis.output.out_dir,
  output:
    plot = "pyclone/{tumor}/{tumor}_plot_loci_{plot_loci_type}.png",
  params:
    bioinfo_workflows_path = config["bioinfo_workflows_path"],
  conda: "envs_dir/pyclone.yaml"
  log: "pyclone/log/{tumor}_plot_loci_{plot_loci_type}.log"
  benchmark: "pyclone/benchmarks/{tumor}_plot_loci_{plot_loci_type}.benchmark"
  resources:
    mem_mb = lambda wildcards, attempt: attempt * 5000,
    time_min = lambda wildcards, attempt: attempt * 24 * 60,	# time in minutes
  shell:
    """
    PyClone plot_loci \
      --config_file {input.yaml_run} \
      --plot_file {output.plot} \
      --plot_type {wildcards.plot_loci_type} \
      --burnin 1000 \
      --thin 1 &> {log}
    """

rule pyclone_plot_clusters:
  input:
    yaml_run = rules.parse_inputs.output.yaml_run,
    trace_dir = rules.pyclone_run_analysis.output.out_dir,
  output:
    plot = "pyclone/{tumor}/{tumor}_plot_clusters_{plot_type}.png",
  params:
    bioinfo_workflows_path = config["bioinfo_workflows_path"],
  conda: "envs_dir/pyclone.yaml"
  log: "pyclone/log/{tumor}_plot_clusters_{plot_type}.log"
  benchmark: "pyclone/benchmarks/{tumor}_plot_clusters_{plot_type}.benchmark"
  resources:
    mem_mb = lambda wildcards, attempt: attempt * 5000,
    time_min = lambda wildcards, attempt: attempt * 24 * 60,	# time in minutes
  shell:
    """
    PyClone plot_clusters \
      --config_file {input.yaml_run} \
      --plot_file {output.plot} \
      --plot_type {wildcards.plot_type} \
      --burnin 1000 \
      --thin 1 &> {log}
    """

rule pyclone_build_tables:
  input:
    yaml_run = rules.parse_inputs.output.yaml_run,
    trace_dir = rules.pyclone_run_analysis.output.out_dir,
  output:
    table = "pyclone/{tumor}/{tumor}_build_tables_{table_type}.tsv",
  params:
    bioinfo_workflows_path = config["bioinfo_workflows_path"],
  conda: "envs_dir/pyclone.yaml"
  log: "pyclone/log/{tumor}_build_tables_{table_type}.log"
  benchmark: "pyclone/benchmarks/{tumor}_build_tables_{table_type}.benchmark"
  resources:
    mem_mb = lambda wildcards, attempt: attempt * 5000,
    time_min = lambda wildcards, attempt: attempt * 24 * 60,	# time in minutes
  shell:
    """
    PyClone build_table \
      --config_file {input.yaml_run} \
      --out_file {output.table} \
      --table_type {wildcards.table_type} \
      --burnin 1000 \
      --thin 1 &> {log}
    """
