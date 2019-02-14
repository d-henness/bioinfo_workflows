configfile: "{}/ref.yaml".format(workflow.basedir)

include: "mutect2.snakefile"
include: "Strelka.snakefile"
#include: "TitanCNA.snakefile"

rule phylowgs_all:
  input:
#    expand("prepro_vfcs/{tumor}/pre_pro_mutect/mutect.split.filter.vcf.gz", tumor = config["pairs"]),
#    expand("prepro_vfcs/{tumor}/pre_pro_strelka/strelka.split.filter.vcf.gz", tumor = config["pairs"]),
#    expand("prepro_vfcs/{tumor}/overlap/0002.vcf", tumor = config["pairs"])
    expand("run_phylo/{tumor}", tumor = config["pairs"])
    
    
rule pre_pro_mutect:
  input:
    vcf = lambda wildcards : "runs/" + wildcards.tumor + "/FilterByOrientationBias_" + config["pairs"][wildcards.tumor] + "/out.vcf",
  output:
    split_filter_vfc_zip = "prepro_vfcs/{tumor}/pre_pro_mutect/mutect.split.filter.vcf.gz",
    split_filter_vfc = temp("prepro_vfcs/{tumor}/pre_pro_mutect/out.split.filter.vcf"),
    split_vfc = temp("prepro_vfcs/{tumor}/pre_pro_mutect/out.split.vcf"),
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
    split_filter_vfc_zip = "prepro_vfcs/{tumor}/pre_pro_strelka/strelka.split.filter.vcf.gz",
    split_filter_vfc = temp("prepro_vfcs/{tumor}/pre_pro_strelka/out.split.filter.vcf"),
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

rule make_overlap:
  input:
    rules.pre_pro_mutect.output.split_filter_vfc_zip,
    rules.pre_pro_strelka.output.split_filter_vfc_zip
  output:
    overlap_vcf = "prepro_vfcs/{tumor}/overlap/0003.vcf"
  conda:
    "envs_dir/phylowgs.yaml"
  params:
    out_dir = "prepro_vfcs/{tumor}/overlap"
  resources:
    mem_mb = 4000
  shell:
    """
      bcftools isec -p {params.out_dir} {input}
    """

rule parse_cnvs:
  input:
    vcf = rules.make_overlap.output.overlap_vcf,
    cnv = "results/titan/hmm/optimalClusterSolution.txt"
  output:
    pre_filter = "pre_pro_cnv/{tumor}/cnvs.txt",
    post_filter = "pre_pro_cnv/{tumor}/cnvs_filter.txt"
  conda:
    "envs_dir/phylowgs.yaml"
  shell:
    '''
      PURITY=$(cat {input.cnv} | sed 's/\,[^\t]*//' | grep {wildcards.tumor} | awk '{{print $6}}')
      TITANFILE1=$(grep {wildcards.tumor} {input.cnv} | sed 's/.*results\//results\//; s/$/\.segs.txt/')
      TITANFILE2=$(echo $TITANFILE1 | sed 's/.*\///; s/^/pre_pro_cnv\/{wildcards.tumor}\//')
      if [[ !(-d pre_pro_cnv/{wildcards.tumor}) ]]; then
        mkdir -p pre_pro_cnv/{wildcards.tumor}
      fi
      cp $TITANFILE1 pre_pro_cnv/{wildcards.tumor}
      sed -i 's/\.bp\./\(bp\)/g' $TITANFILE2
      sed -i 's/Clonal_Cluster/Clonal_Frequency/' $TITANFILE2
      python2 $CONDA_PREFIX/share/phylowgs/parser/parse_cnvs.py -f titan --cnv-output {output.pre_filter} -c $PURITY $TITANFILE2
      python3 /usr/local/bioinfo_workflows/scripts_dir/filter_cnvs.py {output.pre_filter} {output.post_filter}
    '''

rule parse_cnv_and_vcf:
  input:
    cnv = rules.parse_cnvs.output.post_filter,
    vcf = rules.make_overlap.output.overlap_vcf
  output:
    out_cnv = "parse_cnv_and_vcf/{tumor}/cnv_data.txt",
    out_var = "parse_cnv_and_vcf/{tumor}/ssm_data.txt",
    out_param = "parse_cnv_and_vcf/{tumor}/params.json"
  conda:
    "envs_dir/phylowgs.yaml"
  shell:
    """
      $CONDA_PREFIX/share/phylowgs/parser/create_phylowgs_inputs.py \
        -s 5000 \
        --output-cnvs {output.out_cnv}  \
        --output-variants {output.out_var}  \
        --output-params {output.out_param}  \
        --cnvs sample1={input.cnv}  \
        --vcf-type sample1=strelka  \
        sample1={input.vcf}
    """

rule run_phylo:
  input:
    cnv = rules.parse_cnv_and_vcf.output.out_cnv,
    var = rules.parse_cnv_and_vcf.output.out_var,
    param = rules.parse_cnv_and_vcf.output.out_param,
  output:
    out_dir = directory("run_phylo/{tumor}")
  conda:
    "envs_dir/phylowgs.yaml"
  threads: 4
  benchmark:
    "benchmarks/{tumor}.run_phylo.benchmark.txt"
  shell:
    """
      python2 /usr/local/bioinfo_workflows/phylowgs/multievolve.py --num-chains 4 --ssms {input.var} --cnvs {input.cnv} -O {output.out_dir}
    """
